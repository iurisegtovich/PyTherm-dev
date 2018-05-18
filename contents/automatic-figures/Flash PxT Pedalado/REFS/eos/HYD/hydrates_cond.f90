!hydrates_2.f90	baseado nos parÂmetros de Munck et al. 1988, expressões de Holder mais extensão até a condição de referencoia de gas idela a partir de eos cubica para agua pura 

module hydrates_cond_mod

	use hydrates_mod! which does use phase_model_mod !which does use comp_list_mod !which does use pureprop_mod
	use ceos_pure_mod
	use pengrob_mod
	use ice_mod
	implicit none
	type, extends(hydrate_model_c) :: hydrate_cond_c
!as variáveis ou rotinas que eu não quiser usar eu deixo tipo deferred e não chamo, refaço o que precisar, e uso o que eu quiser
		!declarations
		class(ceos_pure_c), pointer :: EOSslave
		class(ice_sublim_c), pointer :: ICEslave
		real(8) :: dmi_EL_LW_I, fugrefW
		real(8) :: delmu0, delcp0(2), q(2), delh0(2), delv0(2)
!		real(8) :: T0 = 273.15d0, P0 = 1.d5 !0.d0 !condições de referência
		real(8) :: T0 = 273.17944237908694, P0 = 612.61924977751960 !condições de referência: ponto triplo da água de acordo com as expressões e parâmetros implementadas para eos e gelo

		real(8), allocatable :: A_pSWP(:,:) !nguest,ncav
		real(8), allocatable :: B_pSWP(:,:)!(epsilon/kb) em [Kelvin] !nguest,ncav

			!REVER referência, acho q deveria ser 273.17K/0bar, ou 273.16/600bar, ou 273.15/1bar
			!pois são nessas condições que fug Lw = fug Iw e logo o delta_mi_(T0,P0) seria um só
				!usar 273.15 pois é em relação a essa que os deltaH e deltaCp estão ajustados
		contains
!		// hydrate_cond_c specific methods
		procedure :: calcClang => calcClang_pSWP
		procedure :: calcVEL => calcVEL_88
		procedure :: calc_dmi_EL_refW
		procedure :: calc_fugrefW
		procedure :: calcFugWH => calcFugWH_88
!		// virtual methods
!		procedure :: calc_phi ! sem override
		procedure :: debug_phase_model => debug_hydcond
	end type hydrate_cond_c
!	
	interface hydrate_cond_i
		module procedure init_hyd_cond_fnum_defslavew !baseado em numero de arquivo, modelo default para vapor saturado e gelo
	end interface
		
	contains
!	
	subroutine calcClang_pSWP(z,T,P)
		class(hydrate_cond_c) :: z
		real(8) :: T, P
		!IMPLEMENT
		integer :: i, j
		
		do i = 1, z%nguest
			do j = 1, 2 !ncav
				z%Clang(i,j) = (z%A_pSWP(i,j)/T)*dexp(z%B_pSWP(i,j)/T)/1.d5
				!print*, z%Clang(i,j) ![Pa]
			end do !j = 1, 2 !ncav
		end do !i = 1, nguest
		
	end subroutine calcClang_pSWP
	
	subroutine calc_dmi_EL_refW(z, T, P)
		class(hydrate_cond_c) :: z
		real(8) :: T, P
		integer :: i
		real(8) :: T0, P0, invT0
		real(8) :: h0, v0, cp0, q, dmu0
		
		T0 = z%T0			!REVER referência, acho q deveria ser 173.17K/0bar, ou 273.16/600bar, ou 273.15/1bar
		invT0 = 1.d0/T0		!pois são nessas condições que fug Lw = fug Iw e logo o delta_mi_(T0,P0) seria um só
		P0 = z%P0
		
		if (T <= T0) then
!			dmu0 = 0.d0 !Munck !do jeito dele aparece um descontinuidade no gráfico, fica absurdo para T<T0, aparentemente o artigo explicou errado mas deve ter programado certo
			dmu0 = z%delmu0 !eu
			cp0 = z%delcp0(1)
			q = z%q(1)
			h0 = z%delh0(1)
			v0 = z%delv0(1)
		else
			dmu0 = z%delmu0
			cp0 = z%delcp0(2)
			q = z%q(2)
			h0 = z%delh0(2)
			v0 = z%delv0(2)
		end if
		
		z%dmi_EL_LW_I = dmu0*T*invT0 + & 
(1.d0-T*invT0)*(h0-cp0*T0+(q/2.d0)*(T0**2.0d0)) - &
T*(cp0-q*T0)*DLOG(T*invT0) - (q/2.d0)*T*(T-T0)+ v0*(P-P0) !eu

!		P0 = 0.d0 
!		z%dmi_EL_LW_I = dmu0*T*invT0 + & 
!(1.d0-T*invT0)*(h0-cp0*T0+(q/2.d0)*(T0**2.0d0)) - &
!T*(cp0-q*T0)*DLOG(T*invT0) - (q/2.d0)*T*(T-T0)+ v0*(P-P0)*((2.d0*T)/(T+T0)) !Munck !não muda quase nada, só um pouco em T alta
		
	end subroutine calc_dmi_EL_refW
	
	subroutine calc_fugrefW(z,T,P)
		class(hydrate_cond_c) :: z
		real(8) :: T, P
		
		if (T <= z%T0) then							!porque a fugacidade do gelo é calculada em relação a fugacidade de uma fase de água pura nas condições em questão
														!posso escolher qualquer temperatura desde que use a mesa pra selecionar o dH e o dV e o dCp
														!mas se T0 e P0 não pertencerem a curva de eq gelo-lw então mi(0) não vai ser um único
														!até poderia usar sempre em relação à água liquida nesse programa aqui, pois ele sozinho consegue prever a diferença de fugacidade entre o Lw e o gelo pra prever o equilibrio na zona de hyd-gelo
			z%ICEslave%x(:) = 0.d0
			z%ICEslave%x(1) = 1.d0
			call z%ICEslave%calc_phi(T,P)
			z%fugrefW = z%ICEslave%fug_w_ice
		else
			z%EOSslave%x(:) = 0.d0
			z%EOSslave%x(1) = 1.d0
			call z%EOSslave%calc_phi(T,P)
			z%fugrefW = z%EOSslave%phi(1)*P
		end if
	end subroutine calc_fugrefW
	
	subroutine debug_hydcond(z)
		class(hydrate_cond_c) :: z
		integer :: i ,j
		call z%phase_model_c%debug_phase_model
		print*,
		print*, 'debug_hyd_cond ', z%class_name
		print*, 'z%delmu0', z%delmu0
		print*, 'z%delcp0(:)', z%delcp0(:)
		print*, 'z%q(:)', z%q(:)
		print*, 'z%delh0(:)', z%delh0(:)
		print*, 'z%delv0(:)', z%delv0(:)
		print*, 'z%T0', z%T0, 'z%P0', z%P0 
		do j = 1, 2
			print*, 'z%A_pSWP(:,',j,')', z%A_pSWP(:,j)
		end do !j = 1, ncav
		do j = 1, 2
			print*, 'z%B_pSWP(:,',j,')', z%A_pSWP(:,j)
		end do !j = 1, ncav
		
		call z%EOSslave%debug_phase_model
		call z%ICEslave%debug_phase_model
	end subroutine debug_hydcond
	
	function init_hyd_cond_fnum_defslavew(fnum)
		!declarations
		type(hydrate_cond_c) :: z, init_hyd_cond_fnum_defslavew
		integer :: i, fnum
		!preamble
		z%hydrate_model_c = hydrate_model_i(fnum)
		z%class_name = 'hyd_cond'
		!implementation
		allocate(z%A_pSWP(z%nguest,2), z%B_pSWP(z%nguest,2))
		
!		open(unit=10,file=trim(phase_file_name),status='old',action='read')
!		ler parâmetros
!		close(10)
!		print*, z%phase_cond, 'phase_cond', __FILE__, __LINE__
		
		select case(trim(z%phase_cond))
		case('S1')
			
			do i = 1, z%nguest
				if( comp_list(z%guestindex(i))%comp%name_ == 'methane' ) then
					z%A_pSWP(i,1) = .7228d-3
					z%B_pSWP(i,1) = 3187.d0
					z%A_pSWP(i,2) = 23.35d-3
					z%B_pSWP(i,2) = 2653.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'ethane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = 3.039d-3
					z%B_pSWP(i,2) = 3861.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'propane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = .0d-3
					z%B_pSWP(i,2) = 0.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'i_butane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = .0d-3
					z%B_pSWP(i,2) = 0.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'n_butane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = .0d-3
					z%B_pSWP(i,2) = 0.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'carbon_dioxide' ) then
					z%A_pSWP(i,1) = .2474d-3
					z%B_pSWP(i,1) = 3410.d0
					z%A_pSWP(i,2) = 42.46d-3
					z%B_pSWP(i,2) = 2813.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'THF' ) then
					print*, 'need Clang par for THF @ ', __FILE__, __LINE__; stop
				else
					print*, 'unnexpected comp in hydrate phase model'; stop
				end if
			end do !i = 1, z%nguest
		
			z%delmu0 = 1264.d0 !1264.4136d0 !PP_1972 //	1120 !Holder_1988(apudJPH1985)
			
			!(1):gelo (2):liq
			
			z%delh0(1) = 1151.d0 !1151.37d0 !PP_1972 // 1714 !Holder_1988(apudJPH1985)
			z%delh0(2) = -4858.d0 !z%delh0(1)-6011.d0 !-4862.13084d0 !subtract 6011 !PP_1972 e !Holder_1988
			
			z%delcp0(1) = 0.d0 !-38.141748d0 !PP_1972 //	-37.32 !Holder1980 (T>T0) //	-34.583 !Holder_1988(apudJPH1985) (T>T0)
			z%delcp0(2) = -39.16d0 !							0.565 !Holder1980 (T<T0) //		3.315 !Holder_1988(apudJPH1985) (T<T0)
			z%q(1) = 0.d0 !0.14067648d0 !PP_1972 //	0.179 !Holder1980 (T>T0) //	0.189 !Holder_1988(apudJPH1985) (T>T0)
			z%q(2) = 0.d0 !0.0121d0 !tabela do artigo de 1985
			
			z%delv0(1) = 3.d-6 !3.0d-6 !PP_1972 // 2.9959d-6 !Holder_1988(apudJPH1985)
			z%delv0(2) = z%delv0(1) + 1.6d-6 !4.6d-6 !add 1.6 !Holder_1988(apudJPH1985)
			
		case('S2')

			do i = 1, z%nguest
				if( comp_list(z%guestindex(i))%comp%name_ == 'methane' ) then
					z%A_pSWP(i,1) = .2207d-3
					z%B_pSWP(i,1) = 3453.d0
					z%A_pSWP(i,2) = 100.0d-3
					z%B_pSWP(i,2) = 1916.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'ethane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = 240.0d-3
					z%B_pSWP(i,2) = 2967.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'propane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = 5.455d-3
					z%B_pSWP(i,2) = 4638.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'i_butane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = 189.3d-3
					z%B_pSWP(i,2) = 3800.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'n_butane' ) then
					z%A_pSWP(i,1) = .0d-3
					z%B_pSWP(i,1) = 0.d0
					z%A_pSWP(i,2) = 30.51d-3
					z%B_pSWP(i,2) = 3699.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'carbon_dioxide' ) then
					z%A_pSWP(i,1) = .0845d-3
					z%B_pSWP(i,1) = 3615.d0
					z%A_pSWP(i,2) = 851.d-3
					z%B_pSWP(i,2) = 2025.d0
				else if( comp_list(z%guestindex(i))%comp%name_ == 'THF' ) then
					print*, 'need Clang par for THF @ ', __FILE__, __LINE__; stop
				else
					print*, 'unnexpected comp in hydrate phase model'; stop
				end if
			end do !i = 1, z%nguest
			
			z%delmu0 = 883.d0 !931.d0 !883.4148d0 !PP_1972 //	931 !Holder_1988(apudJPH1985)
			
			!(1):gelo (2):liq
			
			z%delh0(1) = 808.d0 !1400.d0 !tabela do artigo de 1985
			z%delh0(2) = -5201.d0 !z%delh0(1) - 6011.d0 !-5205.44844d0 !subtract 6011 !PP_1972 e !Holder_1988
			
			z%delcp0(1) = 0.d0 !-36.8607d0 !tabela do artigo de 1985
			z%delcp0(2) = -39.16d0 !									0.565 !Holder1980 (T<T0) (sI) //		1.029 !Holder_1988(apudJPH1985) (T<T0)
			z%q(1) = 0.d0 !0.14067648d0 !PP_1972 (sI) //	0.179 !Holder1980 (T>T0) (sI) //	0.1809 !Holder_1988(apudJPH1985) (T>T0)
			z%q(2) = 0.d0 !								0.002 !Holder1980 (T<T0) (sI) //	0.00377 !Holder_1988(apudJPH1985) (T<T0)

			z%delv0(1) = 3.4d-6! 3.4d-6 !PP_1972 // 3.39644d-6 !Holder_1988(apudJPH1985)
			z%delv0(2) = z%delv0(1) + 1.6d-6 !4.6d-6 !add 1.6 !Holder_1988(apudJPH1985)
			!
		case default
			print*, 'unexpected phase_cond argument in ', __FILE__, __LINE__; stop
			
		end select
		
		allocate(z%EOSslave, source = ceos_pure_i(pengrob_i()))

		call z%EOSslave%set_ceos_pure_index(z%water_index)

		z%EOSslave%phase_cond='L'

		allocate(z%ICEslave, source = ice_sublim_i())

		z%ICEslave%phase_cond='ICE'

		init_hyd_cond_fnum_defslavew = z
		
	end function init_hyd_cond_fnum_defslavew
!
	subroutine calcVEL_88(z,T,P)
		!declarations
		class(hydrate_cond_c) :: z
		real(8) :: T, P
		
		call z%calc_fugrefW(T,P)
		
!		select case(trim(z%phase_cond))
!		case('S1')
!			z%V_EL = 2.2588217455095876E-005	!S1		!só pra não ficar vazio por enquanto
!		case('S2')
!			z%V_EL = 2.3029068410175537E-005	!S2		!só pra não ficar vazio por enquanto
!		end select

		if (T <= z%T0) then
			z%V_EL = z%ICEslave%V + z%delv0(1)
		else 
			z%V_EL = z%EOSslave%V + z%delv0(2)
		end if
		
	end subroutine calcVEL_88

	subroutine calcFugWH_88(z,T,P)
		!declarations
		class(hydrate_cond_c) :: z
		real(8) :: T, P
		integer :: i, j
		!IMPLEMENT

		call z%calc_dmi_EL_refW(T,P)
		z%Fug_w_H = z%vdW_P_func*z%fugrefw*dexp(z%dmi_EL_LW_I/(R*T))
!		print*, z%vdW_P_func, __FILE__, __LINE__ !OFF
	end subroutine calcFugWH_88
	
!	subroutine calc_Phi_hyd2(z,T,P)
!		!declarations
!		class(hydrate_cond_c) :: z
!		real(8) :: T, P
!!		real(8) :: Phi(ncomp)
!!		integer :: i
!!		!definitions
!!		!IMPLEMENT
!!		call z%calcClang(T,P)												!this
!!		call z%calcOcc_dadoFUG(T,P) !calcula o x correspondente					!parent
!!		
!!!		call z%calcVEL(T,P)														!parent, (essa correlação não afeta a fugacidade desse modelo); mas para calcular ela eu teria que inicializar os parametros do outro

!!		
!!		call z%calcvdWPfunc(T,P)												!parent
!!		
!!		!metodos novos
!!		call z%calc_dmi_EL_refW(T,P)
!!		call z%calc_fugrefW(T,P)
!!		
!!		call z%calcFugWH(T,P)													!should access this's 
!!		!
!!!!!!		Phi(1) = z%fug_w_H/(z%x(1)*P)
!!!!!!		do i = 1, z%nguest
!!!!!!			if( z%x(i+1) > 1d-22 ) then
!!!!!!				Phi(i+1) = z%fug_g_H(i)/(z%x(i+1)*P)
!!!!!!			else
!!!!!!				z%x(i+1) = 1d-22
!!!!!!				Phi(i+1) = z%fug_g_H(i)/(1d-22*P)
!!!!!!!				print*, dlog(z%fug_g_H(i)/(1d-12*P)); pause
!!!!!!			end if
!!!!!!		end do !i = 1, nguest
!!!!!!		do i = 1+z%nguest+1, ncomp
!!!!!!			Phi(i) = dexp(50.d0)
!!!!!!		end do !i = 1+z%nguest+1, ncomp
!!!!!!		z%phi(:) = phi(:)
!		call z%hydrate_model%calc_phi() !access parent method to conclude the calculation
!	end subroutine calc_Phi_hyd2
end module hydrates_cond_mod

