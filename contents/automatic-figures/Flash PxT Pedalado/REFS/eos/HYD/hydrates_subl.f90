!hydrates_2.f90	baseado em Taiwan 2012
module hydrates_3_mod
	use hydrates_mod
	implicit none
	type, extends(hydrate_model_c) :: hydrate_model3
		!fazer lista de componentes própria
		!declarations
		real(8) :: Psat_w_EL, EL_Poynting
		real(8), allocatable :: V_SWP(:,:), V_SWP_0(:,:) !nguest,ncav
		real(8), allocatable :: depth_SWP(:)!(epsilon/kb) em [Kelvin] !nguest
		!parâmetros
		!por estrutura
		real(8) :: VELpar(6), PSATpar(5)
		!por cavidade
		real(8), allocatable :: VSWPstpar(:,:) !ncav,2
		!por comp e cavidade
		real(8), allocatable :: VSWPgpar(:,:,:) !nguest,ncav,4
		contains
!		// constructor
		procedure, nopass :: init_hydrate_model3
!		// hydrate_model3 specific methods
		procedure :: calcVSWP
		procedure :: calcClang => calcClang_SWP
		procedure :: calcVEL => calcVEL_KST
		procedure :: calcPsatEL
		procedure :: calcPoyntingEL
		procedure :: calcFugWH => calcFugWH_KST
!		// virtual methods
!		procedure :: calc_phi ! sem override
		procedure :: debugeos => debugeoshyd3
	end type hydrate_model3
!	
		interface hydrate_model3_i
			module procedure init_hydrate_model3
		end interface
		
	contains
!	
	subroutine debugeoshyd3(z)
		class(hydrate_model3) :: z
		print*, z%phase_cond!tudo
	end subroutine debugeoshyd3
	
	function init_hydrate_model3(phase_file_name)
		!declarations
		type(hydrate_model3) :: z, init_hydrate_model3
		character(100) :: phase_file_name
		integer :: i
		!preamble
		z%hydrate_model_c = hydrate_model_i(phase_file_name)
		
!		z%eos = eos_i(phase_file_name) !aqui se lê a flag de condição e a composição
		
		!implementation
		
!		call z%findguests()
		!
		allocate(z%V_SWP(z%nguest,2), z%V_SWP_0(z%nguest,2))
		allocate(z%depth_SWP(z%nguest))
		allocate(z%VSWPstpar(2,2),z%VSWPgpar(z%nguest,2,4))
!		open(unit=10,file=trim(phase_file_name),status='old',action='read')
!		read(10,*); do i = 1, ncomp; read(10,*); end do !discard comp
!		read(10,*); read(10,*) !z%phase_cond !discard cond
!		ler parâmetros
!		close(10)
		print*, z%phase_cond, 'phase_cond'
		select case(trim(z%phase_cond))
		case('S1')
			print*, 'param s1'
			!
			do i = 1, z%nguest
				if( comp_list(z%guestindex(i))%comp%name_ == 'methane' ) then
					z%V_SWP_0(i,1) = 44.592d-3; z%V_SWP_0(i,2) = 501.136d-3; !angstrom cubico
					z%depth_SWP(i) = 2548.81d0
					!
					z%VSWPgpar(i,1,1) = -11.171d0; z%VSWPgpar(i,1,2) = 14.983d-1; z%VSWPgpar(i,1,3) = -62.222d-3; z%VSWPgpar(i,1,4) = 1.197d-3;
					z%VSWPgpar(i,2,1) = -15.048d0; z%VSWPgpar(i,2,2) = 18.031d-1; z%VSWPgpar(i,2,3) = -59.931d-3; z%VSWPgpar(i,2,4) = 0.739d-3;
					
				else if( comp_list(z%guestindex(i))%comp%name_ == 'ethane' ) then
					z%V_SWP_0(i,1) = 0.d0; z%V_SWP_0(i,2) = 72.0280d-3; !angstrom cubico
					z%depth_SWP(i) = 3683.4000d0
					!
					z%VSWPgpar(i,1,1) = 0.d0; z%VSWPgpar(i,1,2) = 0.d0; z%VSWPgpar(i,1,3) = 0.d0; z%VSWPgpar(i,1,4) = 0.d0;
					z%VSWPgpar(i,2,1) = 0.d0; z%VSWPgpar(i,2,2) = 0.d0; z%VSWPgpar(i,2,3) = 0.d0; z%VSWPgpar(i,2,4) = 0.d0;
					
				else if( comp_list(z%guestindex(i))%comp%name_ == 'propane' ) then
					z%V_SWP_0(i,1) = 0.d-3; z%V_SWP_0(i,2) = 0.d-3;
					z%depth_SWP(i) = 0.d0
					!
					z%VSWPgpar(i,1,1) = 0.d0; z%VSWPgpar(i,1,2) = 0.d0; z%VSWPgpar(i,1,3) = 0.d0; z%VSWPgpar(i,1,4) = 0.d0;
					z%VSWPgpar(i,2,1) = 0.d0; z%VSWPgpar(i,2,2) = 0.d0; z%VSWPgpar(i,2,3) = 0.d0; z%VSWPgpar(i,2,4) = 0.d0;
				
				else if( comp_list(z%guestindex(i))%comp%name_ == 'carbon_dioxide' ) then
					z%V_SWP_0(i,1) = 0.1840d-3; z%V_SWP_0(i,2) = 160.5620d-3; !angstrom cubico
					z%depth_SWP(i) = 3222.9400d0
					!
					z%VSWPgpar(i,1,1) = 0.d0; z%VSWPgpar(i,1,2) = 0.d0; z%VSWPgpar(i,1,3) = 0.d0; z%VSWPgpar(i,1,4) = 0.d0;
					z%VSWPgpar(i,2,1) = 0.d0; z%VSWPgpar(i,2,2) = 0.d0; z%VSWPgpar(i,2,3) = 0.d0; z%VSWPgpar(i,2,4) = 0.d0;
				
				elseif( comp_list(z%guestindex(i))%comp%name_ == 'n_butane' ) then
					z%V_SWP_0(i,1) = 0.d-3; z%V_SWP_0(i,2) = 0.d-3; !angstrom cubico
					z%depth_SWP(i) = 0.d0
					!
					z%VSWPgpar(i,1,1) = 0.d0; z%VSWPgpar(i,1,2) = 0.d-1; z%VSWPgpar(i,1,3) = 0.d-3; z%VSWPgpar(i,1,4) = 0.d-3;
					z%VSWPgpar(i,2,1) = 0.d0; z%VSWPgpar(i,2,2) = 0.d-1; z%VSWPgpar(i,2,3) = 0.d-3; z%VSWPgpar(i,2,4) = 0.d-3;
					
				elseif( comp_list(z%guestindex(i))%comp%name_ == 'i_butane' ) then
					z%V_SWP_0(i,1) = 0.d-3; z%V_SWP_0(i,2) = 0.d-3; !angstrom cubico
					z%depth_SWP(i) = 0.d0
					!
					z%VSWPgpar(i,1,1) = 0.d0; z%VSWPgpar(i,1,2) = 0.d-1; z%VSWPgpar(i,1,3) = 0.d-3; z%VSWPgpar(i,1,4) = 0.d-3;
					z%VSWPgpar(i,2,1) = 0.d0; z%VSWPgpar(i,2,2) = 0.d-1; z%VSWPgpar(i,2,3) = 0.d-3; z%VSWPgpar(i,2,4) = 0.d-3;
					
				else
					print*, 'unnexpected comp in hydrate phase model'; stop
				end if
			end do !i = 1, z%nguest
			!
			!str s1
			z%VELpar(1) = 11.835d0; z%VELpar(2) = 2.217d-5; z%VELpar(3) = 2.242d-6; z%VELpar(4) = 0.d0; z%VELpar(5) = -8.006d-3; z%VELpar(6) = 5.448d0;
!			z%PSATpar(1) = 28.966d0; z%PSATpar(2) = 0.12879d0; z%PSATpar(3) = 2.1434d-3; z%PSATpar(4) = -6003.9d0; z%PSATpar(5) = 10002.d0;
	!mudança de sinal em 2 parÂmetros
			z%PSATpar(1) = 28.966d0; z%PSATpar(2) = 0.12879d0; z%PSATpar(3) = -2.1434d-3; z%PSATpar(4) = -6003.9d0; z%PSATpar(5) = -10002.d0;
			!str s1 (2 cavidades x 2 parameteros)
			z%VSWPstpar(1,1) = 2.786d-6; z%VSWPstpar(1,2) = 0.72d-6
			z%VSWPstpar(2,1) = 0.668d-6; z%VSWPstpar(2,2) = 49979.d-6

		case('S2')
			!
			do i = 1, z%nguest
				if( comp_list(z%guestindex(i))%comp%name_ == 'methane' ) then
					z%V_SWP_0(i,1) = 0.d-3; z%V_SWP_0(i,2) = 0.d-3; !angstrom cubico
					z%depth_SWP(i) = 0.d0
					!
					z%VSWPgpar(i,1,1) = 0.d0; z%VSWPgpar(i,1,2) = 0.d-1; z%VSWPgpar(i,1,3) = 0.d-3; z%VSWPgpar(i,1,4) = 0.d-3;
					z%VSWPgpar(i,2,1) = 0.d0; z%VSWPgpar(i,2,2) = 0.d-1; z%VSWPgpar(i,2,3) = 0.d-3; z%VSWPgpar(i,2,4) = 0.d-3;
					print*, __FILE__, __LINE__; stop
				else if( comp_list(z%guestindex(i))%comp%name_ == 'ethane' ) then
					
				else if( comp_list(z%guestindex(i))%comp%name_ == 'propane' ) then
!					z%V_SWP_0(i,1) = 0.d-3; z%V_SWP_0(i,2) = 14.222d-3;
					z%V_SWP_0(i,1) = 0.d-3; z%V_SWP_0(i,2) = 150.d-3;
!					z%depth_SWP(i) = 5128.5d0 !5.3d3 !GRÁFICO BOM !5128.5d0 !PARAMETRO ORIGINAL
					z%depth_SWP(i) = 4628.5d0
					!
					z%VSWPgpar(i,1,1) = 0.d0; z%VSWPgpar(i,1,2) = 0.d0; z%VSWPgpar(i,1,3) = 0.d0; z%VSWPgpar(i,1,4) = 0.d0;
					z%VSWPgpar(i,2,1) = 0.d0; z%VSWPgpar(i,2,2) = 0.d0; z%VSWPgpar(i,2,3) = 0.d0; z%VSWPgpar(i,2,4) = 0.d0;
					
				else if( comp_list(z%guestindex(i))%comp%name_ == 'n_butane' ) then
				
				else if( comp_list(z%guestindex(i))%comp%name_ == 'i_butane' ) then
				
				else if( comp_list(z%guestindex(i))%comp%name_ == 'carbon_dioxide' ) then
				
				else
					print*, 'unnexpected comp in hydrate phase model'; stop
				end if
			end do !i = 1, z%nguest
			
			!str s2																	z%VELpar(4) = 1.009d-9 na klauda e sandler 2000
			z%VELpar(1) = 17.13d0; z%VELpar(2) = 2.249d-4; z%VELpar(3) = 2.013d-6; z%VELpar(4) = 0.d0; z%VELpar(5) = -8.006d-3; z%VELpar(6) = 5.448d0;
			z%PSATpar(1) = 28.858d0; z%PSATpar(2) = 0.d0; z%PSATpar(3) = 0.d0; z%PSATpar(4) = -6017.6d0; z%PSATpar(5) = 0.d0;
			!str s2
			z%VSWPstpar(1,1) = 0.d0; z%VSWPstpar(1,2) = 0.d0
			z%VSWPstpar(2,1) = -47850.d-6; z%VSWPstpar(2,2) = 635315.d-6
			
		case default
			print*, 'unexpected phase_cond argument in ', __FILE__, __LINE__; stop
		end select
		
		init_hydrate_model3 = z
		
	end function init_hydrate_model3
!
	subroutine calcVSWP(z,T,P)
		!declarations
		class(hydrate_model3) :: z
		real(8) :: T, P
		!IMPLEMENT
		integer :: i, j
		real(8) :: Predw
		Predw = P/comp_list(1)%comp%Pc	!água é o comp 1 nessa versão
		do i = 1, z%nguest
			do j = 1, 2 !ncav
!			print*, z%VSWPstpar(j,1)*Predw, z%VSWPstpar(j,2)*Predw*Predw, 'dexp', z%VSWPgpar(i,j,1), z%VSWPgpar(i,j,2)*Predw, z%VSWPgpar(i,j,3)*Predw*Predw, z%VSWPgpar(i,j,4)*Predw*Predw*Predw
!			pause !OFF
				z%V_SWP(i,j) = z%V_SWP_0(i,j)/ &  !em angstrom cubico
(1.d0+z%VSWPstpar(j,1)*Predw+z%VSWPstpar(j,2)*Predw*Predw + &
dexp(z%VSWPgpar(i,j,1)+z%VSWPgpar(i,j,2)*Predw+z%VSWPgpar(i,j,3)*Predw*Predw+z%VSWPgpar(i,j,4)*Predw*Predw*Predw))
			end do !j = 1, 2 !ncav
		end do !i = 1, nguest
	end subroutine calcVSWP
!
	subroutine calcClang_SWP(z,T,P)
		!declarations
		class(hydrate_model3) :: z
		real(8) :: T, P
		!IMPLEMENT
		real(8) :: Kb = 1.3806488d-23
		integer :: i, j
		
		call z%calcVSWP(T,P)
		
		do i = 1, z%nguest
			do j = 1, 2 !ncav
				z%Clang(i,j) = ((4.d0*3.14159d0)/(Kb*T))*1d-30*z%V_SWP(i,j)*dexp(z%depth_SWP(i)/T)
				!print*, z%Clang(i,j) ![Pa]
			end do !j = 1, 2 !ncav
		end do !i = 1, nguest
	end subroutine calcClang_SWP
!
	subroutine calcPsatEL(z,T,P)
		!declarations
		class(hydrate_model3) :: z
		real(8) :: T, P
		integer :: i, j
!		real(8) :: 
		!IMPLEMENT
		!z%Psat_w_EL [Pa]
!		z%Psat_w_EL = dexp(z%PSATpar(1)+z%PSATpar(2)*dlog(T)+z%PSATpar(3)*T+z%PSATpar(4)/(T)+z%PSATpar(5)/(T*T))
		select case (trim(z%phase_cond))
		case ('S1')
!		z%Psat_w_EL = dexp(4.6477d0*dlog(T)-5242.979d0/T+2.7789d0-8.7156d-3*T) 
!!!!		print*, z%Psat_w_EL, 'klaudac1'
!!!!		z%Psat_w_EL = dexp(4.6766d0*dlog(T)-5263.9565d0/T+2.7789d0-9.0154d-3*T)
!!!!		print*, z%Psat_w_EL, 'klaudac2'
!!!!		z%Psat_w_EL = dexp(4.6652d0*dlog(T)-5424.1108d0/T+2.7789d0-8.8658d-3*T)
!!!!		print*, z%Psat_w_EL, 'klaudacc3h6'
!!!!		z%Psat_w_EL = dexp(4.6188d0*dlog(T)-5020.8289d0/T+2.7789d0-8.3455d-3*T)
!!!!		print*, z%Psat_w_EL, 'klaudaco2'
!!!!		z%Psat_w_EL = dexp(4.6446d0*dlog(T)-5150.3690d0/T+2.7789d0-8.7553d-3*T)
!!!!		print*, z%Psat_w_EL, 'klaudah2s'
			z%Psat_w_EL = dexp(z%PSATpar(1)+z%PSATpar(2)*dlog(T)+z%PSATpar(3)*T+z%PSATpar(4)/(T)+z%PSATpar(5)/(T*T))
!print*, z%PSATpar(1), z%PSATpar(2)*dlog(T), z%PSATpar(3)*T, z%PSATpar(4)/(T), z%PSATpar(5)/(T*T)
!pause !OFF
!!!!		print*, z%Psat_w_EL, 'taiwan s1'
		case ('S2')
!			z%Psat_w_EL = dexp(5.2578d0*dlog(T)-5650.5584d0/T+2.7789d0-16.2021d-3*T)
!			print*, z%Psat_w_EL, 'klauda s2 C3'
			z%Psat_w_EL = dexp(z%PSATpar(1)+z%PSATpar(2)*dlog(T)+z%PSATpar(3)*T+z%PSATpar(4)/(T)+z%PSATpar(5)/(T*T))
!print*, z%PSATpar(1), z%PSATpar(2)*dlog(T), z%PSATpar(3)*T, z%PSATpar(4)/(T), z%PSATpar(5)/(T*T)
!			print*, z%Psat_w_EL, 'taiwan s2'
!			pause
		end select
		
!		t=223.d0 !			debug, PLOT PSAT_EL
!		do i = 1, 100
!			write(91,*)	dexp(4.6477d0*dlog(T)-5242.979d0/T+2.7789d0-8.7156d-3*T), &
!						dexp(4.6766d0*dlog(T)-5263.9565d0/T+2.7789d0-9.0154d-3*T), &
!						dexp(4.6652d0*dlog(T)-5424.1108d0/T+2.7789d0-8.8658d-3*T), &
!						dexp(4.6188d0*dlog(T)-5020.8289d0/T+2.7789d0-8.3455d-3*T), &
!						dexp(4.6446d0*dlog(T)-5150.3690d0/T+2.7789d0-8.7553d-3*T)

!!			write(91,*)	dexp(5.1511d0*dlog(T)-5595.4346d0/T+2.7789d0-16.0445d-3*T), &
!!						dexp(5.2578d0*dlog(T)-5650.5584d0/T+2.7789d0-16.2021d-3*T), &
!!						dexp(4.6818d0*dlog(T)-5455.2664d0/T+2.7789d0-8.9678d-3*T), &
!!						dexp(5.1449d0*dlog(T)-5544.3272d0/T+2.7789d0-14.8446d-3*T)

!!			write(91,*) dexp(z%PSATpar(1)+z%PSATpar(2)*dlog(T)+z%PSATpar(3)*T+z%PSATpar(4)/(T)+z%PSATpar(5)/(T*T))
!!!			write(91,*) t, dexp(4.6477d0*dlog(T)-5242.979d0/T+2.7789d0-8.7156d-3*T), 'klauda', dexp(z%PSATpar(1)+z%PSATpar(2)*dlog(T)+z%PSATpar(3)*T+z%PSATpar(4)/(T)+z%PSATpar(5)/(T*T)), 'taiwan'
!!!			write(91,*) t, 
!!!			write(91,*) t, dexp(4.6477d0*dlog(T)-5242.9790d0/T+2.7789d0-8.7156d-3*T) !klauda CH4 s1
!			t=t+1.d0
!		end do !1 = 1, 100
!		stop
		
	end subroutine calcPsatEL

	subroutine calcVEL_KST(z,T,P)
		!declarations
		class(hydrate_model3) :: z
		real(8) :: T, P
		integer :: i, j
		real(8) :: g_w_ratio
!		real(8) :: dummy(3)
		!IMPLEMENT
!		z%V_EL = 1.92d-5+8.387d-10*T+4.016d-12*T*T; print*, z%V_EL, 'gelo'
		!z%V_EL [m3/mol_water]:
!!!!		T = 273.15d0			teste para um plot
!!!!		P = 10.d5			teste para um plot
		z%V_EL = ((z%VELpar(1)+z%VELpar(2)*T+z%VELpar(3)*T*T+z%VELpar(4)*T*T*T)**3)*(6.023d-7/z%UCnw) + z%VELpar(5)*1.d-6*1.d-6*P + z%VELpar(6)*1.d-6*1.d-6*1.d-6*1.d-6*P*P
!		print*, z%V_EL, 'hidrato vazio'
!!!!		dummy(1) = z%V_EL			teste para um plot
!!!!		call calcPoyntingH(z,T,P)			teste para um plot
!!!!		dummy(1) = z%EL_Poynting			teste para um plot
!!!!		z%velpar(4) = 1.009d-9			teste para um plot
!!!!!		z%V_EL = ((z%VELpar(1)+z%VELpar(2)*T+z%VELpar(3)*T*T+z%VELpar(4)*T*T*T)**3)*(6.023d-7/z%UCnw) + z%VELpar(5)*1.d-6*1.d-6*P + z%VELpar(6)*1.d-6*1.d-6*1.d-6*1.d-6*P*P
!!!!		dummy(2) = z%V_EL			teste para um plot
!!!!		call calcPoyntingH(z,T,P)			teste para um plot
!!!!		dummy(2) = z%EL_Poynting			teste para um plot
!!!!		print*, dummy(1), dummy(2); stop			teste para um plot

	end subroutine calcVEL_KST
	
	subroutine calcPoyntingEL(z,T,P)
		!declarations
		class(hydrate_model3) :: z
		real(8) :: T, P
		integer :: i, j
		real(8) :: term1, term2, term3
		!IMPLEMENT
		term1 = ((z%VELpar(1)+z%VELpar(2)*T+z%VELpar(3)*T*T+z%VELpar(4)*T*T*T)**3)*(6.023d-7/z%UCnw)
		term2 = z%VELpar(5)*1.d-6*1.d-6/2.d0
		term3 = z%VELpar(6)*1.d-6*1.d-6*1.d-6*1.d-6/3.d0
		
		z%EL_Poynting = (term1*P+term2*P*P+term3*P*P*P) - &
(term1*z%Psat_w_EL+term2*z%Psat_w_EL*z%Psat_w_EL+term3*z%Psat_w_EL*z%Psat_w_EL*z%Psat_w_EL)
		
		z%EL_Poynting = dexp(z%EL_Poynting/(8.31446*T))
		
	end subroutine calcPoyntingEL
	
	subroutine calcFugwH_KST(z,T,P)
		!declarations
		class(hydrate_model3) :: z
		real(8) :: T, P
		integer :: i, j
		!IMPLEMENT
		call z%calcPsatEL(T,P)
		call z%calcPoyntingEL(T,P)
		z%Fug_w_H = z%vdW_P_func*z%Psat_w_EL*z%EL_Poynting !*phi_sat?
	end subroutine calcFugwH_KST
	
	
!	subroutine calc_Phi_hyd(z,T,P)
!		!declarations
!		class(hydrate_model3) :: z
!		real(8) :: T, P
!!		real(8) :: Phi(ncomp)
!!		integer :: i
!!		!definitions
!!		!IMPLEMENT
!!		call z%calcVSWP(T,P)
!!		call z%calcClang(T,P)
!!		call z%calcOcc_dadoFUG(T,P) !calcula o x correspondente
!!		call z%calcVEL(T,P)			!chama depois do occ q ele calcula o V ocupado tb
!!		call z%calcPsatEL(T,P)
!!		call z%calcPoyntingEL(T,P)
!!		call z%calcvdWPfunc(T,P)
!!		call z%calcFugWH(T,P)
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
!		call z%hydrate_model%calc_phi !access parent method to conclude the calculation
!	end subroutine calc_Phi_hyd
end module hydrates_3_mod

