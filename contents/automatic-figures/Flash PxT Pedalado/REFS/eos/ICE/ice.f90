module ice_mod
	use phase_model_mod !which does use comp_list_mod !which does use pureprop_mod
!	default da slave water vapor model
	use ceos_pure_mod
	use pengrob_mod
	implicit none
	type, extends(phase_model_c) :: ice_sublim_c
		!fazer lista de componentes própria
		class(ceos_pure_c), pointer :: EOSslave
		!declarations
		real(8) :: fug_w_ice, Psat_w_ice, V_ice, ice_Poynting, phisat
		!parâmetros
		!fixos, geométricos
		!por estrutura
		real(8) :: Vpar(3), PSATpar(5)
		integer :: water_index
		contains
!		// ice_sublim specific methods
		procedure :: calcV_ice
		procedure :: calcPsatice
		procedure :: calcPoyntingice
		procedure :: calcFugWice
!		// virtual methods
		procedure :: calc_phi => calc_phi_ice
		procedure :: debug_phase_model => debug_ice
	end type ice_sublim_c
!	
		interface ice_sublim_i
			module procedure init_ice_sublim_def !default tudo
			module procedure init_ice_sublim_fnum_defslavewV !baseado em numero de arquivo, modelo default para vapor saturado
		end interface
		
	contains
!	
	subroutine debug_ice(z)
		class(ice_sublim_c) :: z
		call z%phase_model_c%debug_phase_model
		print*, ' '
		print*, 'debug_ice ', z%class_name
		print*, 'z%fug_w_ice, z%Psat_w_ice, z%V_ice, z%ice_Poynting, z%phisat', z%fug_w_ice, z%Psat_w_ice, z%V_ice, z%ice_Poynting, z%phisat
		print*, 'z%Vpar(:), z%PSATpar(:)', z%Vpar(:), z%PSATpar(:) !***legendar
		print*, 'z%water_index', z%water_index
		call z%EOSslave%debug_phase_model
	end subroutine debug_ice
	
	function init_ice_sublim_fnum_defslavewV(fnum_in)
		!declarations
		type(ice_sublim_c) :: z, init_ice_sublim_fnum_defslavewV
		integer :: i, file_length, phase_file_number, fnum_in
		!preamble
		z%phase_model_c = phase_model_i(fnum_in)
		!implementation
		z%phase_cond='ICE'
		z%class_name = 'ice_sublim'
		!detect water_index
		z%water_index = 1 !***falta generalizar
!		print*, z%phase_cond, 'phase_cond'

		z%Vpar(1) = 1.912d-5; z%Vpar(2) = 8.387d-10; z%Vpar(3) = 4.016d-12;
!psat:
!parametros antigos
!		z%PSATpar(1) = 2.9446d0; z%PSATpar(2) = 4.6056d0; z%PSATpar(3) = -8.1431d-3; z%PSATpar(4) = -5501.1243d0; z%PSATpar(5) = 0.d0;
!parametros meus
		z%PSATpar(1) = 0.179896d2; z%PSATpar(2) = 0.201646d1; z%PSATpar(3) = -.441051d-2; z%PSATpar(4) = -.592242d4; z%PSATpar(5) = 0.d0;
		
		allocate(z%EOSslave, source = ceos_pure_i(pengrob_i()))
		call z%EOSslave%set_ceos_pure_index(z%water_index)
		
		z%EOSslave%phase_cond='V'
		
		init_ice_sublim_fnum_defslavewV = z
	end function init_ice_sublim_fnum_defslavewV
	
	function init_ice_sublim_def()
		!declarations
		type(ice_sublim_c) :: z, init_ice_sublim_def
		character(100) :: phase_file_name, Lwaterphase_file_name
		integer :: i, file_length, phase_file_number
		!preamble
		z%phase_model_c = phase_model_i()
		z%phase_cond='ICE'
		!detect water_index
		z%water_index = 1
		
		!initialize parameters
		z%Vpar(1) = 1.912d-5; z%Vpar(2) = 8.387d-10; z%Vpar(3) = 4.016d-12;
!psat:
!parametros antigos
!		z%PSATpar(1) = 2.9446d0; z%PSATpar(2) = 4.6056d0; z%PSATpar(3) = -8.1431d-3; z%PSATpar(4) = -5501.1243d0; z%PSATpar(5) = 0.d0;
!parametros meus
		z%PSATpar(1) = 0.179896d2; z%PSATpar(2) = 0.201646d1; z%PSATpar(3) = -.441051d-2; z%PSATpar(4) = -.592242d4; z%PSATpar(5) = 0.d0;
		
		!allocate slave eos
		allocate(z%EOSslave, source = ceos_pure_i(pengrob_i()))
		call z%EOSslave%set_ceos_pure_index(z%water_index)
		
		z%EOSslave%phase_cond='V'
!		print*, icephase_file_name; pause
		
		init_ice_sublim_def = z
	end function init_ice_sublim_def
!
	subroutine calcPsatice(z,T,P)
		!declarations
		class(ice_sublim_c) :: z
		real(8) :: T, P
		integer :: i, j
!		real(8) :: 
		!IMPLEMENT
		!z%Psat_w_EL [Pa]
		z%Psat_w_ice = dexp(z%PSATpar(1)+z%PSATpar(2)*dlog(T)+z%PSATpar(3)*T+z%PSATpar(4)/(T)+z%PSATpar(5)/(T*T))
		!o ponto triplo tá saindo um pouquinho errado pq aqui tá dando 613.49127586775921 a 273.16K, quando o experimental é 611.657
!		z%Psat_w_ice = 0.9973d0*z%Psat_w_ice !correção empírica pra 'cravar' o ponto triplo com os parametros antigos

	end subroutine calcPsatice

	subroutine calcV_ice(z,T,P)
		!declarations
		class(ice_sublim_c) :: z
		real(8) :: T, P
		integer :: i, j
		real(8) :: g_w_ratio
		!IMPLEMENT
!		z%V_EL = 1.92d-5+8.387d-10*T+4.016d-12*T*T; print*, z%V_EL, 'gelo'
		!z%V_EL [m3/mol_water]:
		z%V_ice = (z%Vpar(1)+z%Vpar(2)*T+z%Vpar(3)*T*T)

		
		z%V=z%V_ice
	end subroutine calcV_ice
	
	subroutine calcPoyntingice(z,T,P)
		!declarations
		class(ice_sublim_c) :: z
		real(8) :: T, P
		integer :: i, j
		real(8) :: term1, term2, term3
		!IMPLEMENT
		term1 = (z%Vpar(1)+z%Vpar(2)*T+z%Vpar(3)*T*T)
		
		z%ice_Poynting = (term1*P) - (term1*z%Psat_w_ice)
		
		z%ice_Poynting = dexp(z%ice_Poynting/(8.31446*T))
		
	end subroutine calcPoyntingice
	
	subroutine calcFugwice(z,T,P)
		!declarations
		class(ice_sublim_c) :: z
		real(8) :: T, P
		integer :: i, j
		!IMPLEMENT
		call z%EOSslave%set_ceos_pure_index(z%water_index)
		
		call z%EOSslave%calc_phi(T,z%Psat_w_ice)
		z%phisat = z%EOSslave%phi(1)
!		print*, T, z%phisat
		z%Fug_w_ice = z%Psat_w_ice*z%ice_Poynting*z%phisat
	end subroutine calcFugwice
	
	subroutine calc_Phi_ice(z,T,P)
		!declarations
		class(ice_sublim_c) :: z
		real(8) :: T, P
		real(8) :: Phi(ncomp)
		integer :: i
		!definitions
		!IMPLEMENT
		call calcV_ice(z,T,P)
		call calcPsatice(z,T,P)
		call calcPoyntingice(z,T,P)
		call calcFugWice(z,T,P)
		!
		z%x(1) = 1.d0
		!não mexo nos outros xizes, eles podem ser números suficientemente pequenos que satisfaçam a igualdade de fugacidade
		Phi(1) = z%fug_w_ice/P
		do i = 2, ncomp
			Phi(i) = dexp(20.d0) !com dexp(50 apareceram infinity em alguns momentos, e ao mudar para 20 convergiu)
!			if (z%x(i)>1.d-7) then
!				z%x(i) = 1.d-12
!				print*, 'gambiarra at', __LINE__, __FILE__;
!			end if
		end do !i = 1+z%nguest+1, ncomp
		z%phi(:) = phi(:)
	end subroutine calc_Phi_ice
end module ice_mod

