module pur_liq_sat_model_mod
	use phase_model_mod !which does use comp_list_mod !which does use pureprop_mod
!	default da slave comp vapor model
	use ceos_pure_mod
	use pengrob_mod
	implicit none
	type, extends(phase_model_c) :: pur_liq_sat_model_c
		!fazer lista de componentes própria
		class(ceos_pure_c), pointer :: EOSslave
		!declarations
		real(8) :: Psat_pur_liq_sat, pur_liq_sat_Poynting, phisat
		real(8) :: fug_pur_liq_sat
		!parâmetros
		!fixos, geométricos
		!por estrutura
		real(8), allocatable :: Vpar(:,:) !3,ncomp,3
		real(8), allocatable :: PSATpar(:,:) !5,ncomp,5
		!:os valores vão depender da substância e devem ser lidos em arquivo na init, baseado na variável pure_index
										!ou eu faço matrix (nparmax,ncomp)
		integer :: pure_index
		contains
!		// pur_liq_sat specific methods
		procedure :: calcV
		procedure :: calcPsatpur_liq_sat
		procedure :: calcPoyntingpur_liq_sat
		procedure :: calcFugWpur_liq_sat
		
		procedure :: set_pur_liq_sat_index
!		// virtual methods
		procedure :: calc_phi => calc_phi_pur_liq_sat
		procedure :: debug_phase_model => debug_pur_liq_sat
	end type pur_liq_sat_model_c
!	
		interface pur_liq_sat_model_i
			module procedure allocate_pur_liq_sat !default tudo
			module procedure init_pur_liq_sat_fnum_index_defslaveV !baseado em numero de arquivo, e compindex, modelo default para vapor saturado
		end interface
		
	contains
!	
	subroutine debug_pur_liq_sat(z)
		class(pur_liq_sat_model_c) :: z
		call z%phase_model_c%debug_phase_model
		print*, ' '
		print*, 'debug_pur_liq_sat ', z%class_name
		print*, 'z%fug_pur_liq_sat, z%Psat_pur_liq_sat, z%V, z%pur_liq_sat_Poynting, z%phisat', z%fug_pur_liq_sat, z%Psat_pur_liq_sat, z%V, z%pur_liq_sat_Poynting, z%phisat
		print*, 'z%pure_index', z%pure_index
		if( associated(z%EOSslave) ) then
			call z%EOSslave%debug_phase_model
		else
			print*, 'z%EOSslave not associated'
		end if
	end subroutine debug_pur_liq_sat
	
	subroutine set_pur_liq_sat_index(z,pure_index_in)
		class(pur_liq_sat_model_c) :: z
		integer :: pure_index_in
		z%pure_index = pure_index_in
		z%x(:) = 0.d0
		z%x(z%pure_index) = 1.d0
		
		call z%EOSslave%set_ceos_pure_index(z%pure_index)
	end subroutine set_pur_liq_sat_index

	function init_pur_liq_sat_fnum_index_defslaveV(fnum_in)
		!declarations
		type(pur_liq_sat_model_c) :: z, init_pur_liq_sat_fnum_index_defslaveV
		integer :: fnum_in
		!preamble
		init_pur_liq_sat_fnum_index_defslaveV=allocate_pur_liq_sat()
		z%phase_model_c = phase_model_i(fnum_in)
		!implementation
		
		init_pur_liq_sat_fnum_index_defslaveV = z
	end function init_pur_liq_sat_fnum_index_defslaveV
	
	function allocate_pur_liq_sat()
		!declarations
		type(pur_liq_sat_model_c) :: z, allocate_pur_liq_sat
		character(100) :: phase_file_name, Lcompphase_file_name
		integer :: i, file_length, phase_file_number
		!preamble
		z%phase_model_c = phase_model_i()
		z%class_name = 'pur_liq_sat'
		z%phase_cond = 'L'
		z%pure_index = 0 !não inicializado
		
		allocate(z%Vpar(ncomp,3))
		allocate(z%PSATpar(ncomp,5))
		
		print*, 'buscar parametros', __FILE__, __LINE__
		
		!water
		z%Vpar(1,1) = 1.912d-5; z%Vpar(1,2) = 8.387d-10; z%Vpar(1,3) = 4.016d-12;
		z%PSATpar(1,1) = 0.179896d2; z%PSATpar(1,2) = 0.201646d1; z%PSATpar(1,3) = -.441051d-2; z%PSATpar(1,4) = -.592242d4; z%PSATpar(1,5) = 0.d0;

		!ethanol
		z%Vpar(2,1) = 1.912d-5; z%Vpar(2,2) = 8.387d-10; z%Vpar(2,3) = 4.016d-12;
!		z%PSATpar(2,1) = 0.179896d2; z%PSATpar(2,2) = 0.201646d1; z%PSATpar(2,3) = -.441051d-2; z%PSATpar(2,4) = -.592242d4; z%PSATpar(2,5) = 0.d0;

!		!limonene
		z%Vpar(3,1) = 1.912d-5; z%Vpar(3,2) = 8.387d-10; z%Vpar(3,3) = 4.016d-12;
!		z%PSATpar(3,1) = 0.179896d2; z%PSATpar(3,2) = 0.201646d1; z%PSATpar(3,3) = -.441051d-2; z%PSATpar(3,4) = -.592242d4; z%PSATpar(3,5) = 0.d0;

!		!THF (Antoine eq.)
!		z%Vpar(2,1) = 1.912d-5; z%Vpar(2,2) = 8.387d-10; z%Vpar(2,3) = 4.016d-12;
!		z%PSATpar(2,1) = 0.179896d2; z%PSATpar(2,2) = 0.201646d1; z%PSATpar(2,3) = -.441051d-2; z%PSATpar(2,4) = -.592242d4; z%PSATpar(2,5) = 0.d0;

!		!n-heptane (Antoine eq.)
!		z%Vpar(3,1) = 1.912d-5; z%Vpar(3,2) = 8.387d-10; z%Vpar(3,3) = 4.016d-12;
!		z%PSATpar(3,1) = 0.179896d2; z%PSATpar(3,2) = 0.201646d1; z%PSATpar(3,3) = -.441051d-2; z%PSATpar(3,4) = -.592242d4; z%PSATpar(3,5) = 0.d0;

		allocate(z%EOSslave, source = ceos_pure_i(pengrob_i()))
		z%EOSslave%phase_cond = 'V'
		
		!allocate slave eos
		allocate(z%EOSslave, source = ceos_pure_i(pengrob_i()))
		z%EOSslave%phase_cond='V'
!		print*, pur_liq_satphase_file_name; pause
		
		allocate_pur_liq_sat = z
	end function allocate_pur_liq_sat
!
	subroutine calcPsatpur_liq_sat(z,T)
		!declarations
		class(pur_liq_sat_model_c) :: z
		real(8) :: T
		integer :: i, j
!		real(8) :: 
		!IMPLEMENT
		!z%Psat_EL [Pa]
		z%Psat_pur_liq_sat = dexp(z%PSATpar(1,1)+z%PSATpar(1,2)*dlog(T)+z%PSATpar(1,3)*T+z%PSATpar(1,4)/(T)+z%PSATpar(1,5)/(T*T))
		select case (z%pure_index)
!!		z%Psat_pur_liq_sat = 133.322*(10**((A)-((B)/(T-273.15+C)))) !Antoine eq
		case(1)
!		!water
		z%Psat_pur_liq_sat = 1.d5*133.322*(10.d0**((5.11564)-((1687.537)/(T-273.15+230.17))))
		case(2)
		!THF
		z%Psat_pur_liq_sat = 1.d5*133.322*(10.d0**((4.12142)-((1203.11)/(T-273.15+226.355)))) 
		!print*, 'z%Psat_pur_liq_sat', z%Psat_pur_liq_sat, __FILE__, __LINE__
		case(3)
		!n_heptane
		z%Psat_pur_liq_sat = 1.d5*133.322*(10.d0**((4.02023)-((1263.909)/(T-273.15+216.432))))  
!		case default
!			print*, 'sub ainda não implementada genericamente', __FILE__, __LINE__; stop
		end select
!		
!!		print*, 		z%Psat_pur_liq_sat
!		z%Psat_pur_liq_sat = dexp(z%PSATpar(z%pure_index,1)+z%PSATpar(z%pure_index,2)*dlog(T)+z%PSATpar(z%pure_index,3)*T+z%PSATpar(z%pure_index,4)/(T)+z%PSATpar(z%pure_index,5)/(T*T))

!		print*, 		z%Psat_pur_liq_sat
!		pause


	end subroutine calcPsatpur_liq_sat

	subroutine calcV(z,T,P)
		!declarations
		class(pur_liq_sat_model_c) :: z
		real(8) :: T, P
		integer :: i, j
		!IMPLEMENT
!		z%V_EL = 1.92d-5+8.387d-10*T+4.016d-12*T*T; print*, z%V_EL, 'gelo'
		!z%V_EL [m3/mol_comp]:
!		z%V = z%Vpar(z%pure_index,1)+z%Vpar(z%pure_index,2)*T+z%Vpar(z%pure_index,3)*T*T !para o compindex, z%Vpar(z%pure_index,1)
		select case (z%pure_index)
!		z%Psat_pur_liq_sat = 133.322*(10**((A)-((B)/(T-273.15+C)))) !Antoine eq
		case(1)
		!water
		z%V = (55.95d0)**(-1)*1.d-3 !(55.95*1.d-6)
		case(2)
		!THF
		z%V = (224.d0)**(-1)*1.d-3 !224*1.d-6
		case(3)
		!n_heptane
		z%V = (428.d0)**(-1)*1.d-3 !428*1d-6
		case default
			print*, 'sub ainda não implementada genericamente', __FILE__, __LINE__; stop
		end select
!		print*, 'V', z%V, __FILE__, _LINE__

	end subroutine calcV
	
	subroutine calcPoyntingpur_liq_sat(z,T,P)
		!declarations
		class(pur_liq_sat_model_c) :: z
		real(8) :: T, P
		integer :: i, j
		real(8) :: term1, term2, term3
		!IMPLEMENT
!		term1 = (z%Vpar(1)+z%Vpar(2)*T+z%Vpar(3)*T*T)
		term1 = z%V !já calculado
		
		z%pur_liq_sat_Poynting = (term1*P) - (term1*z%Psat_pur_liq_sat)
		
		z%pur_liq_sat_Poynting = dexp(z%pur_liq_sat_Poynting/(8.31446*T))
		
	end subroutine calcPoyntingpur_liq_sat
	
	subroutine calcFugwpur_liq_sat(z,T,P)
		!declarations
		class(pur_liq_sat_model_c) :: z
		real(8) :: T, P
		integer :: i, j
		!IMPLEMENT

		call z%EOSslave%calc_phi(T,z%Psat_pur_liq_sat)
		z%phisat = z%EOSslave%phi(z%pure_index)
!		print*, z%pure_index, z%EOSslave%pure_index, z%Psat_pur_liq_sat, z%phisat, z%EOSslave%V
!		print*, T, z%phisat
		z%Fug_pur_liq_sat = z%Psat_pur_liq_sat*z%pur_liq_sat_Poynting*z%phisat
	end subroutine calcFugwpur_liq_sat
	
	subroutine calc_Phi_pur_liq_sat(z,T,P)
		!declarations
		class(pur_liq_sat_model_c) :: z
		real(8) :: T, P
		real(8) :: Phi(ncomp)
		integer :: i
		!definitions
		!IMPLEMENT
		call calcV(z,T,P)
		call calcPsatpur_liq_sat(z,T)
		call calcPoyntingpur_liq_sat(z,T,P)
		call calcFugWpur_liq_sat(z,T,P)
		!
		!não mexo nos outros xizes, eles podem ser números suficientemente pequenos que satisfaçam a igualdade de fugacidade
		do i = 1, ncomp
			if(i==z%pure_index) then
				Phi(i) = z%fug_pur_liq_sat/P
			else
				Phi(i) = dexp(20.d0) !com dexp(50 apareceram infinity em alguns momentos, e ao mudar para 20 convergiu)
			endif
!			if (z%x(i)>1.d-7) then
!				z%x(i) = 1.d-12
!				print*, 'gambiarra at', __LINE__, __FILE__;
!			end if
		end do !i = 1+z%nguest+1, ncomp
		z%phi(:) = phi(:)
	end subroutine calc_Phi_pur_liq_sat
end module pur_liq_sat_model_mod

