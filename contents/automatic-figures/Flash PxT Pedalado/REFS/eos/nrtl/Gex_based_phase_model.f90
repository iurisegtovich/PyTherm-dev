module Gex_based_phase_model_mod
	use phase_model_mod !which does use comp_list_mod !which does use pureprop_mod
	use pur_liq_sat_model_mod
	use Gex_model_mod
	implicit none
	type, extends(phase_model_c) :: Gex_based_phase_model_c
!		character(5) :: phase_cond !L(liquida), V(vapor), ou S(mais estavel) para cúbicas, ...
		class(Gex_model_c), pointer :: Gex_model_obj !pra ter qualquer filho aqui tem que ser class, e para ser class tem que ser pointer[.]
		class(pur_liq_sat_model_c), pointer :: pur_liq_sat_model_slave !pra ter qualquer filho aqui tem que ser class, e para ser class tem que ser pointer[.]
!		// EoS parameters
		!// parameters Theta, b of mixture eos
!		real(8) :: Gex
		!// calcPhi method
		!//partial molar parameters (?)
		real(8), allocatable :: phi0(:)
		real(8), allocatable :: V0(:)
		contains
!		// constructor
!		// Gex_based_phase_model specific methods
		procedure :: calc_Volume_L
!		procedure :: calc_Gex

!		// pure virtual methods implemented
!		virtual void calcPhi()=0;
		procedure :: calc_phi => calc_phi_Gex_based_phase_model !não pode dar deferred se o tipo não for abstrato, e o tipo não pode ser abstrato se tiver alguma função normal (?)
		!misc
		procedure :: debug_phase_model => debug_Gex_based_phase_model
	end type Gex_based_phase_model_c
!
		interface Gex_based_phase_model_i
			module procedure allocate_Gex_based_phase_model
			module procedure init_Gex_based_phase_model_Gex_model_obj !apenas aloca e inicializa com valores neutros, requer Gex_model como argumento
			module procedure init_Gex_based_phase_model_fnum_Gex_model_obj !usa arquivo por número, requer Gex_model como argumento
		end interface
		
	contains
!	
	subroutine debug_Gex_based_phase_model(z)
		class(Gex_based_phase_model_c) :: z
		call z%phase_model_c%debug_phase_model
		print*, ''
		print*, 'debug_Gex_based_phase_model ', z%class_name
		print*, 'z%phi0(:)', z%phi0(:)
		print*, 'z%V0(:)', z%V0(:)
		if( associated(z%Gex_model_obj) ) then
			call z%Gex_model_obj%debug_Gex_model
		else
			print*, 'z%Gex_model_obj%debug_Gex_model not associated'
		end if
		if( associated(z%pur_liq_sat_model_slave) ) then
			call z%pur_liq_sat_model_slave%debug_phase_model
		else
			print*, 'z%pur_liq_sat_model_slave not associated'
		end if
	end subroutine debug_Gex_based_phase_model
	
	function allocate_Gex_based_phase_model()
		!declarations
		type(Gex_based_phase_model_c) :: allocate_Gex_based_phase_model
		!preamble
		allocate_Gex_based_phase_model%phase_model_c = allocate_phase_model()
		!implementation
!		// adapting vector lengths to number of components
		allocate_Gex_based_phase_model%class_name = 'Gex_based_phase_model'
		allocate(allocate_Gex_based_phase_model%phi0(ncomp))
		allocate(allocate_Gex_based_phase_model%V0(ncomp))
		allocate_Gex_based_phase_model%phi0(:) = 1.d0
		allocate_Gex_based_phase_model%V0(:) = 0.d0
		nullify(allocate_Gex_based_phase_model%Gex_model_obj)
!SLAVE
		allocate(allocate_Gex_based_phase_model%pur_liq_sat_model_slave, source = pur_liq_sat_model_i())
		return
		!***destroy z ao final ou não usar z e usar direto a variavel inicializada?
!			z não é allocatable então não deve dar memory leak
	end function allocate_Gex_based_phase_model
!
	function init_Gex_based_phase_model_Gex_model_obj(Gex_model_obj_in)
		type(Gex_based_phase_model_c) :: z, init_Gex_based_phase_model_Gex_model_obj
		!arguments
		class(Gex_model_c) :: Gex_model_obj_in
		!preamble and allocation
		z=allocate_Gex_based_phase_model()
		z%phase_model_c = phase_model_i()
		!implementation
		z%class_name = 'Gex_based_phase_model'
		z%phase_cond = 'L'
!OBJ
		allocate(z%Gex_model_obj, source = Gex_model_obj_in)

!TRANSLATE
		init_Gex_based_phase_model_Gex_model_obj = z !*** isso vai funcionar ou precisarei fazer tudo pointer e dar allocate init source = z ???
!		call init_Gex_based_phase_model_fnum_Gex_model_obj%Gex_model_obj%debug_pur_ceos_par !funcionou
	end function init_Gex_based_phase_model_Gex_model_obj
!
	function init_Gex_based_phase_model_fnum_Gex_model_obj(fnum_in,Gex_model_obj_in)
		type(Gex_based_phase_model_c) :: z, init_Gex_based_phase_model_fnum_Gex_model_obj
		!arguments
		class(Gex_model_c) :: Gex_model_obj_in
		integer :: fnum_in
		!LOCAL
		!preamble and allocation
		z=allocate_Gex_based_phase_model()
		z%phase_model_c = phase_model_i(fnum_in)
		!implementation
		z%class_name = 'Gex_based_phase_model'
!OBJ
		allocate(z%Gex_model_obj, source = Gex_model_obj_in)

!TRANSLATE
		init_Gex_based_phase_model_fnum_Gex_model_obj = z !*** isso vai funcionar ou precisarei fazer tudo pointer e dar allocate init source = z ???
!		call init_Gex_based_phase_model_fnum_Gex_model_obj%Gex_model_obj%debug_pur_ceos_par !funcionou
	end function init_Gex_based_phase_model_fnum_Gex_model_obj

!DEFERRED

!void Gex_based_phase_model::calcVolume()
	subroutine calc_Volume_L(z,T,P)
!declarations
		class(Gex_based_phase_model_c) :: z
!ARGUMENTS
		real(8) :: T, P
!LOCAL
		integer :: i
!IMPLEMENT
		z%V = 0.d0
		do i = 1, ncomp
			call z%pur_liq_sat_model_slave%set_pur_liq_sat_index(i)
			call z%pur_liq_sat_model_slave%calcV(T,P)
			z%V0(i) = z%pur_liq_sat_model_slave%V
			z%V = z%V + z%x(i)*z%V0(i)
		end do !i = 1, ncomp
		return
	end subroutine calc_Volume_L
	
!!!!!//void Gex_based_phase_model::calcPhi()
	subroutine calc_Phi_Gex_based_phase_model(z,T,P)
		!ARGUMENTS
		class(Gex_based_phase_model_c) :: z
		real(8) :: T, P
		integer :: i
		!LOCAL
		do i = 1, ncomp
			call z%pur_liq_sat_model_slave%set_pur_liq_sat_index(i)
			call z%pur_liq_sat_model_slave%calc_phi(T,P)
			z%phi0(i) = z%pur_liq_sat_model_slave%phi(i)
		end do !i = 1, ncomp
		
		call z%Gex_model_obj%calc_Gammai(T,P,z%x)
		z%phi(:) = z%Gex_model_obj%Gammai(:)*z%phi0(:)
	end subroutine calc_Phi_Gex_based_phase_model
	
	!
end module Gex_based_phase_model_mod

