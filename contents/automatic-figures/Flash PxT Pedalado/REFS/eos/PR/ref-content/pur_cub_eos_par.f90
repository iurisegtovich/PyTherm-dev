module pur_cub_eos_par_mod
	use comp_list_mod
!	use phase_model_mod
	implicit none
	type pur_cub_eos_par_c
		character(100) :: class_name = 'pur_cub_eos_par'
		!// initialized in the children classes' constructors
		real(8) :: sigma_, epsilon_ !cubic eos definers
		!// parameters alpha, a, b, theta, of pure compounds eos
		real(8), allocatable :: alphai(:), ai(:), bi(:), Thetai(:)
		contains
		procedure :: debug_pur_ceos_par
		procedure :: calc_T_dep => calc_T_dep_deferred
	end type
	
	interface pur_cub_eos_par_i
		module procedure init_pur_cub_eos_par
	end interface
	
	contains
	
		function init_pur_cub_eos_par()
		!declarations
		type(pur_cub_eos_par_c) :: z, init_pur_cub_eos_par
		!implementation
		allocate(z%alphai(ncomp), z%ai(ncomp), z%bi(ncomp), z%thetai(ncomp))
		z%alphai(:)=0.d0; z%ai(:)=0.d0; z%bi(:)= 0.d0; z%thetai(:)=0.d0
		init_pur_cub_eos_par = z
	end function init_pur_cub_eos_par
	
	subroutine debug_pur_ceos_par(z)
		class(pur_cub_eos_par_c) :: z
		print*, ''
		print*, 'debug_pur_ceos_par', z%class_name
		print*, 'z%sigma_', z%sigma_
		print*, 'z%epsilon_', z%epsilon_
		print*, 'z%alphai(:)', z%alphai(:)
		print*,  'z%ai(:)', z%ai(:)
		print*, 'z%bi(:)', z%bi(:)
		print*, 'z%Thetai(:)', z%Thetai(:)
	end subroutine debug_pur_ceos_par
	
	subroutine calc_T_dep_deferred(z,T)
		class(pur_cub_eos_par_c) :: z
		real(8) :: T
		!DEFERRED
		print*, 'DEFERRED', __FILE__, __LINE__; stop
	end subroutine calc_T_dep_deferred
end module
