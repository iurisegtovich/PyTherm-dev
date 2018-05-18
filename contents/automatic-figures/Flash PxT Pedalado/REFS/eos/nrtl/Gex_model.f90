module Gex_model_mod
	use comp_list_mod
	implicit none
	type Gex_model_c
		character(100) :: class_name = 'Gex_model'
		real(8), allocatable :: Gammai(:)
		!// initialized in the children classes' constructors
		contains
		procedure :: debug_Gex_model
		procedure :: calc_Gammai => calc_Gammai_deferred
	end type
	
	interface Gex_model_i
		module procedure init_Gex_model
	end interface
	
	contains
	
	function init_Gex_model()
		!ARGUMENTS
		type(Gex_model_c) :: z, init_Gex_model
		!LOCAL
		
		!implementation
		allocate(z%Gammai(ncomp))
		z%Gammai(:) = 1.d0
		init_Gex_model = z
	end function init_Gex_model
	
	subroutine debug_Gex_model(z)
		class(Gex_model_c) :: z
		print*, ''
		print*, 'debug_Gex_model ', z%class_name
		print*, 'z%Gammai(:)', z%Gammai(:)
	end subroutine debug_Gex_model
	
	subroutine calc_Gammai_deferred(z,T,P,x)
		!ARGUMENTS
		class(Gex_model_c) :: z
		real(8) :: T, P
		real(8) :: x(ncomp)
		!LOCAL
		print*, 'DEFERRED', __FILE__, __LINE__; stop
		return
	end subroutine calc_Gammai_deferred
end module
