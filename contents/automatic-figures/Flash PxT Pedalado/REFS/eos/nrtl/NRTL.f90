module NRTL_mod
	use Gex_model_mod !which does use eos_mod !which does use comp_list_mod !which does use pureNRTLop_mod
	implicit none
	type, extends(Gex_model_c) :: NRTL_c
		real(8) :: alpha = 0.2d0
		real(8), allocatable :: Aij(:,:)

		contains
		procedure :: debug_Gex_model => debug_NRTL
		procedure :: calc_Gammai => calc_Gammai_NRTL
	end type NRTL_c
!	
	interface NRTL_i
		module procedure init_NRTL !parâmetros registrados no codigo fonte
	end interface
!
	contains
!
	function init_NRTL()
		!declarations
		!ARGUMENTS
		type(NRTL_c) :: z, init_NRTL
		!LOCAL
		integer :: i
		character(100) :: NRTLpar_fname
		!Preamble
		z%Gex_model_c = Gex_model_i()
		!implementation
		allocate(z%Aij(ncomp,ncomp))
		z%class_name = 'NRTL'

		!leitura de parâmetros em arquivo:
		NRTLpar_fname = 'input/parameters/'//'NRTL'//'/'//'Aij.dat'
		call read_asym_binary_matrix(NRTLpar_fname,z%Aij)

		init_NRTL = z
	end function init_NRTL
!
	subroutine debug_NRTL(z)
		class(NRTL_c) :: z
		call z%Gex_model_c%debug_Gex_model
		print*, ''
		print*, 'debug_NRTL ', z%class_name
		print*, 'z%alpha', z%alpha
		print*, 'z%Aij(:,:)', z%Aij(:,:)
	end subroutine debug_NRTL
	
	subroutine calc_Gammai_NRTL(z,T,P,x)
		!ARGUMENTS
		class(NRTL_c) :: z
		real(8) T, P
		real(8) :: x(ncomp)
		!LOCAL
		real(8) :: G(ncomp,ncomp)
		real(8) :: Lambda(ncomp,ncomp)
		real(8) :: E0(ncomp,ncomp)
		real(8) :: L0(ncomp,ncomp)

!		regra de indices do reshape:
!		Aij = transpose(reshape((/11,12,13,21,22,23,31,32,33/),(shape(Aij))))

		!exemplo
!		z%Aij = transpose(reshape((/0.d0, 107.99d0, 1011.98d0, 555.81d0, 0.d0, -1113.1d0, 2277.37d0, 1217.37d0, 0.d0/),(shape(z%Aij))))

		!CALCULATIONS
		G(:,:) = dexp(-1*z%alpha*z%Aij(:,:)/T)
		Lambda(:,:) = z%Aij(:,:)*G(:,:)/T
		E0 = matmul(Lambda,diag(1.d0/(matmul(transpose(G),x))))
		L0 = matmul(G,diag(1.d0/(matmul(transpose(G),x))))

		z%Gammai(:) = dexp( &		!1
		matmul( &						!2
		(E0+transpose(E0)-matmul( &		!3,4
		L0,matmul( &					!5
		diag(x),transpose(E0) &
		) &								!5
		)) &							!4,3
		,x) &							!2
		)								!1

!		print*, 'x', x(:), __FILE__, __LINE__
!		print*, 'gamma', z%Gammai(:), __FILE__, __LINE__; stop

		return
	
	end subroutine calc_Gammai_NRTL
	
	function diag(vector_in)
		real(8), intent(in) :: vector_in(:)
		real(8), allocatable :: diag(:,:)
		integer :: i, size
!		Msize = size(vector_in)
!		allocate(diag(Msize,Msize))
		allocate(diag(size(vector_in),size(vector_in)))
		diag = 0.d0
		do i = 1, size(vector_in)
			diag(i,i) = vector_in(i)
		end do !i = 1, size(vector_in)
		return
	end function diag
	
end module NRTL_mod
