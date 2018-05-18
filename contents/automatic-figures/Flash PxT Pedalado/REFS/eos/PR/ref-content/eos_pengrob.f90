module pengrob_mod
	use pur_cub_eos_par_mod !which does use eos_mod !which does use comp_list_mod !which does use pureprop_mod
	implicit none
	type, extends(pur_cub_eos_par_c) :: pengrob_c
		real(8), allocatable :: k(:)
		!parâmetros especiais para a água:
		real(8) :: bparam = 0.848496d0 !fugacidade a 272.d0K, 100.d5Pa dá 607.856Pa
		real(8) :: aparam = 0.960392d0
		real(8) :: alphaparam(3) = (/0.816384d0, 0.128223d1, 0.538174d-3/)

		contains
		procedure :: debug_pur_ceos_par => debug_PR
		procedure :: calc_T_dep => calc_T_dep_PR
	end type pengrob_c
!	
	interface pengrob_i
		module procedure init_pengrob !parâmetros registrados no codigo fonte
	end interface
!
	contains
!
	function init_pengrob()
		!declarations
		type(pengrob_c) :: z, init_pengrob
		integer :: i
		real(8), parameter :: sqrt2 = dsqrt(2.d0)
		!preamble
		z%pur_cub_eos_par_c = pur_cub_eos_par_i()
		!implementation
		z%class_name = 'PR'
		allocate(z%k(ncomp))
		z%sigma_ = 1.d0 + sqrt2	!material charlles
		z%epsilon_ = 1.d0 - sqrt2	!material charlles
		do i = 1, ncomp
			if( trim(comp_list(i)%comp%name_) == 'water' ) then
				z%ai(i) = z%aparam*0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
				z%bi(i) = z%bParam*0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc) !; print*, 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc), 'b'; pause
			else
				z%ai(i) = 0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
				z%bi(i) = 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc)
			end if
			z%k(i) = 0.37464d0 + 1.54226d0*comp_list(i)%comp%w-0.26992d0*(comp_list(i)%comp%w)**2
		end do
		init_pengrob = z
	end function init_pengrob
!
	subroutine calc_T_dep_PR(z,T)
		class(pengrob_c) :: z
		real(8) :: T
		integer :: i
		!preamble
		!implementation
		do i = 1, ncomp
			if( trim(comp_list(i)%comp%name_) == 'water' ) then
				if (T < comp_list(i)%comp%Tc)then
					z%alphai(i) = dexp(z%alphaparam(1)*(1.d0-T/comp_list(i)%comp%Tc) &		
* (dabs(1.d0-T/comp_list(i)%comp%Tc))**(z%alphaparam(2)-1.d0) &							
+z%alphaparam(3)*(comp_list(i)%comp%Tc/T-1.d0))											
				else
					z%alphai(i) = 1.d0
				end if
!			print*, z%alphai(i), 'it is water'; pause
			else
				z%alphai(i) = (1.d0 +z%k(i)*(1.d0-dsqrt(T/comp_list(i)%comp%Tc)))**2
			end if
			z%thetai(i) = z%ai(i)*z%alphai(i)
		end do
		!
	end subroutine calc_T_dep_PR
	
	subroutine debug_PR(z)
		class(pengrob_c) :: z
		print*, ' '
		print*, 'debug_PR'
		print*, 'k(:)', z%k(:)
		print*, 'z%aparam', z%aparam
		print*, 'z%bparam', z%bparam
		print*, 'z%alphaparam(:)', z%alphaparam(:)
		call z%pur_cub_eos_par_c%debug_pur_ceos_par
	end subroutine debug_PR
end module pengrob_mod

module pengrob_w_mod
	use pengrob_mod
	implicit none
	type, extends(pengrob_c) :: pengrob_w_c
!		real(8), allocatable :: k(:)
!		real(8) :: bparam = 0.848496d0 !fugacidade a 272.d0K, 100.d5Pa dá 607.856Pa
!		real(8) :: aparam = 0.960392d0
!		real(8) :: alphaparam(3) = (/0.816384d0, 0.128223d1, 0.538174d-3/)
!!		real(8) :: bparam = 1.d0 !fugacidade a 272.d0K, 100.d5Pa dá 615.404Pa
!!		real(8) :: aparam = 1.d0
!!		real(8) :: alphaparam(3) = (/0.81473d0, 0.96611d0, 0.02707d0/)

!		contains
!		procedure, nopass :: init_pengrob
!		procedure :: debug_phase_model => debug_PR
!		procedure :: calc_T_dep => calc_T_dep_PR
	end type pengrob_w_c
!!	
!	interface pengrob_i
!		module procedure init_pengrob, init_pengrob_0
!	end interface
!!
!	contains
!!
!	function init_pengrob(phase_file_name)
!		!declarations
!		type(pengrob_c) :: z, init_pengrob
!		character(100) :: phase_file_name !, mixingrule_file_name
!		integer :: i
!		!preamble
!		z%pur_cub_eos_par_c = pur_cub_eos_par_i(phase_file_name)
!		!implementation
!		z%class_name = 'PRw'
!		!ncomp = size(z%x) !ou size(comp_list)
!		allocate(z%k(ncomp))
!!		print*, ncomp, 'ncomp'
!		z%sigma_ = 1.d0 + sqrt2	!material charlles
!		z%epsilon_ = 1.d0 - sqrt2	!material charlles
!		do i = 1, ncomp
!			if( trim(comp_list(i)%comp%name_) == 'water' ) then
!				z%ai(i) = z%aparam*0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
!				z%bi(i) = z%bParam*0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc) !; print*, 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc), 'b'; pause
!			else
!				z%ai(i) = 0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
!				z%bi(i) = 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc)
!			end if
!			z%k(i) = 0.37464d0 + 1.54226d0*comp_list(i)%comp%w-0.26992d0*(comp_list(i)%comp%w)**2
!		end do
!		init_pengrob = z
!	end function init_pengrob
!!
!	function init_pengrob_0()
!		!declarations
!		type(pengrob_c) :: z, init_pengrob_0
!		integer :: i
!		!preamble
!		z%pur_cub_eos_par_c = pur_cub_eos_par_i()
!		!implementation
!		!ncomp = size(z%x) !ou size(comp_list)
!		allocate(z%k(ncomp))
!!		print*, ncomp, 'ncomp'
!		z%sigma_ = 1.d0 + sqrt2	!material charlles
!		z%epsilon_ = 1.d0 - sqrt2	!material charlles
!		do i = 1, ncomp
!			if( trim(comp_list(i)%comp%name_) == 'water' ) then
!				z%ai(i) = z%aparam*0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
!				z%bi(i) = z%bParam*0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc) !; print*, 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc), 'b'; pause
!			else
!				z%ai(i) = 0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
!				z%bi(i) = 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc)
!			end if
!			z%k(i) = 0.37464d0 + 1.54226d0*comp_list(i)%comp%w-0.26992d0*(comp_list(i)%comp%w)**2
!		end do
!		init_pengrob_0 = z
!	end function init_pengrob_0
!!
!	subroutine calc_T_dep_PR(z,T,P)
!		class(pengrob_c) :: z
!		real(8) :: T, P
!		integer :: i
!		!preamble
!		!implementation
!		do i = 1, ncomp
!			if( trim(comp_list(i)%comp%name_) == 'water' ) then
!				if (T < comp_list(i)%comp%Tc)then
!					z%alphai(i) = dexp(z%alphaparam(1)*(1.d0-T/comp_list(i)%comp%Tc) &		
!* (dabs(1.d0-T/comp_list(i)%comp%Tc))**(z%alphaparam(2)-1.d0) &							
!+z%alphaparam(3)*(comp_list(i)%comp%Tc/T-1.d0))											
!				else
!					z%alphai(i) = 1.d0
!				end if
!!			print*, z%alphai(i), 'it is water'; pause
!			else
!				z%alphai(i) = (1.d0 +z%k(i)*(1.d0-dsqrt(T/comp_list(i)%comp%Tc)))**2
!			end if
!			z%thetai(i) = z%ai(i)*z%alphai(i)
!		end do
!		!
!		!recycle parent procedure
!!		print*, 'passo1'
!		
!	end subroutine calc_T_dep_PR
!	
!	
!	subroutine debug_PR(z)
!		class(pengrob_c) :: z
!		print*, 'debugPR'
!		print*, z%sigma_
!		
!		print*, z%ai, 'debug ai' !OK
!		print*, z%bi !OK
!		
!	end subroutine debug_PR
end module pengrob_w_mod
