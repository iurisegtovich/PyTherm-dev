module hydrates_mod
	use phase_model_mod !which does use comp_list_mod !which does use pureprop_mod
	implicit none
	type, extends(phase_model_c) :: hydrate_model_c
		!fazer lista de componentes própria
		integer :: nguest
		!declarations
		real(8) :: vdW_P_func, fug_w_H, V_EL
		real(8), allocatable :: occ(:,:), Clang(:,:) !nguest,ncav
		real(8), allocatable :: fug_g_H(:) !nguest
		integer, allocatable :: guestindex(:) !tem tamanho máximo ncomp, mas só pretendo usar os espaços de 1 até z%nguest
		integer :: water_index !host
		!parâmetros
		!fixos, geométricos
		real(8) :: UCnw, ni(2) !ncav
		contains
!		// constructor
!		// hydrate_model specific methods
		procedure :: findguests
		procedure :: calcClang => calcClang_deferred !depende do modelo filho
		procedure :: calcOcc_dadoxiz !calcula o fug_g correspondente
		procedure :: calcOcc_dadoFUG !calcula o x correspondente
		procedure :: calcVEL => calcVEL_deferred !depende do modelo filho_deferred !método essencial para conhecer a variável volume da fase, mas depende do modelo filho
		procedure :: calcV => calcV_hyd
		procedure :: calcvdWPfunc
		procedure :: calcFugWH => calcFugWH_deferred !depende do modelo filho
!		// virtual methods
		procedure :: calc_phi => calc_phi_hyd !conclusão do calculo de phi, dadas as fug e os xizes
		procedure :: debug_phase_model => debug_hyd
	end type hydrate_model_c
!	
		interface hydrate_model_i
			module procedure :: init_hydrate_model_fnum !le condição fisica em arquivo pelo numero, usa default para agua liquida e gelo de referencias
!			module procedure :: allocate_hydrate_model
		end interface
		
	contains
!	
	subroutine debug_hyd(z)
		class(hydrate_model_c) :: z
		call z%phase_model_c%debug_phase_model
		print*,
		print*, 'debug_hyd ', z%class_name
		print*, 'z%nguest', z%nguest
!		vdW_P_func, fug_w_H, V_EL
!		occ(:,:), Clang(:,:)
!		fug_g_H(:)
!		guestindex(:)
	end subroutine debug_hyd
	
	subroutine findguests(z)
		class(hydrate_model_c) :: z
		character(100), allocatable :: possibleguestsnames(:)
		integer :: i, j, k
		integer :: listsize

		open(unit=10,file='input/parameters/hyd/guestlist.dat',status='old',action='read')
		read(10,*) !cabeçalho
		read(10,*) listsize
		allocate(possibleguestsnames(listsize))
		read(10,*) !cabeçalho
		do i = 1, listsize
			read(10,*) possibleguestsnames(i)
		end do !i = 1, listsize
		close(10)
		
		allocate(z%guestindex(ncomp))
		
		z%nguest = 0
		z%guestindex(:)=0
		j=1
		i=1
!		print*, listsize
		do
			if (j>listsize) then
				j = 1
				i = i + 1
				if (i>ncomp) exit
			end if
			
!			print*, comp_list(i)%comp%name_, possibleguestsnames(j)
			if (possibleguestsnames(j)==comp_list(i)%comp%name_) then
				z%nguest = z%nguest + 1
				z%guestindex(z%nguest) = i
				j = 0
				i = i + 1
				
				if (i>ncomp) exit
			end if
			j = j + 1
		end do
		
!		print*, z%nguest, guestindex(:); pause
	end subroutine findguests
	
	function init_hydrate_model_fnum(fnum_in)
		!declarations
		type(hydrate_model_c) :: z, init_hydrate_model_fnum
		integer :: i, fnum_in
		!preamble
		z%phase_model_c = phase_model_i(fnum_in) !aqui se lê a flag de condição e a composição
		
		!implementation
		
		!detect guests
		call z%findguests()
		
		!detect water_index
		z%water_index = 1 !***falta generalizar
		
		allocate(z%occ(z%nguest,2), z%Clang(z%nguest,2))
		allocate(z%fug_g_H(z%nguest))
		
!		print*, z%phase_cond, 'phase_cond'
		select case(trim(z%phase_cond))
		case('S1')
			z%ni(1) = 2.d0/46.d0; z%ni(2) = 6.d0/46.d0
			z%UCnw = 46.d0
			!
		case('S2')
			z%ni(1) = 16.d0/136.d0; z%ni(2) = 8.d0/136.d0
			z%UCnw = 136.d0
			!
		end select
		init_hydrate_model_fnum = z
		
	end function init_hydrate_model_fnum
	
!
	subroutine calcClang_deferred(z,T,P)
		!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		!IMPLEMENT
		print*, 'DEFERRED @', __FILE__, __LINE__; stop
	end subroutine calcClang_deferred
!
	subroutine calcOcc_dadoFUG(z,T,P)
!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		real(8) :: sum1, sum2, sum3
		integer :: i, j
!IMPLEMENT
		do i = 1, z%nguest
			z%fug_g_H(i) = z%fug_aux(z%guestindex(i))
!			print*, z%fug_g_H(i), __FILE__, __LINE__ !OFF
		end do
		
		do j = 1, 2
			sum1 = 0.d0
			do i = 1, z%nguest
				sum1 = sum1+z%Clang(i,j)*z%fug_g_H(i)
			end do !i = 1, nguest
			do i = 1, z%nguest
				z%Occ(i,j) = z%Clang(i,j)*z%fug_g_H(i)/(1.d0+sum1)
!				print*, z%Occ(i,j), z%Clang(i,j), z%fug_g_H(i), sum1, __FILE__, __LINE__
			end do !i = 1, nguest
		end do !j = 1, 2
		!calcular os xizes também:
		sum2 = 0.d0
		do i = 1, z%nguest
			do j = 1, 2
				sum2 = sum2 + z%ni(j)*z%occ(i,j)
			end do !j = 1, 2
		end do !i = 1, nguest

		sum3 = 0.d0
		do i = 1, z%nguest
		sum1 = 0.d0
			do j = 1, 2
				sum1 = sum1 + z%ni(j)*z%occ(i,j)
			end do !j = 1, 2
		z%x(z%guestindex(i)) = sum1/(1.d0+sum2)
		sum3 = sum3 + z%x(z%guestindex(i))
		end do !i = 1, nguest
		z%x(1) = 1.d0 - sum3
		
	end subroutine calcOcc_dadoFUG

	subroutine calcOcc_dadoxiz(z,T,P)
!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		real(8) :: sum1
		integer :: i, j
!IMPLEMENT
!		calcular as fug resolvendo o sistema não linear
		!!!
		print*, 'men at work', __FILE__, __LINE__; stop
		!!!
	end subroutine calcOcc_dadoxiz

	subroutine calcVEL_deferred(z,T,P)
		!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		!IMPLEMENTATION
		print*, 'DEFERRED @', __FILE__, __LINE__
	end subroutine calcVEL_deferred

	subroutine calcV_hyd(z,T,P)
		!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		!IMPLEMENTATION
		integer :: i
		real(8) :: g_w_ratio
		
		g_w_ratio = 0.d0
		do i = 1, z%nguest
			g_w_ratio = g_w_ratio + z%x(z%guestindex(i))
		end do !i = 1, z%nguest
		g_w_ratio = g_w_ratio/z%x(1)
		!z%V [m3/mol_phase]:
		z%V=(((1.d0/z%V_EL)+(g_w_ratio/z%V_EL))**(-1))
!		print*, z%V, 'hidrato'
	end subroutine calcV_hyd
	
	subroutine calcFugwH_deferred(z,T,P)
		!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		!IMPLEMENT
		print*, 'DEFERRED @', __FILE__, __LINE__
	end subroutine calcFugwH_deferred
	
	subroutine calcvdWPfunc(z,T,P)
		!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		integer :: i, j
		!IMPLEMENT
		real(8) :: arg, sum1
		arg = 0.d0
		do j = 1, 2
			sum1 = 0.d0
			do i = 1, z%nguest
				sum1 = sum1 + z%Occ(i,j)
			end do !i = 1, z%nguest
!			print*, sum1, 'sum1'
			arg = arg + z%ni(j)*dlog(1-sum1)
		end do !j = 1, 2
		z%vdW_P_func = dexp(arg)
!		print*, z%Occ(:,:), z%ni(:), arg, dexp(arg), 'arg'; pause
	end subroutine calcvdWPfunc
	
	subroutine calc_Phi_hyd(z,T,P)
		!declarations
		class(hydrate_model_c) :: z
		real(8) :: T, P
		!definitions
		real(8) :: Phi(ncomp)
		integer :: i
		!IMPLEMENT
		call z%calcClang(T,P)
		call z%calcOcc_dadoFUG(T,P) !calcula o x correspondente
		call z%calcVEL(T,P)
		call z%calcV(T,P)			!chama depois do occ e do VEL q ele calcula o V ocupado
		call z%calcvdWPfunc(T,P)
		call z%calcFugWH(T,P)
		!conclusion
		Phi(1) = z%fug_w_H/(z%x(1)*P)
		
		!CONFERIR O USO DO GUEST INDEX AQUI
		
		do i = 2, ncomp
			Phi(i) = dexp(20.d0)
			!para os componentes não presentes:
			!phi 'muito grande' mas não infinity
		end do
		
		do i = 1, z%nguest
		
			if( z%x(z%guestindex(i)) > 1.d-20 ) then !cuidado para não dividir por zero
				Phi(z%guestindex(i)) = z%fug_g_H(i)/(z%x(z%guestindex(i))*P)
			else
				!se o comp for guest mas não aparecer na fase
				!ele ganha um x quase zero a partir do phi 'muito grande' mas não infinity
!				z%x(z%guestindex(i)) = 1d-22
!				Phi(z%guestindex(i)) = z%fug_g_H(i)/(1d-22*P)
!				z%x(z%guestindex(i)) = Phi(z%guestindex(i))*P/z%fug_g_H(i)
!				z%x(z%guestindex(i)) = Phi(z%guestindex(i))*P/max(z%fug_g_H(i),1.d-20)
!				print*, z%x(z%guestindex(i)), Phi(z%guestindex(i)), P, z%fug_g_H(i), __FILE__, __LINE__
!				print*, dlog(z%fug_g_H(i)/(1d-12*P)); pause
			end if
		end do !i = 1, nguest
		
		
!			if( z%x(i+1) > 1d-22 ) then
!				Phi(i+1) = z%fug_g_H(i)/(z%x(i+1)*P)
!			else
!				z%x(i+1) = 1d-22
!				Phi(i+1) = z%fug_g_H(i)/(1d-22*P)
!!				print*, dlog(z%fug_g_H(i)/(1d-12*P)); pause
!			end if
!		end do !i = 1, nguest
!		do i = 1+z%nguest+1, ncomp
!			Phi(i) = dexp(50.d0)
!		end do !i = 1+z%nguest+1, ncomp
		z%phi(:) = phi(:)
	end subroutine calc_Phi_hyd
end module hydrates_mod

