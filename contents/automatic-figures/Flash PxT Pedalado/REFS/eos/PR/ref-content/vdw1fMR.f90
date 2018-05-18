!esse .f90 deverá conter os modules MR_vdw1f, MR_vdw1f_symmetric e MR_pure

module MR_vdw1f_mod !modulo das diversas regras de mistura para ede cubicas
	use comp_list_mod !precisa usar o modulo lista de componentes, (e, assim, o modulo pure prop é usado indiretamente)
	use cub_eos_mod
	implicit none
	type, extends(cub_eos_c) :: MR_vdw1f !regra de van der Waals para 1 fluido - vdw1f (regra de mistura quadrática para ambos os parametros e com regra de combinação clássica
		!(REF:Kontogeorgis e Folas, 2010, p. 43.)
		!válida para misturas aleatórias
		!a regra de vdw1f tem 2 matrizes de parametros, Kij e lij
		real(8), allocatable :: kij(:,:), lij(:,:)
		contains
!		procedure, nopass :: init_vdw1f_fil
		procedure :: mix_Theta => mix_Theta_vdw1f
!		procedure :: mix_b => mix_b_vdw1f
		procedure :: calc_Thetai_bar => calc_Thetai_bar_vdw1f
!		procedure :: calc_bi_bar => calc_bi_bar_vdw1f !***falta
!		procedure :: calc_phi => calc_phi_vdw1f !***falta
		procedure :: debug_phase_model => debug_MR_vdw1f
	end type MR_vdw1f
!	
	interface MR_vdw1f_i
		module procedure :: init_MR_vdw1f
	end interface
!	
	!agora que já declaramos todas as variáveis e definimos todas as classes do módulo vamos implementar cada função, seja para inicializações, seja para cálculos
	contains
!	
	subroutine debug_MR_vdw1f(z)
		class(MR_vdw1f) :: z
		call z%cub_eos_c%debug_phase_model
		print*, ' '
		print*, 'debug_MR_vdw1f', z%class_name
		print*, 'z%kij(:,:)', z%kij(:,:)
		print*, 'z%lij(:,:)', z%lij(:,:)
	end subroutine debug_MR_vdw1f
	
	function init_MR_vdw1f(filename)
		!declarations
		type(MR_vdw1f) :: z, init_MR_vdw1f
		character(100) :: filename
		allocate(z%kij(ncomp,ncomp), z%lij(ncomp,ncomp))
		z%kij(ncomp,ncomp) = 0.d0
		z%lij(ncomp,ncomp) = 0.d0
		init_MR_vdw1f = z
		print*, 'men at work ', __FILE__, __LINE__
		stop
	end function init_MR_vdw1f
	
!IMPLEMENTAÇÃO DAS EXPRESSÕES DE CÁLCULO DAS REGRAS DE MISTURA
!
	real(8) function mix_Theta_vdw1f(z) !,x,Thetai)
	class(MR_vdw1f) :: z
	real(8) :: a
!	real(8) :: x(ncomp), Thetai(ncomp)
	integer :: i, j
	a = 0.d0
	do i = 1, ncomp
		do j = 1, ncomp
			a = a + z%x(i)*z%x(j)*dsqrt(z%pur_cub_eos_par_obj%Thetai(i)*z%pur_cub_eos_par_obj%Thetai(j))*(1.d0-z%kij(i,j))
		end do
	end do
	mix_Theta_vdw1f = a
	end  function mix_Theta_vdw1f
		
	function calc_Thetai_bar_vdw1f(z) !,x,Thetai)
	class(MR_vdw1f) :: z
!	real(8) :: x(ncomp), Thetai(ncomp), ,
	real(8) :: calc_Thetai_bar_vdw1f(ncomp), Thetai_bar(ncomp)
	real(8) :: sum1
	integer :: i, j
	Thetai_bar(:) = 0.d0
	sum1 = 0.d0
	do i = 1,ncomp
		sum1 = 0.d0
		do j=1,ncomp
			sum1 = sum1+z%x(j)*dsqrt(z%pur_cub_eos_par_obj%Thetai(j))*(1-z%kij(i,j))
		end do
		Thetai_bar(i)=dsqrt(z%pur_cub_eos_par_obj%Thetai(i))*sum1
	end do
	calc_Thetai_bar_vdw1f = Thetai_bar
	end function calc_Thetai_bar_vdw1f
!	
!	Expressão para cálculo de phi
!	***
end module MR_vdw1f_mod

module MR_vdw1f_symmetric_mod
	use MR_vdw1f_mod
	implicit none
	type, extends(MR_vdw1f) :: MR_vdw1f_symmetric !essa extensão é na verdade uma simplificação, então não vão ser reutilizadas muitas funções da classe mãe como preambulo
		!regra de mistura de vdw1f com lij igual a zero e expressões simplificadas analiticamente
		!válida para misturas aproximadamente simétricas
		contains
!		procedure :: mix_Theta  => mix_Theta_vdw1f !a mesma que a das misturas assimétricas
		procedure :: mix_b => mix_b_vdw1f_symmetric
!		procedure :: calc_Thetai_bar => calc_Thetai_bar_vdw1f !a mesma que a das misturas assimétricas
		procedure :: calc_bi_bar => calc_bi_bar_vdw1f_symmetric
		procedure :: calc_phi => calc_phi_vdw1f_symmetric
		procedure :: debug_phase_model => debug_MR_vdw1f_symmetric
	end type MR_vdw1f_symmetric
	!
	interface MR_vdw1f_symmetric_i
		!a primeira opção será implementada para ler arquivo, a segunda para vallor padrão igual a zero, a terceira para usar uma matriz pré-preparada
		module procedure :: init_MR_vdw1f_symmetric_def
		module procedure :: init_MR_vdw1f_symmetric_fnum_purobj
	end interface
	
	contains
	
	subroutine debug_MR_vdw1f_symmetric(z)
		class(MR_vdw1f_symmetric) :: z
		integer :: i
		call z%cub_eos_c%debug_phase_model
		print*, ' '
		print*, 'debug_MR_vdw1f_symmetric', z%class_name
		do i = 1, ncomp
			print*, 'z%kij(',i,',:)', z%kij(i,:)
		end do !i = 1, ncomp
	end subroutine debug_MR_vdw1f_symmetric
	
	function init_MR_vdw1f_symmetric_fnum_purobj(fnum_in,purobj_in)
		type(MR_vdw1f_symmetric) :: z, init_MR_vdw1f_symmetric_fnum_purobj
		!arguments
		integer :: fnum_in
		character(100) :: mixingrule_par_fname
		class(pur_cub_eos_par_c) :: purobj_in
		!declarations
!		print*, 'men at work', __FILE__, __LINE__
!		print*, 'parte pura alocada, proximo passo é alocar os vetores da mistura, ler os dados de inicilização no arquivo de simulação'
!		print*, 'e ler os parametros kij na base de parametros'
!		print*, 'para isso, preciso ver primeiro a classe mãe cub_eos, e a classe mãe phase_model'
!		print*, 'obs:quem abre o arquivo de simulação da fase vai ser o alocador: (o main ou phase list), e daí passa só o file number'
!		print*, 'lembrar de dar o close nesse próprio alocador'
		!preamble
		z%cub_eos_c = cub_eos_i(fnum_in,purobj_in)
!		print*, 'phasemodel, cubeos verificados, estamos de volta aqui'
		!implementation
		z%class_name = 'MR_vdw1f_sym'
		allocate(z%kij(ncomp,ncomp))
		z%kij = 0.d0
		mixingrule_par_fname = 'input/parameters/'//trim(z%class_name)//'/'//trim(z%pur_cub_eos_par_obj%class_name)//'/'//'kij.dat'
		call read_sym_binary_matrix(mixingrule_par_fname,z%kij)
		init_MR_vdw1f_symmetric_fnum_purobj = z
	end function init_MR_vdw1f_symmetric_fnum_purobj
	
!
	function init_MR_vdw1f_symmetric_def() !default
		!declarations
		type(MR_vdw1f_symmetric) :: z, init_MR_vdw1f_symmetric_def
		character(100) :: mixingrule_file
		
		!definitions
		allocate(z%kij(ncomp,ncomp))
		z%kij(:,:) = 0.d0
		init_MR_vdw1f_symmetric_def = z
		deallocate(z%kij)
!		print*, kij; pause
!		print*, 'mixingrule initialized'
		!pause
	end function init_MR_vdw1f_symmetric_def
	
	!IMPLEMENTAÇÃO DAS EXPRESSÕES DE CÁLCULO DAS REGRAS DE MISTURA
	
	real(8) function mix_b_vdw1f_symmetric(z) !,x,bi)
	class(MR_vdw1f_symmetric) :: z
	real(8) :: b
!	real(8) :: x(ncomp), bi(ncomp)
	integer :: i
	b = 0.d0
	do i = 1, ncomp
		b = b + z%x(i)*z%pur_cub_eos_par_obj%bi(i)
	end do
	mix_b_vdw1f_symmetric = b
	end  function mix_b_vdw1f_symmetric
	
	function calc_bi_bar_vdw1f_symmetric(z) !,x,bi)
	class(MR_vdw1f_symmetric) :: z
!	real(8) :: x(ncomp), bi(ncomp), ,
	real(8) :: calc_bi_bar_vdw1f_symmetric(ncomp), bi_bar(ncomp)
	integer :: i, j
	bi_bar(:) = 0.d0
	do i = 1,ncomp
		bi_bar(i)=z%pur_cub_eos_par_obj%bi(i)
	end do
	calc_bi_bar_vdw1f_symmetric = bi_bar
	end function calc_bi_bar_vdw1f_symmetric
	
!	Expressão para cálculo de phi
	subroutine calc_phi_vdw1f_symmetric(z,T,P) !,x,T,P,V,phi)
		class(MR_vdw1f_symmetric) :: z
		real(8) :: lnphi(ncomp)
		real(8) :: T, P
		integer :: i
		real(8) :: qsi
		!tipo preamble
		call z%pur_cub_eos_par_obj%calc_T_dep(T) !da classe da pura
		!essa classe não possui KIJ dependente de T então na precisa de z%calc_T_dep(T), só da parte pura mesmo.
		call z%calc_MixingParam() !rotina implementada na classe de mistura (filha) !depende apenas do x da propria classe e dos parametro de purto da classe objeto, logo nao precisa de argumentos.
		!aqui deu erro pois uma vez que eu subo para reciclar a função da classe mãe, qualquer chamada por aqui vai cair nas funções da mãe
		call z%calc_Volume(T,P) ! da classe cub_eos
		!agora retorna para a classe filha para o procedimento final
		!implementation
		if (z%pur_cub_eos_par_obj%sigma_ == z%pur_cub_eos_par_obj%epsilon_) then
			qsi = 1.d0/(z%V+z%pur_cub_eos_par_obj%sigma_*z%b_m)
		else
			qsi = (1.d0/(z%b_m*(z%pur_cub_eos_par_obj%epsilon_-z%pur_cub_eos_par_obj%sigma_)))*dlog((z%V+z%pur_cub_eos_par_obj%epsilon_*z%b_m)/(z%V+z%pur_cub_eos_par_obj%sigma_*z%b_m))
		!print*, qsi, 'qsi'; pause !OK
		end if
		lnPhi(:) = 0.d0
		do i = 1,ncomp
			lnPhi(i) = (z%bi_bar(i)/z%b_m)*((P*z%V)/(R*T)-1) &
				-dlog(P*(z%V-z%b_m)/(R*T)) &
				-(z%Theta_m/(R*T))*qsi &
				*((2.d0*z%thetai_bar(i)/z%Theta_m)- &
				(z%bi_bar(i)/z%b_m))
		end do
		!
!		print*, 'passo2'
		z%phi(:) = dexp(lnphi(:))
		!
!		print*, __FILE__, __LINE__; stop
	end subroutine calc_phi_vdw1f_symmetric
	
end module MR_vdw1f_symmetric_mod


