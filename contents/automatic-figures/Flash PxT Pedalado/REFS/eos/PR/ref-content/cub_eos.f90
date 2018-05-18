module cub_eos_mod
	use phase_model_mod !which does use comp_list_mod !which does use pureprop_mod
	use pur_cub_eos_par_mod
	implicit none
	type, extends(phase_model_c) :: cub_eos_c
!		character(5) :: phase_cond !L(liquida), V(vapor), ou S(mais estavel) para cúbicas, ...
		class(pur_cub_eos_par_c), pointer :: pur_cub_eos_par_obj !pra ter qualquer filho aqui tem que ser class, e para ser class tem que ser pointer[.]
!		// EoS parameters
		!// parameters Theta, b of mixture eos
		real(8) :: Theta_m, b_m
		!// calcPhi method
		!//partial molar parameters (?)
		real(8), allocatable :: bi_bar(:), Thetai_bar(:)
		contains
!		// constructor
!		// cub_eos specific methods
		procedure :: calc_MixingParam !chama as funções específicas das regras de mistura implementadas e alocadas
		procedure :: calc_Volume
!		// pure virtual methods implemented
!		virtual void calcPhi()=0;
		procedure :: calc_phi => calc_phi_cubic !não pode dar deferred se o tipo não for abstrato, e o tipo não pode ser abstrato se tiver alguma função normal (?)
		!mixingrule procedures
		procedure :: mix_Theta => mix_Theta_deferred
		procedure :: mix_b => mix_b_deferred
		procedure :: calc_bi_bar => calc_bi_bar_deferred
		procedure :: calc_Thetai_bar => calc_Thetai_bar_deferred
		!misc
		procedure :: debug_phase_model => debug_cub_eos
		procedure :: Plot_PvsV
	end type cub_eos_c
!
		interface cub_eos_i
			module procedure allocate_cub_eos
			module procedure init_cub_eos_fnum_purobj !usa arquivo por número, requer pur_cub_eos_par como argumento
			module procedure init_cub_eos_purobj !apenas aloca e inicializa com valores neutros, requer pur_cub_eos_par como argumento
		end interface
		
	contains
!	
	subroutine debug_cub_eos(z)
		class(cub_eos_c) :: z
		call z%phase_model_c%debug_phase_model
		print*, ''
		print*, 'debug_cub_eos ', z%class_name
		print*, 'z%Theta_m', z%Theta_m
		print*, 'z%b_m', z%b_m
		print*, 'z%Thetai_bar(:)', z%Thetai_bar(:)
		print*, 'z%bi_bar(:)', z%bi_bar(:)
	end subroutine debug_cub_eos
	
	function allocate_cub_eos()
		!declarations
		type(cub_eos_c) :: z, allocate_cub_eos
		!preamble
		z%phase_model_c = allocate_phase_model()
		!implementation
!		// adapting vector lengths to number of components
		z%class_name = 'cub_eos'
		allocate(z%thetai_bar(ncomp), z%bi_bar(ncomp))
		z%thetai_bar(:) = 0.d0
		z%bi_bar(:) = 0.d0
		z%Theta_m = 0.d0
		z%b_m = 0.d0
		allocate_cub_eos = z
		!destroy z ao final ou não usar z e usar direto a variavel resposta initcubeos?
	end function allocate_cub_eos
	
!
	function init_cub_eos_purobj(purobj_in)
		type(cub_eos_c) :: z, init_cub_eos_purobj
		!arguments
		class(pur_cub_eos_par_c) :: purobj_in
		!preamble and allocation
		z=allocate_cub_eos()
		z%phase_model_c = phase_model_i()
		!implementation
		z%class_name = 'cub_eos'
		allocate(z%pur_cub_eos_par_obj, source = purobj_in)
		init_cub_eos_purobj = z !*** isso vai funcionar ou precisarei fazer tudo pointer e dar allocate init source = z ???
!		call init_cub_eos_fnum_purobj%pur_cub_eos_par_obj%debug_pur_ceos_par !funcionou
	end function init_cub_eos_purobj
!
	function init_cub_eos_fnum_purobj(fnum_in,purobj_in)
		type(cub_eos_c) :: z, init_cub_eos_fnum_purobj
		!arguments
		integer :: fnum_in
		class(pur_cub_eos_par_c) :: purobj_in
		!preamble and allocation
		z=allocate_cub_eos()
		z%phase_model_c = phase_model_i(fnum_in)
		!implementation
		z%class_name = 'cub_eos'
		allocate(z%pur_cub_eos_par_obj, source = purobj_in)
		init_cub_eos_fnum_purobj = z !*** isso vai funcionar ou precisarei fazer tudo pointer e dar allocate init source = z ???
!		call init_cub_eos_fnum_purobj%pur_cub_eos_par_obj%debug_pur_ceos_par !funcionou
	end function init_cub_eos_fnum_purobj

!DEFERRED
	real(8) function mix_Theta_deferred(z) !,x,Thetai)
		class(cub_eos_c) :: z
		print*, 'DEFERRED', __FILE__, __LINE__; stop
	end  function mix_Theta_deferred
!	
	real(8) function mix_b_deferred(z) !,x,bi)
		class(cub_eos_c) :: z
		print*, 'DEFERRED', __FILE__, __LINE__; stop
	end  function mix_b_deferred
!	
	function calc_bi_bar_deferred(z) !,x,bi)
		real(8) :: calc_bi_bar_deferred(ncomp)
		class(cub_eos_c) :: z
		print*, 'DEFERRED', __FILE__, __LINE__; stop
	end function calc_bi_bar_deferred
	
	function calc_Thetai_bar_deferred(z) !,x,Thetai)
		real(8) :: calc_Thetai_bar_deferred(ncomp)
		class(cub_eos_c) :: z
		print*, 'DEFERRED', __FILE__, __LINE__; stop
	end function calc_Thetai_bar_deferred

!void cub_eos::calcMixingPar()
	subroutine calc_MixingParam(z) 
		!declarations
		class(cub_eos_c) :: z
		integer :: i
		!definitions
		!IMPLEMENT
		do i = 1, ncomp
			z%Theta_m = z%mix_Theta() !(z%x,z%thetai)
			z%b_m = z%mix_b() !(z%x,z%bi)
			z%Thetai_bar = z%calc_Thetai_bar() !(z%x,z%thetai)
			z%bi_bar = z%calc_bi_bar() !(z%x,z%bi)
		end do
	end subroutine calc_MixingParam
	
!void cub_eos::calcVolume()		
	subroutine calc_Volume(z,T,P)
!declarations
		class(cub_eos_c) :: z
		real(8) :: T, P
!IMPLEMENT
		real(8) :: alpha_cub, beta_cub, gama_cub ! parametros do polinomio de 3o grau
		real(8) :: ni_quad, phi_quad, delta_bhaskara !parametros do polinomio de 2o grau
		integer :: i, j, maxnit, frac_maxnit
		real(8) :: tol
		real(8) :: Z_1, Z_old, Z_2, Z_3 !raízes do compress_fac
		frac_maxnit = 500
		maxnit = 3*frac_maxnit
		tol = 1.d-6

!print*, z%ai, z%bi, z%alphai, z%thetai, 'debug1'

!print*, 'debug@calc_volume'
!print*, z%theta_m, z%b_m, 'theta_m, b_m'

	!alpha, beta e gama, 'dependem de T'
	alpha_cub = (z%pur_cub_eos_par_obj%sigma_+z%pur_cub_eos_par_obj%epsilon_-1.d0)*z%b_m*(P/(R*T))-1.d0
	beta_cub = (z%pur_cub_eos_par_obj%sigma_*z%pur_cub_eos_par_obj%epsilon_*(z%b_m**2) &
		-((R*T/P)+z%b_m)*(z%pur_cub_eos_par_obj%sigma_+z%pur_cub_eos_par_obj%epsilon_)*z%b_m &
		+(z%theta_m/P))*((P/(R*T))**2)
	gama_cub = (-((R*T/P)+(z%b_m))*(z%pur_cub_eos_par_obj%sigma_*z%pur_cub_eos_par_obj%epsilon_*(z%b_m**2))-(z%b_m*z%theta_m/P))*((P/(R*T))**3)
!print*, alpha_cub, beta_cub, gama_cub, 'polinomio em Z'

		Z_1 = z%compress_fac !estimativa inicial com o Z antigo da eos
		do i = 1,maxnit
			Z_old=Z_1
			Z_1 = (2.d0*(Z_old**3)+alpha_cub*(Z_old**2)-gama_cub)/(3.d0*(Z_old**2)+2.d0*alpha_cub*Z_old+beta_cub)
			if(i>2)then !minimum of 2 iterations desired
				if((abs(Z_1-Z_old))<tol) then
!					print*, i, __FILE__, __LINE__
					exit !do
				end if
				if(i==frac_maxnit) then
					Z_1 = 1.d0 !estimativa de gás ideal
					cycle !do
				end if
				if(i==2*frac_maxnit) then
					Z_1 = z%b_m*P/(R*T) !estimativa de pressão infinita
					cycle !do
				end if
				if(i==maxnit) then
					print*, 'too many steps in calc_Volume(T=', T,', P=', P,', x(:)=', z%x, 'cond', z%phase_cond,')'
					z%V = 1.d3*((R*T/P)+z%b_m)
					z%compress_fac = z%V*P/(R*T)
					stop
					return
				end if
			end if
		end do
!		print*, T, z%x(1), z%phase_cond, Z_1, 'Z_1'
!
!problemas observados:
!too many steps in calc_Volume(T=   276.48285085134154      , P=   1576332.8789174831      , x(:)=   3.5613883543010472E-003  0.99643861164569891      condL    )
!água/propano
!ao mudar a estimativa inicial de Z para 2, ele converge
!
!Após desenvolver a eos da água fria ele converge sem prob com zinicial igual a 1
!
!		calcular demais raízes para realizar teste de estabilidade e retornar a mais estável ou pré-definida L/V
		ni_quad = Z_1
		phi_quad = -(alpha_cub+ni_quad)/2.d0
		delta_bhaskara = phi_quad**2 + 2.d0*phi_quad*ni_quad - beta_cub
		if(delta_bhaskara<0) then !raiz unica
			if (Z_1>0.d0) then
				z%V = (R*T/P)*Z_1
				z%compress_fac = Z_1
				return
			else
				print*, 'no Z roots in calc_Volume(T=', T,', P=', P,', x(:)=', z%x, 'cond', z%phase_cond,')', __FILE__, __LINE__
				stop
			end if
		else
		Z_2 = phi_quad+dsqrt(delta_bhaskara)
		Z_3 = phi_quad-dsqrt(delta_bhaskara)
		end if
		if (Z_1>0.d0) then
			if (Z_2>0.d0) then
				if (Z_3>0.d0) then
						if(z%phase_cond=='V') then
							z%V=(R*T/P)*max(Z_1,Z_2,Z_3)
							z%compress_fac = max(Z_1,Z_2,Z_3)
						elseif(z%phase_cond=='L') then
							z%V=(R*T/P)*min(Z_1,Z_2,Z_3)
							z%compress_fac = min(Z_1,Z_2,Z_3)
						else
							print*, 'unexpected phase_cond argument in ', __FILE__, __LINE__, z%phase_cond; stop
						end if
				else
					z%V=(R*T/P)*Z_1 !SE V_3 FOR NEGATIVO E V_2 POSITIVO, NÃO POSSO PEGAR V_3 POR SER NEGATIVO,
										!MAS TB NÃO POSSO PEGAR V_2 POR SER INTERMEDIÁRIA
					z%compress_fac = Z_1
				end if
			end if
		else
			print*, 'no Z roots in calc_Volume(T=', T,', P=', P,', x(:)=', z%x, 'cond', z%phase_cond,')', __FILE__, __LINE__
			stop
		end if

!		print*, z%V, (z%V*P)/(R*T), 'v e z'
	end subroutine calc_Volume
!!!!!//void cub_eos::calcPhi()
	subroutine calc_Phi_cubic(z,T,P)
		!declarations
		class(cub_eos_c) :: z
		real(8) :: T, P
		!DEFERRED por não saber usar regras de mistura ainda
		print*, 'DEFERRED', __FILE__, __LINE__; stop
!		real(8) :: lnphi(ncomp)
!		real(8) :: qsi
!		integer :: i
		!definitions
		!IMPLEMENT
		!
!		call z%pur_cub_eos_par_obj%calc_T_dep(T) !da classe da pura
!		!essa classe não possui KIJ dependente de T então na precisa de z%calc_T_dep(T), só da parte pura mesmo.
!		call z%calc_MixingParam() !rotina implementada na classe de mistura (filha) !depende apenas do x da propria classe e dos parametro de purto da classe objeto, logo nao precisa de argumentos.
!		!aqui deu erro pois uma vez que eu subo para reciclar a função da classe mãe, qualquer chamada por aqui vai cair nas funções da mãe
!		call z%calc_Volume(T,P) ! dessa mesma (cub_eos)
!		!agora retorna para a classe filha para o procedimento final
		!
!		print*, __FILE__, __LINE__; stop
	end subroutine calc_Phi_cubic
	
	subroutine  Plot_PvsV(z,T,fnum_in)
		!declarations
		class(cub_eos_c) :: z
		real(8) :: T
		integer, intent(in) :: fnum_in !***gambiarra-improviso
		integer :: i
		real(8) :: V, P !temporárias para o gráfico
		!definitions
		!IMPLEMENT
		!
		call z%calc_MixingParam()
		V = z%b_m*1.01d0
		write(fnum_in,*) 'phase_cond        ', '(', '   z%x(           1 )', ('  / ','        z%x(',i,')',i=2,ncomp), ' )'
		write(fnum_in,*) z%phase_cond,'        ', '(', z%x(1), (' / ',z%x(i),i=2,ncomp), ')'
		write(fnum_in,*) '                        V                    P'
		do i = 1,1000
			V = V*1.02d0
			P = R*T/(V-z%b_m) - z%Theta_m/((V+z%pur_cub_eos_par_obj%sigma_*z%b_m)*(V+z%pur_cub_eos_par_obj%epsilon_*z%b_m))
			write(fnum_in,*) V, P
		end do
	end subroutine Plot_PvsV
	!
end module cub_eos_mod

module ceos_pure_mod !talvez esse modulo pudesse sumir se eu desse essas rotinas de cálculo para o próprio pur_cub_eos_par, e chamaria ele de pur_cub_eos apenas
					!mas lembrando que, o pur_cub_eos não sendo extensão de phase model, eu não poderia usar pra modelar fase do flash por exemplo.
					!esse poderi ficar simplificado, fazendo	ceos_pure_obj%phi = pur_cub_eos_par%calc_phi(T,P,x)
!																pur_cub_eos_based_phase_model%phi = pur_cub_eos_model%calc_phi(T,P,x)
	use cub_eos_mod
	implicit none
	type, extends(cub_eos_c) :: ceos_pure_c !solução ideal
		integer :: pure_index !índice de qual substancia na lista de componentes será a substÂncia modelada como pura
		contains
		!mixingrule procedures
		procedure :: mix_Theta => mix_Theta_pure
		procedure :: mix_b => mix_b_pure
		procedure :: calc_Thetai_bar => calc_Thetai_bar_pure
		procedure :: calc_bi_bar => calc_bi_bar_pure
		
		procedure :: calc_phi => calc_phi_pure

		procedure :: set_ceos_pure_index
		!misc
		procedure :: debug_phase_model => debug_ceos_pure
	end type ceos_pure_c

!INTERFACE
	interface ceos_pure_i
		module procedure :: init_ceos_pure !allocate_ceos_pure desnecessário pois não tem nada pra alocar
	end interface
	
	contains
!INIT
	function init_ceos_pure(purobj_in)
		!arguments
		type(ceos_pure_c) :: init_ceos_pure, z
		character(5) :: phase_cond_in
		class(pur_cub_eos_par_c) :: purobj_in
		!preamble
		z%cub_eos_c = cub_eos_i(purobj_in) !***ou manda arquivo pela lista de args ou usa uma função diferente, com a phase cond apenas, por exemplo
		!implementation
		z%class_name = 'ceos_pure'
		z%pure_index = 0 !não inicializado
!		print*, 'men at work', __FILE__, __LINE__; stop
		init_ceos_pure = z
	end function init_ceos_pure
	
	subroutine set_ceos_pure_index(z,pure_index_in)
		class(ceos_pure_c) :: z
		integer :: pure_index_in
		z%pure_index = pure_index_in
		z%x(:) = 0.d0
		z%x(z%pure_index) = 1.d0
	end subroutine set_ceos_pure_index
	
!IMPLEMENTAÇÃO DAS EXPRESSÕES DE CÁLCULO DAS REGRAS DE MISTURA

	real(8) function mix_Theta_pure(z)
		class(ceos_pure_c) :: z
		mix_Theta_pure = z%pur_cub_eos_par_obj%Thetai(z%pure_index)
	end function mix_Theta_pure
	
	function calc_Thetai_bar_pure(z) !,x,Thetai)
		class(ceos_pure_c) :: z
		real(8) :: calc_Thetai_bar_pure(ncomp)
		calc_Thetai_bar_pure(:) = 0.d0
		calc_Thetai_bar_pure(z%pure_index) = z%pur_cub_eos_par_obj%Thetai(z%pure_index)
	end function calc_Thetai_bar_pure
	
	real(8) function mix_b_pure(z) !,x,bi)
		class(ceos_pure_c) :: z
		mix_b_pure = z%pur_cub_eos_par_obj%bi(z%pure_index)
	end  function mix_b_pure
	
	function calc_bi_bar_pure(z) !,x,bi)
		class(ceos_pure_c) :: z
		real(8) :: calc_bi_bar_pure(ncomp)
		calc_bi_bar_pure(:) = 0.d0
		calc_bi_bar_pure(z%pure_index) = z%pur_cub_eos_par_obj%bi(z%pure_index)
	end function calc_bi_bar_pure

!calculo de phi
	subroutine calc_phi_pure(z,T,P)
		class(ceos_pure_c) :: z
		real(8) :: lnphi(ncomp)
		real(8) :: T, P
		integer :: i
		real(8) :: qsi
		!tipo preamble
		
		call z%pur_cub_eos_par_obj%calc_T_dep(T) !da classe da pura
		!essa classe não possui KIJ dependente de T então na precisa de z%calc_T_dep(T), só da parte pura mesmo.
!		call z%calc_MixingParam() !rotina implementada na classe de mistura (filha) !depende apenas do x da propria classe e dos parametro de purto da classe objeto, logo nao precisa de argumentos.
		!aqui deu erro pois uma vez que eu subo para reciclar a função da classe mãe, qualquer chamada por aqui vai cair nas funções da mãe
			!necessária no pure apenas para gerar os parametros globais para reciclar a função de calc volume:
		z%Theta_m = z%mix_Theta() !(z%x,z%thetai)
		z%b_m = z%mix_b() !(z%x,z%bi)
		!não chamo o thetaibar e o bibar pois não usarei-os nos cálculos a seguir.
		call z%calc_Volume(T,P) ! da classe cub_eos
		!agora retorna para a classe filha para o procedimento final
		!implementation
		!simplificação analitica de que B e theta parcial molares são iguais aos dos puros e iguais aos da mistura aplicada:
		if (z%pur_cub_eos_par_obj%sigma_ == z%pur_cub_eos_par_obj%epsilon_) then
			qsi = 1.d0/(z%V+z%pur_cub_eos_par_obj%sigma_*z%b_m)
		else
			qsi = (1.d0/(z%b_m*(z%pur_cub_eos_par_obj%epsilon_-z%pur_cub_eos_par_obj%sigma_)))*dlog((z%V+z%pur_cub_eos_par_obj%epsilon_*z%b_m)/(z%V+z%pur_cub_eos_par_obj%sigma_*z%b_m))
		!print*, qsi, 'qsi'; pause !OK
		end if
		lnPhi(:) = 20.d0 !arbitrariamente grande
		lnPhi(z%pure_index) = ((P*z%V)/(R*T)-1.d0) &
				-dlog(P*(z%V-z%b_m)/(R*T)) &
				-(z%Theta_m/(R*T))*qsi
		z%phi(:) = dexp(lnphi(:))
!		print*, __FILE__, __LINE__; stop
	end subroutine calc_phi_pure
!MISC
	subroutine debug_ceos_pure(z)
		class(ceos_pure_c) :: z
		call z%cub_eos_c%debug_phase_model
		print*, ' '
		print*, 'debug_ceos_pure ', trim(z%class_name)
		print*, 'z%pure_index', z%pure_index
	end subroutine debug_ceos_pure
	
end module ceos_pure_mod
