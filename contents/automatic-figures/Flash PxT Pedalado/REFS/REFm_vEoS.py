import numpy as np

#objetivo - código simples, sem polimorfismo, implementação de uma única EoS, toma Peng e Robinson como base, inclui modificações, fazer consistency checks nessa rotina (tipos e valores)
class c_vEoS(): #Peng Robinson
  hardcodeparam=1
  def __init__(self,component_list_by_name,iglobal,iprops,pure_substance_DATA): #roda uma vez para carregar as propriedades por componentes
    #Array dimensioning info
    self.ncomp = np.size(component_list_by_name)
    #vEoS specific parameters
    self.sigma = 1.0 + np.sqrt(2.)
    self.epsilon = 1.0 - np.sqrt(2.)
    self.ac = np.zeros(self.ncomp)
    self.bc = np.zeros(self.ncomp)
    self.k = np.zeros([self.ncomp,self.ncomp])
    #Extracted pure component properties
    self.Tc = np.zeros(self.ncomp)
    self.Pc = np.zeros(self.ncomp)
    self.Vc = np.zeros(self.ncomp)
    self.PitzerAf = np.zeros(self.ncomp)
    for i, name_i in enumerate(component_list_by_name):
      self.Tc[i] = pure_substance_DATA[iglobal[component_list_by_name[i]],iprops["Tc"]]
      self.Pc[i] = pure_substance_DATA[iglobal[component_list_by_name[i]],iprops["Pc"]]
      self.Vc[i] = pure_substance_DATA[iglobal[component_list_by_name[i]],iprops["Vc"]]
      self.PitzerAf[i] = pure_substance_DATA[iglobal[component_list_by_name[i]],iprops["PitzerAf"]]


    
    
    
    
    
    
#call like Fluid=vEoS(comp_list(comp_names_list))

def akbPR(Tc,Pc,w):
  '''do i = 1, ncomp
  if( trim(comp_list(i)%comp%name_) == 'water' ) then
  z%ai(i) = z%aparam*0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
  z%bi(i) = z%bParam*0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc) !; print*, 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc), 'b'; pause
  else
  z%ai(i) = 0.45724d0*((R)**2)*((comp_list(i)%comp%Tc)**2)/(comp_list(i)%comp%Pc)
  z%bi(i) = 0.07780d0*(R)*(comp_list(i)%comp%Tc)/(comp_list(i)%comp%Pc)
  end if
  z%k(i) = 0.37464d0 + 1.54226d0*comp_list(i)%comp%w-0.26992d0*(comp_list(i)%comp%w)**2
  end do'''
  return

#water
##faw
#faw = 0.960392
##fbw
#fbw = 0.848496 #!fugacidade a 272.d0K, 100.d5Pa dá 607.856Pa
##palphaw
#palphaw = np.array([0.816384e0, 0.128223e1, 0.538174e-3])

def phiTPx(T,P,x,a,b,k,n):
    print("men@work")
    return 0 

def phi(T,V,P,x,n):
    print("men@work")
    return 0 

def mix_Theta_vdw1f(x,Theta,k,n):
  ThetaM = 0
  for i in range(n):
    for j in range(n):
      ThetaM += x[i]*x[j]*sqrt(Theta[i]*Theta[j])*(1.-k[i,j])
  return ThetaM


####################tradução

'''

function calc_Thetai_bar_vdw1f(z) !,x,Thetai)
  class(MR_vdw1f) :: z
!  real(8) :: x(ncomp), Thetai(ncomp), ,
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

  real(8) function mix_b_vdw1f_symmetric(z) !,x,bi)
  class(MR_vdw1f_symmetric) :: z
  real(8) :: b
!  real(8) :: x(ncomp), bi(ncomp)
  integer :: i
  b = 0.d0
  do i = 1, ncomp
    b = b + z%x(i)*z%pur_cub_eos_par_obj%bi(i)
  end do
  mix_b_vdw1f_symmetric = b
  end  function mix_b_vdw1f_symmetric

  function calc_bi_bar_vdw1f_symmetric(z) !,x,bi)
  class(MR_vdw1f_symmetric) :: z
!  real(8) :: x(ncomp), bi(ncomp), ,
  real(8) :: calc_bi_bar_vdw1f_symmetric(ncomp), bi_bar(ncomp)
  integer :: i, j
  bi_bar(:) = 0.d0
  do i = 1,ncomp
    bi_bar(i)=z%pur_cub_eos_par_obj%bi(i)
  end do
  calc_bi_bar_vdw1f_symmetric = bi_bar
  end function calc_bi_bar_vdw1f_symmetric



subroutine calc_phi_vdw1f_symmetric(z,T,P) !,x,T,P,V,phi)
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
!    print*, 'passo2'
    z%phi(:) = dexp(lnphi(:))
    !
!    print*, __FILE__, __LINE__; stop
  end subroutine calc_phi_vdw1f_symmetric


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
!      print*, z%alphai(i), 'it is water'; pause
      else
        z%alphai(i) = (1.d0 +z%k(i)*(1.d0-dsqrt(T/comp_list(i)%comp%Tc)))**2
      end if
      z%thetai(i) = z%ai(i)*z%alphai(i)
    end do
    !
  end subroutine calc_T_dep_PR

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
!          print*, i, __FILE__, __LINE__
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
!    print*, T, z%x(1), z%phase_cond, Z_1, 'Z_1'
!    calcular demais raízes para realizar teste de estabilidade e retornar a mais estável ou pré-definida L/V
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

!    print*, z%V, (z%V*P)/(R*T), 'v e z'
  end subroutine calc_Volume
'''
