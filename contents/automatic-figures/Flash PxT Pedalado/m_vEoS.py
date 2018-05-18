#this is a compact implementation of an PR eos
#look for some eos notebook to see commented development leading to this
import numpy as np
from scipy.constants import R
#objetivo - código simples, sem polimorfismo, implementação de uma única EoS, toma Peng e Robinson como base, inclui modificações, fazer consistency checks nessa rotina (tipos e valores)
class c_vEoS(): #Peng Robinson
  def __init__(self,ncomp,Tc,Pc,acentric,k): #roda uma vez para carregar as propriedades por componentes
    #Array dimensioning info
    self.ncomp = ncomp
    #vEoS specific parameters
    self.sigma = 1.0 + np.sqrt(2.)
    self.epsilon = 1.0 - np.sqrt(2.)
    self.ac = np.zeros(self.ncomp)
    self.bc = np.zeros(self.ncomp)
    self.k = np.zeros([self.ncomp,self.ncomp])
    
    #Extracted pure component properties
    self.Tc = np.zeros(self.ncomp) #needed at every alpha updating
    self.Pc = np.zeros(self.ncomp) #not really needed after initialization
    self.acentric = np.zeros(self.ncomp) #needed at every alpha updating
    
    self.kPR = np.zeros(self.ncomp) #needed at every alpha updating
    
    for i in range(self.ncomp):
      self.Tc[i]           = Tc[i]
      self.Pc[i]           = Pc[i]
      self.acentric[i]     = acentric[i]
    
    for i in range(self.ncomp):
      self.ac[i]           = 0.45724*((R)**2)*((self.Tc[i])**2)/(self.Pc[i])
      self.bc[i]           = 0.07780*(R)*(self.Tc[i])/(Pc[i])

      for j in range(self.ncomp):
        self.k[i,j]        = k[i,j]
        
      self.kPR[i]          = 0.37464 + 1.54226*acentric[i]-0.26992*(acentric[i])**2
        
    return #NoneTypeObj

  def Pressure(self,T,V,x):
    bm=self.f_bmix(x)
    Aalpham,Aalpha=self.f_Aalphamix(T,x)
    P = (R*T)/(V-bm) - Aalpham/(V**2 + 2*bm*V - bm**2) #sigma & epsilon hardcoded here
    return P

  def f_Aalpha(self,T,x):
    alpha=np.zeros(self.ncomp)
    Aalpha=np.zeros(self.ncomp)
    for i in range(self.ncomp):
      alpha[i] = (1. +self.kPR[i]*(1.-np.sqrt(T/self.Tc[i])))**2
      Aalpha[i] = self.ac[i]*alpha[i]
    return Aalpha

  def f_bmix(self,x):
    bm = 0.
    for i in range(self.ncomp):
      bm += x[i]*self.bc[i]
    return bm
    
  def f_Aalphamix(self,T,x):
    Aalpha=self.f_Aalpha(T,x)
    Aalpham = 0.
    for i in range(self.ncomp):
      for j in range(self.ncomp):
        Aalpham += x[i]*x[j]*np.sqrt(Aalpha[i]*Aalpha[j])*(1.-self.k[i,j])
    return Aalpham, Aalpha
   
  def f_dbdn(self,x):
    bm=self.f_bmix(x)
    dbdn=np.zeros(self.ncomp)
    for i in range(self.ncomp):
      dbdn[i]=self.bc[i]
    return dbdn, bm
    
  def f_dAalphadn(self,T,x):
    Aalpham, Aalpha = self.f_Aalphamix(T,x)
    dAalphadn = np.zeros(self.ncomp)
    sum1 = 0.
    for i in range(self.ncomp):
      sum1 = 0.
      for j in range(self.ncomp):
        sum1 += x[j]*np.sqrt(Aalpha[j])*(1.-self.k[i,j])
      dAalphadn[i]=np.sqrt(Aalpha[i])*sum1
    return dAalphadn, Aalpham

  def Volume(self,T,P,x):
  # T em unidade K
  # P em unidade Pa
  # x array normalizado

    bm=self.f_bmix(x)
    Aalpham,_=self.f_Aalphamix(T,x)

    c3 = 1.                                # Coeficiente para V^3 para EoS PR
    c2 = bm - R*T/P                        # Coeficiente para V^2 para EoS PR
    c1 = Aalpham/P - 3.*(bm**2) - 2.*bm*R*T/P   # Coeficiente para V^1 para EoS PR
    c0 = (R*T*bm**2)/P + bm**3 - Aalpham*bm/P  # Termo independente para EoS PR
    
    Vs=np.roots([c3,c2,c1,c0])
    Vs[np.logical_not(np.isreal(Vs))]=0.
    Vs=np.real(Vs)
    return np.array([np.nanmin(Vs[Vs>bm]),np.nanmax(Vs[Vs>bm])])

  #phase equilibrium common
  def fugacity_coeff(self,T,V,x): #for a vdw1f mixrule cubic eos with sigma!=epsilon
    P=self.Pressure(T,V,x)
    dbdn,bm = self.f_dbdn(x)
    dAalphadn, Aalpham = self.f_dAalphadn(T,x)
    qsi = (1./(bm*(self.epsilon-self.sigma)))*np.log((V+self.epsilon*bm)/(V+self.sigma*bm))
    lnPhi = np.zeros(self.ncomp)
    for i in range(self.ncomp):
      lnPhi[i] = ( #multiline
        (dbdn[i]/bm)*((P*V)/(R*T)-1) #&
        -np.log(P*(V-bm)/(R*T)) #&
        -(Aalpham/(R*T))*qsi*((2.*dAalphadn[i]/Aalpham) #&
        -(dbdn[i]/bm))
                 )#done
    phi = np.exp(lnPhi)
    return phi

  def f_dAalphadT(self,T,x):
    dalphadT=np.zeros(self.ncomp)
    dAalphadT=np.zeros(self.ncomp)
    for i in range(self.ncomp):
      #alpha[i] = (1. + self.kPR[i]*(1.-np.sqrt(T/self.Tc[i])))**2
      dalphadT[i] = ( #multiline
        2.*(1. + self.kPR[i]*(1.-np.sqrt(T/self.Tc[i])))*
        self.kPR[i]*(-1.)*(
        1./(2.*(np.sqrt(T/self.Tc[i])))
        )*(1./self.Tc[i])
      )
      dAalphadT[i] = self.ac[i]*dalphadT[i]
    return dAalphadT

  def f_dAalphamdT(self,T,x):
    _, Aalpha=self.f_Aalphamix(T,x)
    dAalphadT=self.f_dAalphadT(T,x)

    dAalphamdT = 0.
    for i in range(self.ncomp):
      for j in range(self.ncomp):
        dAalphamdT += ( #multiline
          x[i]*x[j]*
          (1./(2.*np.sqrt(Aalpha[i]*Aalpha[j])))*
          (Aalpha[i]*dAalphadT[j]+dAalphadT[i]*Aalpha[j])*
          (1.-self.k[i,j])
        )
    #numerical
    #AalphammaisT,_=self.f_Aalphamix(T+1e-3,x)
    #AalphammenosT,_=self.f_Aalphamix(T-1e-3,x)
    #dAalphamdT=(AalphammaisT-AalphammenosT)/(2.*1e-3)
    return dAalphamdT

  #other spec flashes
  def f_H_res(self,T,V,x):
    P=self.Pressure(T,V,x)
    bm=self.f_bmix(x)
    dAalphamdT = self.f_dAalphamdT(T,x)
    Aalpham, _ = self.f_Aalphamix(T,x)
    qsi = (1./(bm*(self.epsilon-self.sigma)))*np.log((V+self.epsilon*bm)/(V+self.sigma*bm)) #again
    H_res = P*V - R*T + (T*dAalphamdT-Aalpham)*qsi
    return H_res

  def f_S_res(self,T,V,x):
    P=self.Pressure(T,V,x)
    bm=self.f_bmix(x)
    dAalphamdT = self.f_dAalphamdT(T,x)
    Aalpham, Aalpha = self.f_Aalphamix(T,x)
    qsi = (1./(bm*(self.epsilon-self.sigma)))*np.log((V+self.epsilon*bm)/(V+self.sigma*bm)) #again
    S_res = R*np.log((P*(V-bm))/(R*T)) + dAalphamdT*qsi
    return S_res

  #optimization issues
  #i will recalc bmix and amix when calling phi and when calling v
  #optimizing this would require that either bmix were public and required prior to calc V
  # or that there were both a public calcv(t,p,x) and a private calv(t,p,x,am,bm) and either class variable bm and am wih status checking at every call or combo calls calcv(tpx) calcphi(tpx) calcv_and_phi(tpx) combinatorially.
  
