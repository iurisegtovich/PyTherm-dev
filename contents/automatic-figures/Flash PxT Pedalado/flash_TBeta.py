import numpy as np

from m_vEoS import c_vEoS

#def load_flash_PBeta(Tc,Pc,w,kij,Ncomp):
#    eos = c_vEoS(Ncomp,Tc,Pc,w,kij)
#    return 

def flash_TBeta(T,P,z,BETA,K,Ncomp,eos):
    RES_flash=1
    TOL=1e-10
    j=0
    while (np.abs(RES_flash)>TOL):
        K_ol=1.*K
        x,y=update_x(z,K,BETA,Ncomp)
        P=Newton_P(z,K,BETA,P,x,y,Ncomp,T,eos)
        K = update_K(T,P,x,y,eos)
        RES_flash=np.linalg.norm(K_ol-K)
        j+=1
    return x,y,K,P,j

def Newton_P(z,K,BETA,P,x,y,Ncomp,T,eos):
    RES=1
    TOL=1e-8
    MAXi=100
    i=0
    while (np.abs(RES)>TOL and i < MAXi):
        RES=RES_RR_P(z,K,BETA,P,x,y,Ncomp,T,eos)
        step=1
        JAC=(RES_RR_P(z,K,BETA,P+step,x,y,Ncomp,T,eos)-RES_RR_P(z,K,BETA,P-step,x,y,Ncomp,T,eos))/(2*step)
        P-=RES/JAC
        i+=1
    return P

def update_x(z,K,BETA,Ncomp):
    x=np.zeros(Ncomp)
    y=np.zeros(Ncomp)
    for i in range(Ncomp):
        x[i] = z[i]*( (1.) / (1+BETA*(K[i]-1.)) )            
        y[i] = K[i]*x[i]
    return x/np.sum(x), y/np.sum(y)

def update_K(T,P,x,y,eos):
    VL=eos.Volume(T=T,P=P,x=x)[0]
    VV=eos.Volume(T=T,P=P,x=y)[1] 

    phiL=eos.fugacity_coeff(T=T,V=VL,x=x)
    phiV=eos.fugacity_coeff(T=T,V=VV,x=y)

    K=phiL/phiV
    return K

def RES_RR_P(z,K,BETA,P,x,y,Ncomp,T,eos):   
    RES = 0.
    K = update_K(T,P,x,y,eos)
    for i in range(Ncomp):
        RES += z[i]*( (K[i]-1.) / (1.+BETA*(K[i]-1.)) )            
    return RES