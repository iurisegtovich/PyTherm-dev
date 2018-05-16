import numpy as np
from scipy import optimize as opt

def calc_P_sat(vEoS_obj,T,iguess_P):#,index):
    
    #variáveis de método numérico
    RES=1
    TOL=1e-9
    MAX=1000

    #estimativa inicial
    P=iguess_P
    
    #contagem de iterações
    i=0
    while(RES>TOL and i<MAX): #Kernel > Interrupt (console > Ctrl+C)

        x=np.array([1.]) #zeros(ncomp), x[index]=1.
        fluidV=vEoS_obj.Volume(T=T, P=P, x=x)

        if (len(fluidV) == 1): #caso so tenha sido encontrada uma raiz, não será possível calcular P_sat para essa temperatura por esse método
            return np.nan #return sem variável pois não há solução, o código q chama essa função precisa saber lidar com esse resultado excepcional!
        [V_L,V_V] = fluidV
        phiL=vEoS_obj.fugacity_coeff(T=T, V=fluidV[0], x=x)
        phiV=vEoS_obj.fugacity_coeff(T=T, V=fluidV[1], x=x)
        P=P*(phiL/phiV)
        RES=np.abs(phiL/phiV-1.)
        i=i+1
    return P
    
def calc_Psat_curve(vEoS_obj,gridT,guess_P=100.):
    #grid feito usando a função linspace, gera pontos igualmente espaçados em escala linear. No caso: 100 pontos, entre 100 e Tc, inclusive.
    n=len(gridT)
    gridP=np.zeros(n)

    #efetua primeiro calculo
    gridP[0]=calc_P_sat(vEoS_obj,gridT[0],guess_P) #primeiro ponto

    #efetua demais cálculos
    for i in range(1,n): #demais pontos
        gridP[i]=calc_P_sat(vEoS_obj,gridT[i],gridP[i-1])
        #print(grid_P[i]) #habilitar essa linha - removendo o símbolo de comentário -- # -- - faz com que os resultados de cada iteração sejam exibidos na seção de impressão da célula
    return gridT, gridP
