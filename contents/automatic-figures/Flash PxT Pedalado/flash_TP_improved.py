#with improvements presented by Martin Cismondi, MissingREF!

import numpy
#
R = 83.14462175 #; //R = cm³.bar/(mol.K);
#
def fugacity(Ncomp, T, P, a, b, kij, frac, root, index):
    bmix = 0.
    amix = 0.
    Samix = 0.
    
    for i in range(Ncomp):
        print(i)
        bmix = bmix + frac[i]*b[i]

    for i in range(Ncomp):
        for j in range(Ncomp):
            amix = amix + frac[i]*frac[j]*numpy.sqrt(a[i]*a[j])*(1-kij[i*Ncomp+j]);
            Samix = Samix + 2.*frac[i]*numpy.sqrt(a[index]*a[i])*(1-kij[index*Ncomp+i]);
#
    coefCubic = [
    1,
    -(1.-(bmix*P/(R*T)))*(P*P/(R*T*R*T))/(P*P*P/(R*T*R*T*R*T)),
    ((amix*P/(R*T*R*T))-2.*(bmix*P/(R*T))-3.*(bmix*P/(R*T))*(bmix*P/(R*T)))*(P/(R*T))/(P*P*P/(R*T*R*T*R*T)),
    -((amix*P/(R*T*R*T))*(bmix*P/(R*T))-(bmix*P/(R*T))*(bmix*P/(R*T))-(bmix*P/(R*T))*(bmix*P/(R*T))*(bmix*P/(R*T)))/(P*P*P/(R*T*R*T*R*T))
    ]
    V=numpy.concatenate(([0],numpy.roots(coefCubic)))
    if (root == 0 or (root == 1 and V[3] == 0. and V[2] == 0.)):
        V[0] = V[1];
    elif (root == 1 and V[3] != 0. and V[2] != 0.):
        V[0] = V[3];

    Z = (P*V[0])/(R*T);

    FugCoef = numpy.exp((b[index]/bmix)*(Z-1)-numpy.log(P*(V[0]-bmix)/(R*T))-(amix/(2.*numpy.sqrt(2)*bmix*R*T))*(Samix/amix-b[index]/bmix)*numpy.log((V[0]+(1+numpy.sqrt(2))*bmix)/(V[0]+(1-numpy.sqrt(2))*bmix)))
 
    return FugCoef

#//********MAIN PROGRAM**************//
T=298
P=10

Ncomp=3
      
z = numpy.array([.1,.2,.7])
Tc = numpy.array([300,270,500])
Pc = numpy.array([10, 50, 80])
w = numpy.array([0,0.1,.05])
kij = numpy.array([0,0,0,0,0,0,0,0,0])


def lo_w(w):
        return 0.3764+1.5423*w-0.2699*w*w
    
def hi_w(w):
    return 0.3796+1.4850*w-0.1644*w*w+0.01667*w*w*w
    
#def fa(aux):
#    alfa=(1.+aux*(1.-numpy.sqrt(T/Tc)))*(1.+aux*(1.-numpy.sqrt(T/Tc)));
#    return 0.45724*R*R*Tc*Tc*alfa/Pc
#    
#def fb(aux):
#    return 0.07780*R*Tc/Pc

def fab(aux):
    alfa=(1.+aux*(1.-numpy.sqrt(T/Tc)))*(1.+aux*(1.-numpy.sqrt(T/Tc)));
    return [0.45724*R*R*Tc*Tc*alfa/Pc, 0.07780*R*Tc/Pc]
    
#[a,b]=numpy.where(w<=0.09,
#[fa(lo_w(w)),fb(lo_w(w))],
#[fa(hi_w(w)),fb(hi_w(w))]
#)
#print(a,b)

[a,b]=numpy.where(w<=0.09,
fab(lo_w(w)),
fab(hi_w(w))
)

x = numpy.zeros(Ncomp)

y = numpy.zeros(Ncomp)

K = numpy.exp(numpy.log(Pc/P)+5.373*(1.+w)*(1.-Tc/T))

while True:
    g = -1.;
    for i in range(Ncomp):
        g += K[i]*z[i];
    if(g < 0.):
        for i in range(Ncomp):
            K[i] = 1.1*K[i];
        g = -1.;
        for i in range(Ncomp):
            g += K[i]*z[i];
    if(g < 0.):
        break
   
while True:
    g = 1.
    for i in range(Ncomp):
        g -= z[i]/K[i]
    if(g > 1):
        for i in range(Ncomp):
            K[i] = 0.9*K[i]
        g = 1.;
        for i in range(Ncomp):
            g -= z[i]/K[i]
    if(g > 1.):
        break
#
#   //Determining Beta range
#   //Light Ncomponent
BetaIndex = 0;
for i in range(Ncomp):
    BetaIndex = (K[BetaIndex] > K[i] ? BetaIndex : i);

if(K[BetaIndex]*z[BetaIndex] > 1.)
    BetaMin = (K[BetaIndex]*z[BetaIndex]-1.)/(K[BetaIndex]-1.);
else
    BetaMin = 0.;

BetaMin = (BetaMin < 0. ? 0. : BetaMin);

#   //Heavy Ncomponent
BetaIndex = 0;

for i in range(Ncomp): BetaIndex = (K[BetaIndex] < K[i] ? BetaIndex : i);

if(K[BetaIndex] < z[BetaIndex]) BetaMax = (1.-z[BetaIndex])/(1.-K[BetaIndex]);
else BetaMax = 1.;

BetaMax = (BetaMax > 1. ? 1. : BetaMax);
#
beta = (BetaMin+BetaMax/2.);
#   
SBK = calloc(Ncomp,sizeof(double));

  
#   //Main loop
tol = 1.E-6;
maxit = 50;
itNRtot = itSStot = 0;

while True:
#   {
    itSStot++;
    it = 0;
    while True:
#       {
        itNRtot++;
        it++;

#         //Newton Method
        SXIR = SXIJ = drdb1 = drdb2 = F = F1 = F2 = dF = 0.;
        for i in range(Ncomp):
#         {
            SBK[i] = 1. + beta*(K[i]-1);
             xir = z[i]/SBK[i];
             SXIR += xir;
             SXIJ += xir*K[i];
#         }
        for i in range(Ncomp):
#         {
#             //F += z[i]*(K[i]-1.)/(1.-beta+beta*K[i]); //RR
             xir = z[i]/SBK[i];
             F1 += xir;
             F2 += xir*K[i];
             //dF += z[i]*(K[i]-1.)*(K[i]-1.)/(1.-beta+beta*K[i])/(1.-beta+beta*K[i]); //RR
             dxijdb = z[i]*(K[i]-1.)/(SBK[i]*SBK[i]);
             drdb1 += dxijdb;
             drdb2 += dxijdb*K[i];
#         }
        drdb1 = -drdb1;
        drdb2 = -drdb2;
        dF = (1./SXIJ)*drdb1 - (SXIR/(SXIJ*SXIJ))*drdb2; //mSdC
        F = (F1/F2)-1.; //mSdC
#         //dF = -dF; //RR
#         
        step = F/dF;
        beta -= step;
         
        if(beta < BetaMin)
            beta += 0.5*step;
         
        if(beta > BetaMax)
             beta += 0.5*step;
        if (fabs(step) > tol and it < maxit):
            break
#       }
            
    if(it == maxit)
#       {
        printf("\n\nMaximum number of iterations reached.\n");
        exit(1);
#       }
#
#      
    normx = normy = 0.;
    for i in range(Ncomp):
#      {
        x[i] = z[i]/(1.-beta+beta*K[i]);
        normx += x[i];
        y[i] = K[i]*z[i]/(1.-beta+beta*K[i]);
        normy += y[i];
#      }
#      
    for i in range(Ncomp):
#      {
        x[i] /= normx;
        y[i] /= normy;
#      }
#      
    maxdK = 0.;
    for i in range(Ncomp):
#      {
         fugacity(Ncomp,T,P,&FugCoefLiq,a,b,kij,x,0,i);
         fugacity(Ncomp,T,P,&FugCoefVap,a,b,kij,y,1,i);
         Knew = FugCoefLiq/FugCoefVap;
         maxdK = (fabs(Knew-K[i]) > maxdK ? fabs(Knew-K[i]) : maxdK);
         K[i] = Knew;
#//colocar aqui controle de variação de K[i]////////////////////////////////////////////////////////////////
#      } 
    if(fabs(maxdK) > tol):
        break
#   }
if(beta < 0.)
#   {
      printf("\nSubcooled Liquid");
      beta = 0.;
      for i in range(Ncomp):
#      {
         x[i] = z[i];
         y[i] = 0.;
#      }
#   }
#
if(beta > 1.)
#   {
      printf("\nSuperheated Vapor");
      beta = 1.;
      for i in range(Ncomp):
#      {
         y[i] = z[i];
         x[i] = 0.;
#      }
#   }
#
#
#   //Printing final results
printf("\nTotal number of iterations:\n"
          "Newton-Raphson = %d\n"
          "Successive substitution = %d\n",itNRtot,itSStot);
printf("\nBeta = %7.6lf",beta);
printf("\n\n  Vapor     Liquid   Global Ncomposition\n");

for i in range(Ncomp):
    printf("%7.6lf   %7.6lf      %7.6lf\n",y[i],x[i],z[i]);
printf("\n");
#

#

test_Ncomp=3
test_T=280
test_P=10
test_a=[10,20,30]
test_b=[1e-5,2e-5,3e-5]
test_kij=[0,0,0,0,0,0,0,0,0]
test_frac=[.1,.2,.7]
test_root=1
test_index=1

test_phi = fugacity(test_Ncomp,
test_T,
test_P,
test_a,
test_b,
test_kij,
test_frac,
test_root,
test_index)

print("test_phi=",test_phi)
