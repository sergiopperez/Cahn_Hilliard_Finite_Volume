###############################################################################
###############################################################################
#
# CODE TO SOLVE THE 1D CAHN-HILLIARD EQUATION
#
# AUTHOR OF THE CODE: SERGIO P. PEREZ
#
# FINITE-VOLUME SEMI-IMPLICIT SCHEME

###############################################################################
###############################################################################

##################
# DEFINE FUNCTION: Euler_implicit
##################
import numpy as np

##################
# DEFINE FUNCTION: H CONTRACTIVE DERIVATIVE
##################
def Hc1_con(rho,pot,theta): 
    if pot==1:# Double well
        Hc1=rho**3
        #Hc1=4*rho**3+2*rho
    elif pot==2 and theta!=0:# Logarithmic
        #Hc1=np.zeros(np.shape(rho)[0])
        Hc1=theta/2.*np.log(np.divide(1+rho,1-rho))
    else:
        Hc1=np.zeros(np.shape(rho))  
        
    # Uncomment the following for number==31 and number==32
    greater1_index = np.append(rho,0) > 1
    lowerminus1_index = np.append(rho,0) < -1
    Hc1[greater1_index[:-1]]=99999999999
    Hc1[lowerminus1_index[:-1]]=-99999999999 
    return Hc1

##################
# DEFINE FUNCTION: H CONTRACTIVE 2 DERIVATIVE 
##################
def Hc2_con(rho,pot,theta):
    if pot==1:# Double well
        Hc2=3*rho**2
        #Hc2=12*rho**2+2
    elif pot==2 and theta!=0:# Logarithmic
        if np.abs(rho)<=1:        
            Hc2=theta/(1-rho**2)
        else:
            Hc2=9999999999       
        #greaterabs1_index = np.append(np.abs(rho),0) > 1
        #Hc2[greaterabs1_index[:-1]]=99999999999
    else:
        Hc2=np.zeros(np.shape(rho))
    return Hc2

##################
# DEFINE FUNCTION: H EXPANSIVE 1 DERIVATIVE
##################    
def He1_exp(rho,pot,thetac):
    if pot==1:# Double well
        He1=rho
        #He1=6*rho
    elif pot==2 and thetac!=0:# Logarithmic
        He1=thetac*rho
#        greater1_index = np.append(rho,0) > 1
#        lowerminus1_index = np.append(rho,0) < -1
##        He1[greater1_index[:-1]]=-1.5
##        He1[lowerminus1_index[:-1]]=+1.5
#        He1[greater1_index[:-1]]=thetac*rho[greater1_index[:-1]]**3
#        He1[lowerminus1_index[:-1]]=thetac*rho[lowerminus1_index[:-1]]**3
        
    else:
        He1=np.zeros(np.shape(rho))
    return He1

##################
# DEFINE FUNCTION: MOBILITY
##################
def mobility(rho1, rho2,choicemob): 
    if choicemob==1:
        m=np.zeros(np.shape(rho1))
        m[:]=1
    elif choicemob==2:
        m=(1+rho1)*(1-rho2)
    return m

##################
# DEFINE FUNCTION: FLUX 
##################


def Euler_implicit(rho,rhon,n,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob):
    
    # Define variation of free energy
    
    # a) Hc: contractive (convex) part of the free energy (treated implicitly))
    # Hc1 is the first derivative of Hc
    
    Hc1= Hc1_con(rhon,pot,theta)
    
    # b) He: expansive (concave) part of the free energy (treated explicitly)
    # He1 is the first derivative of He    
    
    He1= He1_exp(rho,pot,thetac)
    
    # c) Laplacian (treated semi-implicitly)
    
    Lap=np.zeros(np.shape(rho)[0])
    Lap[1:-1]=epsilon**2*(rho[0:-2]-2*rho[1:-1]+rho[2:]+rhon[0:-2]-2*rhon[1:-1]+rhon[2:])/dx**2/2.
    Lap[0]=epsilon**2*(-rho[0]+rho[1]-rhon[0]+rhon[1])/dx**2/2.
    Lap[-1]=epsilon**2*(+rho[-2]-rho[-1]+rhon[-2]-rhon[-1])/dx**2/2.
    

    
    # Compute (n-1) velocities
    
    uhalf=-(Hc1[1:]-He1[1:]-Lap[1:]-Hc1[0:-1]+He1[0:-1]+Lap[0:-1])/dx
    
    # Upwind velocities
    
    positive_indices_u = uhalf > 0
    negative_indices_u = uhalf < 0
    
    uhalfplus=np.zeros(n-1)
    uhalfminus=np.zeros(n-1)
    uhalfplus[positive_indices_u]=uhalf[positive_indices_u]
    uhalfminus[negative_indices_u]=uhalf[negative_indices_u]
    
    

    # Compute (n+1) fluxes, including no-flux boundary conditions
    
    Fhalf=np.zeros(n+1)
    
    # 1st order
    Fhalf[1:-1]=uhalfplus*mobility(rhon[:-1],rhon[1:],choicemob)+uhalfminus*mobility(rhon[1:],rhon[:-1],choicemob) 
    

    # Compute function equal to zero

    E_i=rhon-rho+(Fhalf[1:]-Fhalf[:-1])*dt/dx#+penalty
   
    return E_i


