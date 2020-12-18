###############################################################################
###############################################################################
#
# CODE TO SOLVE THE 2D CAHN-HILLIARD EQUATION
#
# AUTHOR OF THE CODE: SERGIO P. PEREZ
#
# FINITE-VOLUME SEMI-IMPLICIT DIMENSIONAL-SPLITTING SCHEME

###############################################################################
###############################################################################

import numpy as np

##################
# DEFINE FUNCTION: H CONTRACTIVE DERIVATIVE
##################
def Hc1_con(rho,pot,theta): 
    if pot==1:# Double well
        Hc1=rho**3
    elif pot==2 and theta!=0:# Logarithmic
        Hc1=theta/2.*np.log(np.divide(1+rho,1-rho))
#        greater1_index = np.append(rho,0) > 1
#        lowerminus1_index = np.append(rho,0) < -1
#        Hc1[greater1_index[:-1]]=9999999999999999999999999
#        Hc1[lowerminus1_index[:-1]]=-9999999999999999999999     
    else:
        Hc1=np.zeros(np.shape(rho))    
    return Hc1


##################
# DEFINE FUNCTION: H EXPANSIVE DERIVATIVE
##################    
def He1_exp(rho,pot,thetac):
    if pot==1:# Double well
        He1=rho
    elif pot==2 and thetac!=0:# Logarithmic
        He1=thetac*rho
#        greater1_index = np.append(rho,0) > 1
#        lowerminus1_index = np.append(rho,0) < -1
#        He1[greater1_index[:-1]]=-9999
#        He1[lowerminus1_index[:-1]]=+9999
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
# DEFINE FUNCTION: Wetting boundary condition
##################
def fw1_con(rho,bc,angle,epsilon): 
    if bc==2: 
        if np.cos(angle)>0:
            return epsilon*np.sqrt(2)/2*np.cos(angle)*(rho**2-1+2*rho)
        else:
            return -epsilon*np.sqrt(2)/2*np.cos(angle)*(2*rho)
    else:
        try:
            return np.zeros(len(rho))
        except:
            return 0
        
def fw1_exp(rho,bc,angle,epsilon): 
    if bc==2: 
        if np.cos(angle)>0:
            return epsilon*np.sqrt(2)/2*np.cos(angle)*(2*rho)       
        else:
            return -epsilon*np.sqrt(2)/2*np.cos(angle)*(rho**2-1+2*rho)
    else:
        try:
            return np.zeros(len(rho))
        except:
            return 0
  



##################
# DEFINE FUNCTION: FLUX 
##################


def Euler_implicit_row(rhor,rhor1,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob,bc,angle,j):
                       

    
    # Define variation of free energy
    
    # a) Hc: contractive (convex) part of the free energy (treated implicitly))
    # Hc1 is the first derivative of Hc
    
    Hc1= Hc1_con(rhor1,pot,theta)
    
    # b) He: expansive (concave) part of the free energy (treated explicitly)
    # He1 is the first derivative of He    
    
    He1= He1_exp(rhor[j,:],pot,thetac)
    
     # c) Laplacian (treated semi-implicitly)
     
    Lap=np.zeros(np.shape(rhor1)[0])
    

     
    if j==0:
        
        Lap[1:-1]=epsilon**2/dx**2/2.*(rhor[j,2:]-2*rhor[j,1:-1]+rhor[j,0:-2]\
           +rhor[j+1,1:-1]-rhor[j,1:-1]\
           +rhor1[2:]-2*rhor1[1:-1]+rhor1[0:-2]\
           +rhor[j+1,1:-1]-rhor1[1:-1])
        
      
        Lap[0]=epsilon**2/dx**2/2.*(rhor[j,1]-rhor[j,0]\
           +rhor[j+1,0]-rhor[j,0]\
           +rhor1[1]-rhor1[0]\
           +rhor[j+1,0]-rhor1[0])
        
        Lap[-1]=epsilon**2/dx**2/2.*(-rhor[j,-1]+rhor[j,-2]\
           +rhor[j+1,-1]-rhor[j,-1]\
           -rhor1[-1]+rhor1[-2]\
           +rhor[j+1,-1]-rhor1[-1])
        
        
         
    elif j==np.shape(rhor)[0]-1:
        
        Lap[1:-1]=epsilon**2/dx**2/2.*(rhor[j,2:]-2*rhor[j,1:-1]+rhor[j,0:-2]\
           -rhor[j,1:-1]+rhor[j-1,1:-1]\
           +rhor1[2:]-2*rhor1[1:-1]+rhor1[0:-2]\
           -rhor1[1:-1]+rhor[j-1,1:-1])
        
        
        Lap[0]=epsilon**2/dx**2/2.*(rhor[j,1]-rhor[j,0]\
           -rhor[j,0]+rhor[j-1,0]\
           +rhor1[1]-rhor1[0]\
           -rhor1[0]+rhor[j-1,0])
        
        Lap[-1]=epsilon**2/dx**2/2.*(-rhor[j,-1]+rhor[j,-2]\
           -rhor[j,-1]+rhor[j-1,-1]\
           -rhor1[-1]+rhor1[-2]\
           -rhor1[-1]+rhor[j-1,-1])
            
         
    else:
         
        Lap[1:-1]=epsilon**2/dx**2/2.*(rhor[j,2:]-2*rhor[j,1:-1]+rhor[j,0:-2]\
           +rhor[j+1,1:-1]-2*rhor[j,1:-1]+rhor[j-1,1:-1]\
           +rhor1[2:]-2*rhor1[1:-1]+rhor1[0:-2]\
           +rhor[j+1,1:-1]-2*rhor1[1:-1]+rhor[j-1,1:-1])
       
        Lap[0]=epsilon**2/dx**2/2.*(rhor[j,1]-rhor[j,0]\
           +rhor[j+1,0]-2*rhor[j,0]+rhor[j-1,0]\
           +rhor1[1]-rhor1[0]\
           +rhor[j+1,0]-2*rhor1[0]+rhor[j-1,0])
        
        Lap[-1]=epsilon**2/dx**2/2.*(-rhor[j,-1]+rhor[j,-2]\
           +rhor[j+1,-1]-2*rhor[j,-1]+rhor[j-1,-1]\
           -rhor1[-1]+rhor1[-2]\
           +rhor[j+1,-1]-2*rhor1[-1]+rhor[j-1,-1])
           
    # Vector with wetting conditions
    Wx=np.zeros((len(rhor1))) 
    Wy=np.zeros((len(rhor1))) 
    
    Wx[0]=fw1_con(rhor1[0],bc,angle,epsilon)-fw1_exp(rhor[j,0],bc,angle,epsilon)
    Wx[-1]=fw1_con(rhor1[-1],bc,angle,epsilon)-fw1_exp(rhor[j,-1],bc,angle,epsilon)
    
    if j==0:       
        Wy=fw1_con(rhor1[:],bc,angle,epsilon)-fw1_exp(rhor[j,:],bc,angle,epsilon)
    elif j==np.shape(rhor)[0]-1:
        Wy=fw1_con(rhor1[:],bc,angle,epsilon)-fw1_exp(rhor[j,:],bc,angle,epsilon)
        
    # Compute (n-1) u velocities
    
    uhalf=-(Hc1[1:]-He1[1:]-Lap[1:]+Wx[1:]/(dx)+Wy[1:]/(dx)-Hc1[0:-1]+He1[0:-1]+Lap[0:-1]-Wx[0:-1]/(dx)-Wy[0:-1]/(dx))/dx
    
    # Upwind u velocities
    
    uhalfplus=np.zeros((len(rhor1)-1))
    uhalfminus=np.zeros((len(rhor1)-1))
    uhalfplus[uhalf > 0]=uhalf[uhalf > 0]
    uhalfminus[uhalf < 0]=uhalf[uhalf < 0]
    
    # Compute (n+1) fluxes, including no-flux boundary conditions
    
    Fhalf=np.zeros(len(rhor1)+1)
    
    # 1st order
    
    Fhalf[1:-1]=uhalfplus*mobility(rhor1[:-1],rhor1[1:],choicemob)+uhalfminus*mobility(rhor1[1:],rhor1[:-1],choicemob) 

    
    # Compute function equal to zero

    E_i=rhor1-rhor[j,:]+(Fhalf[1:]-Fhalf[:-1])*dt/dx#+penalty
   
    return E_i
