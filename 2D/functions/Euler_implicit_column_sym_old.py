# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:10:12 2019

@author: sp3215
"""
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
def fw(rho,bc,angle): 
    if bc==2: 
        return -np.sqrt(2)/2*np.cos(angle)*(1-rho**2)
    else:
        try:
            return np.zeros(len(rho))
        except:
            return 0
  



##################
# DEFINE FUNCTION: FLUX 
##################


def Euler_implicit_column_sym(rhoc,rhoc1,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob,bc,angle,j):
                       

    rhoc1=np.concatenate((rhoc1,np.flip(rhoc1)))
    
    # Define variation of free energy
    
    # a) Hc: contractive (convex) part of the free energy (treated implicitly))
    # Hc1 is the first derivative of Hc
    
    Hc1= Hc1_con(rhoc1,pot,theta)
    
    # b) He: expansive (concave) part of the free energy (treated explicitly)
    # He1 is the first derivative of He    
    
    He1= He1_exp(rhoc[:,j],pot,thetac)
    
     # c) Laplacian (treated semi-implicitly)
     
    Lap=np.zeros(np.shape(rhoc1)[0])
     
    if j==0:
        
        # Wetting condition
        Lap[1:-1]=epsilon**2/dx**2/2.*(rhoc[2:,j]-2*rhoc[1:-1,j]+rhoc[0:-2,j]\
           +rhoc[1:-1,j+1]-rhoc[1:-1,j]-dx*fw(rhoc[1:-1,j],bc,angle)/(epsilon)\
           +rhoc1[2:]-2*rhoc1[1:-1]+rhoc1[0:-2]\
           +rhoc[1:-1,j+1]-rhoc1[1:-1]-dx*fw(rhoc1[1:-1],bc,angle)/(epsilon))
        
        # No need of wetting condition since it's only one cell and in the corner
        Lap[0]=epsilon**2/dx**2/2.*(rhoc[1,j]-rhoc[0,j]\
           +rhoc[0,j+1]-rhoc[0,j]-dx*fw(rhoc[0,j],bc,angle)/(epsilon)\
           +rhoc1[1]-rhoc1[0]\
           +rhoc[0,j+1]-rhoc1[0]-dx*fw(rhoc1[0],bc,angle)/(epsilon))
        
        # No need of wetting condition since it's only one cell and in the corner
        Lap[-1]=epsilon**2/dx**2/2.*(-rhoc[-1,j]+rhoc[-2,j]\
           +rhoc[-1,j+1]-rhoc[-1,j]-dx*fw(rhoc[-1,j],bc,angle)/(epsilon)\
           -rhoc1[-1]+rhoc1[-2]\
           +rhoc[-1,j+1]-rhoc1[-1]-dx*fw(rhoc1[-1],bc,angle)/(epsilon))
        
        
         
    elif j==np.shape(rhoc)[1]-1:
        
        Lap[1:-1]=epsilon**2/dx**2/2.*(rhoc[2:,j]-2*rhoc[1:-1,j]+rhoc[0:-2,j]\
           -rhoc[1:-1,j]+rhoc[1:-1,j-1]\
           +rhoc1[2:]-2*rhoc1[1:-1]+rhoc1[0:-2]\
           -rhoc1[1:-1]+rhoc[1:-1,j-1])
        
        Lap[0]=epsilon**2/dx**2/2.*(rhoc[1,j]-rhoc[0,j]\
           -rhoc[0,j]+rhoc[0,j-1]\
           +rhoc1[1]-rhoc1[0]\
           -rhoc1[0]+rhoc[0,j-1])
        
        Lap[-1]=epsilon**2/dx**2/2.*(-rhoc[-1,j]+rhoc[-2,j]\
           -rhoc[-1,j]+rhoc[-1,j-1]\
           -rhoc1[-1]+rhoc1[-2]\
           -rhoc1[-1]+rhoc[-1,j-1])
            
         
    else:
         
        Lap[1:-1]=epsilon**2/dx**2/2.*(rhoc[2:,j]-2*rhoc[1:-1,j]+rhoc[0:-2,j]\
           +rhoc[1:-1,j+1]-2*rhoc[1:-1,j]+rhoc[1:-1,j-1]\
           +rhoc1[2:]-2*rhoc1[1:-1]+rhoc1[0:-2]\
           +rhoc[1:-1,j+1]-2*rhoc1[1:-1]+rhoc[1:-1,j-1])
        
        Lap[0]=epsilon**2/dx**2/2.*(rhoc[1,j]-rhoc[0,j]\
           +rhoc[0,j+1]-2*rhoc[0,j]+rhoc[0,j-1]\
           +rhoc1[1]-rhoc1[0]\
           +rhoc[0,j+1]-2*rhoc1[0]+rhoc[0,j-1])
        
        Lap[-1]=epsilon**2/dx**2/2.*(-rhoc[-1,j]+rhoc[-2,j]\
           +rhoc[-1,j+1]-2*rhoc[-1,j]+rhoc[-1,j-1]\
           -rhoc1[-1]+rhoc1[-2]\
           +rhoc[-1,j+1]-2*rhoc1[-1]+rhoc[-1,j-1])
           
        
       
    # Compute (n-1) u velocities
    
    uhalf=-(Hc1[1:]-He1[1:]-Lap[1:]-Hc1[0:-1]+He1[0:-1]+Lap[0:-1])/dx
    
    # Upwind u velocities
    
    uhalfplus=np.zeros((len(rhoc1)-1))
    uhalfminus=np.zeros((len(rhoc1)-1))
    uhalfplus[uhalf > 0]=uhalf[uhalf > 0]
    uhalfminus[uhalf < 0]=uhalf[uhalf < 0]
    
    # Compute (n+1) fluxes, including no-flux boundary conditions
    
    Fhalf=np.zeros(len(rhoc1)+1)
    
    # 1st order
    
    Fhalf[1:-1]=uhalfplus*mobility(rhoc1[:-1],rhoc1[1:],choicemob)+uhalfminus*mobility(rhoc1[1:],rhoc1[:-1],choicemob) 

    
    # Compute function equal to zero

    E_i=rhoc1-rhoc[:,j]+(Fhalf[1:]-Fhalf[:-1])*dt/dx#+penalty
   
    return E_i[:(np.shape(rhoc)[0])/2]
