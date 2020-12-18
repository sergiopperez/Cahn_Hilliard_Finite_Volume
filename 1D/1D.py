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
# CALL PACKAGES
##################

import numpy as np
from scipy import optimize
import time
import os
dirpath = os.getcwd()
import sys
sys.path.append(dirpath+'\\functions')
from plots4times import *
from plot_imshow import *
from initial_conditions import *
from Euler_implicit import *
from createplots import *
plt.rc('text',usetex=True)
plt.rc('font',family='serif')



##################
# DEFINE FUNCTION: INEQUALITY OF dFdt
##################

def inequality_dFdt(rho):
    
    H= rho*(rho**2-1)
    
    Lap=np.zeros(np.shape(rho))
    Lap[1:-1,:]=epsilon**2*(rho[0:-2,:]-2*rho[1:-1,:]+rho[2:,:])/dx**2
    Lap[0,:]=Lap[1,:]
    Lap[-1,:]=Lap[-2,:]
    
    uhalf=-(H[1:,:]-Lap[1:,:]-H[0:-1,:]+Lap[0:-1,:])/dx
    
    ineq=-dx*(uhalf**2).sum(axis=0)    
    
    return ineq



    
##################
# DEFINE FUNCTION: MAIN
##################
    
def main():
    
    global n, x, dx, dt,epsilon,ntimes,pot,theta,thetac,choicemob
    
    ##################
    # Select the configuration of the problem (see initial_conditions.py)
    
    choiceinitialcond=32 
    ##################
    
    n, dx, x, rho0, epsilon, dt, tmax, ntimes, pot, theta, thetac, choicemob = initial_conditions(choiceinitialcond)    
    
    rho=np.zeros([n,ntimes+1]) # Density matrix
    rho[:,0]=rho0 # First column of density matrix is initial density
    t=np.zeros(ntimes+1) # Time vector
    tic = time.clock()

    for i in np.arange(ntimes): #Temporal loop
        
        rho[:,i+1],infodict, ier, mesg=optimize.fsolve(lambda rhon: Euler_implicit(rho[:,i],rhon,n,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob), rho[:,i],full_output = True)#,xtol=1e-4)#, col_deriv=True, maxfev=10000*80,factor=1, full_output = False)
        
        t[i+1]=t[i]+dt
        
        print('--------------------')
        print('Time: ',t[i])
        print(['L1 norm of the difference between the new and old state: ',np.linalg.norm(rho[:,i+1]-rho[:,i],1)])
    
    
    # Save final density to compute error

#    rho[:,i].tofile('error/rho%s.dat'%n)
    
 
    # Compute free energy in time
    if pot==1:
        F=dx*((0.25*(rho[:,:]**2-1)**2).sum(axis=0)+(epsilon**2/2*(rho[2:,:]-rho[1:-1,:])**2/dx**2).sum(axis=0))

    elif pot==2:
        
        if theta!=0 and thetac==0:
            F=dx*((theta/2.*(np.multiply(1+rho[:,:],np.log((1+rho[:,:])/2.))+np.multiply(1-rho[:,:],np.log((1-rho[:,:])/2.)))\
                   ).sum(axis=0)+(epsilon**2/2*(rho[2:,:]-rho[1:-1,:])**2/dx**2).sum(axis=0)\
                   +(epsilon**2/2*(-rho[-1,:]+rho[-2,:])**2/dx**2)+(epsilon**2/2*(rho[1,:]-rho[0,:])**2/dx**2))
    
        elif theta<0.01 and thetac!=0:
            F=dx*((thetac/2.*(1-rho[:,:]**2)).sum(axis=0)+(epsilon**2/2*(rho[2:,:]-rho[1:-1,:])**2/dx**2).sum(axis=0)\
                   +(epsilon**2/2*(-rho[-1,:]+rho[-2,:])**2/dx**2)+(epsilon**2/2*(rho[1,:]-rho[0,:])**2/dx**2))
            
        elif theta!=0 and thetac!=0:
            F=dx*((theta/2.*(np.multiply(1+rho[:,:],np.log((1+rho[:,:])/2.))+np.multiply(1-rho[:,:],np.log((1-rho[:,:])/2.)))\
                   +thetac/2.*(1-rho[:,:]**2)).sum(axis=0)+(epsilon**2/2*(rho[2:,:]-rho[1:-1,:])**2/dx**2).sum(axis=0)\
                   +(epsilon**2/2*(-rho[-1,:]+rho[-2,:])**2/dx**2)+(epsilon**2/2*(rho[1,:]-rho[0,:])**2/dx**2))
                
    
    # Compute variation of the free energy
    
    dFdrho=np.zeros([n,ntimes+1])
    for i in np.arange(ntimes+1):
        Lap=np.zeros(n)
        Lap[1:-1]=epsilon**2*(rho[0:-2,i]-2*rho[1:-1,i]+rho[2:,i])/dx**2
        Lap[0]=epsilon**2*(-rho[0,i]+rho[1,i])/dx**2
        Lap[-1]=epsilon**2*(-rho[-2,i]+rho[-1,i])/dx**2
        
        
        dFdrho[:,i]=Hc1_con(rho[:,i],pot,theta)-He1_exp(rho[:,i],pot,thetac)-epsilon**2*Lap
                       
    
    # Compute time derivative free energy 
    
    #dFdt=(F[1:]-F[:-1])/(t[1:]-t[:-1])
    dFdt=(F[2:]-F[:-2])/(t[2:]-t[:-2])
        
    toc = time.clock()
    display(toc-tic)
    
    rho.tofile('data/rho_%s.dat'%choiceinitialcond)
    F.tofile('data/F_%s.dat'%choiceinitialcond)
    t.tofile('data/t_%s.dat'%choiceinitialcond)

    # Create plots
    
    createplots(x,rho,F,dFdrho,t,choiceinitialcond,theta,choicemob)


if __name__ == '__main__':
    main()