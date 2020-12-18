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
from initial_conditions import *
from Euler_implicit_row import *
from Euler_implicit_row_sym import *
from Euler_implicit_column import *
from Euler_implicit_column_sym import *
from createplots import *
plt.rc('text',usetex=True)
plt.rc('font',family='serif')


##################
# DEFINE FUNCTION: MAIN
##################
    
def main():
    
    global nx, ny, x, dx, dt,epsilon,ntimes,pot,theta,thetac,choicemob
    
    
    ##################
    # Select the configuration of the problem (see initial_conditions.py)
    
    choiceinitialcond=33 
    ##################
    
    nx, ny, dx, x, y, rho0, epsilon, dt, tmax, ntimes, pot, theta, thetac, choicemob, bc, angle = initial_conditions(choiceinitialcond)    
    
    rho=np.zeros([nx*ny,ntimes+1]) # Density matrix
    rho[:,0]=rho0 # First column of density matrix is initial density
    t=np.zeros(ntimes+1) # Time vector

    for i in np.arange(ntimes): #Temporal loop
        
        rhor=np.reshape(rho[:,i],(nx,ny)).copy()
  
        # Row by row
        for j in np.arange(nx):
            
            rhor1,infodict, ier, mesg=optimize.fsolve(lambda rhor1: Euler_implicit_row(rhor,rhor1,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob,bc,angle,j), rhor[j,:], full_output = True)
            
            rhor[j,:]=rhor1
            
#            rhor1,infodict, ier, mesg=optimize.fsolve(lambda rhor1: Euler_implicit_row_sym(rhor,rhor1,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob,bc,angle,j), rhor[j,:ny/2], full_output = True)
#            
#            rhor[j,:]=np.concatenate((rhor1, np.flip(rhor1)))
        
        # Column by column     
        rhoc=rhor
        for j in np.arange(ny):
            
#            rhoc1,infodict, ier, mesg=optimize.fsolve(lambda rhoc1: Euler_implicit_column(rhoc,rhoc1,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob,bc,angle,j), rhoc[:,j], full_output = True)
#            
#            rhoc[:,j]=rhoc1
            
            rhoc1,infodict, ier, mesg=optimize.fsolve(lambda rhoc1: Euler_implicit_column_sym(rhoc,rhoc1,dx,dt,epsilon,ntimes,pot,theta,thetac,choicemob,bc,angle,j), rhoc[:nx/2,j], full_output = True)
            
            rhoc[:,j]=np.concatenate((rhoc1, np.flip(rhoc1)))
        
        rho[:,i+1]=np.reshape(rhoc,nx*ny).copy()
        t[i+1]=t[i]+dt
        
        if i>0:
            print('--------------------')
            print('Time: ',t[i])
            print(['L1 norm of the difference between the new and old state: ',np.linalg.norm(rho[:,i]-rho[:,i-1],1)])
    
  
    # Compute free energy in time
    if pot==1:
        F=dx*((0.25*(rho[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
           
    elif pot==2:
        
        if theta!=0 and thetac==0:
            F=dx**2*((theta/2.*(np.multiply(1+rho[:,:],np.log((1+rho[:,:])/2.))+np.multiply(1-rho[:,:],np.log((1-rho[:,:])/2.)))\
                   ).sum(axis=0)+(epsilon**2/2*(rho[2:,:]-rho[1:-1,:])**2/dx**2).sum(axis=0)\
                   +(epsilon**2/2*(-rho[-1,:]+rho[-2,:])**2/dx**2)+(epsilon**2/2*(rho[1,:]-rho[0,:])**2/dx**2))
            
        elif theta==0 and thetac!=0:
            F=dx**2*((thetac/2.*(1-rho[:,:]**2)).sum(axis=0)\
                  +(epsilon**2/2*np.reshape((np.reshape(rho,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
                  +(epsilon**2/2*np.reshape((np.reshape(rho,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        elif theta!=0 and thetac!=0:
            F=dx**2*((theta/2.*(np.multiply(1+rho[:,:],np.log((1+rho[:,:])/2.))+np.multiply(1-rho[:,:],np.log((1-rho[:,:])/2.)))\
                   +thetac/2.*(1-rho[:,:]**2)).sum(axis=0)+(epsilon**2/2*(rho[2:,:]-rho[1:-1,:])**2/dx**2).sum(axis=0)\
                   +(epsilon**2/2*(-rho[-1,:]+rho[-2,:])**2/dx**2)+(epsilon**2/2*(rho[1,:]-rho[0,:])**2/dx**2))
    
    # Add wetting free energy
    rho=np.reshape(rho,(nx,ny,ntimes+1))
    F=F+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho[:,0,:]**3/3-rho[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho[:,-1,:]**3/3-rho[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho[0,:,:]**3/3-rho[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho[-1,:,:]**3/3-rho[-1,:,:])).sum(axis=0))
                
    rho=np.reshape(rho,(nx*ny,ntimes+1))

    
    # Save final density, free energy
    if choiceinitialcond==11 or choiceinitialcond==12:
        rho.tofile('data/rho_%s_%s.dat'%(nx,choiceinitialcond))
    else:
        rho.tofile('data/rho_%s.dat'%choiceinitialcond)
    F.tofile('data/F_%s.dat'%choiceinitialcond)
    t.tofile('data/t_%s.dat'%choiceinitialcond)

    
    # Plot initial density profile
    
    plt.close()       
    figure1 = plt.figure(figsize=(10, 6.5))
    plt.clf()
    plt.title(r'Initial $\phi$ ', fontsize=30)
    plt.imshow(np.flipud(np.reshape(rho[:,0],(nx,ny)).T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
    plt.colorbar()
    plt.ylabel(r'$y$', fontsize=23)
    plt.xlabel(r'$x$', fontsize=23)
    #plt.show()
    plt.tight_layout()
    figure1.savefig('phi_initial.pdf', bbox_inches='tight')
    
    # Plot final density profile
      
    figure2 = plt.figure(figsize=(10, 6.5))
    plt.clf()
    plt.title(r'Final $\phi$ ', fontsize=30)
    plt.imshow(np.flipud(np.reshape(rho[:,-1],(nx,ny)).T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
    plt.colorbar()
    plt.ylabel(r'$y$', fontsize=23)
    plt.xlabel(r'$x$', fontsize=23)
    #plt.show()
    plt.tight_layout()
    figure2.savefig('phi_final.pdf', bbox_inches='tight')
    
    # Plot temporal evolution of the free energy

    figure3 = plt.figure(figsize=(10, 6.5))
    plt.clf()
    plt.title(r'Decay of free energy', fontsize=30)
    plt.plot(t,F);
    plt.ylabel(r'$\mathcal{F}[\phi]$', fontsize=23)
    plt.xlabel(r'$t$', fontsize=23)
    #plt.show()
    plt.tight_layout()
    figure3.savefig('freeenergy.pdf', bbox_inches='tight')
    



if __name__ == '__main__':
    main()