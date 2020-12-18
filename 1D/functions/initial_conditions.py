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

import numpy as np
from scipy.integrate import quad

##################
# DEFINE FUNCTION: initial_conditions
##################
    
def initial_conditions(choice):
    
    if choice==11: #Random initial configuration
        n=400 # Number of cells
        dx=0.4 # Width of cells
        x= np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
        #rho0=0.5*np.random.random_sample([n])-0.25 #Initial density
        #rho0.tofile('rho11.dat')
        rho0=np.fromfile('rho11.dat')
        epsilon=1.#1 # Parameter epsilon   
        dt=0.01 # Time step
        tmax=30.
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.3 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
        # Add modification to Hc1!         
        
    
    elif choice==12: #Random initial configuration
        n=200 # Number of cells
        dx=0.4 # Width of cells
        x= np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
        #rho=0.5*np.random.random_sample([n])-0.25 #Initial density
        rho0=np.fromfile('rho11.dat')
        epsilon=1. # Parameter epsilon   
        dt=0.01 # Time step
        tmax=30
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.5 # Absolut temperature for the logarithmic potential
        thetac=2. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2 
        
        # Add modification to Hc1! 
               
    elif choice==2: # Steady state from paper of Local Discontinuous Galerkin
        
#        l=1
#        dx=2.**(-6-l) # Width of cells
#        n=int(1./dx) # Number of cells
#        x= np.linspace(0,n*dx,n) # Spatial mesh

        
        n=100 # Number of cells
        x= np.linspace(0,1,n) # Spatial mesh  
        dx=x[2]-x[1] # Width of cells         
        
        
        epsilon=0.1 # Parameter epsilon   
        rho0=np.zeros(n)
        for k in np.arange(n):
            if np.abs(x[k]-1./2.)<=np.pi*epsilon/2.:
                #rho0[k]=np.cos((x[k]-0.5)/epsilon)-1.
                #rho0[k],error=quad(lambda y: (np.cos((y-0.5)/epsilon)-1)/dx,x[k]-dx/2.,x[k]+dx/2.,epsabs=1e-15)
                rho0[k]=(np.cos((x[k]-0.5)/epsilon)-1)
            else:
                rho0[k]=-1.0                             
#        dt=40.96*dx**2 # Time step
#        dt=0.0000003        
        dt=0.0001 # Time step for error analysis
        #dt=np.sqrt(dx)
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps      
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.0 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
        
    elif choice==31: # Example of 2 bumps from paper of Local Discontinuous Galerkin
                     # theta=0.3 and mob=1
         
        n=80 # Number of cells
        x= np.linspace(0,1,n) # Spatial mesh  
        dx=x[2]-x[1] # Width of cells         
        
        
        epsilon=0.03162277660168379 # Parameter epsilon   epsilon**2=1e-3
        rho0=np.zeros(n)
        
        for k in np.arange(n):
            if x[k]>=0 and x[k]<=1./3.-1./20.:
                rho0[k]=1.-0.00001
            elif np.abs(x[k]-1./3.)<=1./20.:
                rho0[k],error=quad(lambda y: (20.*(1./3.-y))/dx,x[k]-dx/2.,x[k]+dx/2.)
            elif np.abs(x[k]-41./50.)<=1./20.:
                rho0[k],error=quad(lambda y: (-20.*np.abs(y-41./50.))/dx,x[k]-dx/2.,x[k]+dx/2.)
            else:
                rho0[k]=-1.0+0.00001    
                     
#        dt=40.96*dx**2 # Time step
        dt=0.0001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps      
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.3 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=1 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
        # Add modification to Hc1! You can use first or second order
        
        
    elif choice==32: # Example of 2 bumps from paper of Local Discontinuous Galerkin
                     # theta=0.3 and mob=1-rho**2
         
        n=80 # Number of cells
        x= np.linspace(0,1,n) # Spatial mesh  
        dx=x[2]-x[1] # Width of cells         
        
        
        epsilon=0.03162277660168379 # Parameter epsilon   epsilon**2=1e-3
        rho0=np.zeros(n)
        
        for k in np.arange(n):
            if x[k]>=0 and x[k]<=1./3.-1./20.:
                rho0[k]=1.-0.00001
            elif np.abs(x[k]-1./3.)<=1./20.:
                rho0[k],error=quad(lambda y: (20.*(1./3.-y))/dx,x[k]-dx/2.,x[k]+dx/2.)
            elif np.abs(x[k]-41./50.)<=1./20.:
                rho0[k],error=quad(lambda y: (-20.*np.abs(y-41./50.))/dx,x[k]-dx/2.,x[k]+dx/2.)
            else:
                rho0[k]=-1.0+0.00001    
                     
#        dt=40.96*dx**2 # Time step
        dt=0.0001 # Time step
        tmax=2.0
        ntimes=int(tmax/dt)# 4000 # Number of time steps      
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.3 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
        # Add modification to Hc1 and use second order
        
    
    elif choice==33: # Example of 2 bumps from paper of Local Discontinuous Galerkin
                     # theta=0.0 and mob=1
                     
                     # This one has NOT worked
         
        n=40 # Number of cells
        x= np.linspace(0,1,n) # Spatial mesh  
        dx=x[2]-x[1] # Width of cells         
        
        
        epsilon=0.03162277660168379 # Parameter epsilon   epsilon**2=1e-3
        #epsilon=0.1
        rho0=np.zeros(n)
        
        for k in np.arange(n):
            if x[k]>=0 and x[k]<=1./3.-1./20.:
                rho0[k]=1.-0.01
            elif np.abs(x[k]-1./3.)<=1./20.:
                rho0[k],error=quad(lambda y: (20.*(1./3.-y))/dx,x[k]-dx/2.,x[k]+dx/2.)
            elif np.abs(x[k]-41./50.)<=1./20.:
                rho0[k],error=quad(lambda y: (-20.*np.abs(y-41./50.))/dx,x[k]-dx/2.,x[k]+dx/2.)
            else:
                rho0[k]=-1.0+0.01    
                     
#        dt=40.96*dx**2 # Time step
        dt=0.00001 # Time step
        tmax=0.02
        ntimes=int(tmax/dt)# 4000 # Number of time steps      
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.001 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=1 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
        
    elif choice==34: # Example of 2 bumps from paper of Local Discontinuous Galerkin
                     # theta=0.0 and mob=1-rho**2
         
        n=80 # Number of cells
        x= np.linspace(0,1,n) # Spatial mesh  
        dx=x[2]-x[1] # Width of cells         
        
        
        epsilon=0.03162277660168379
        # Parameter epsilon   epsilon**2=1e-3
        rho0=np.zeros(n)
        
        for k in np.arange(n):
            if x[k]>=0 and x[k]<=1./3.-1./20.:
                rho0[k]=1.-0.00001
            elif np.abs(x[k]-1./3.)<=1./20.:
                rho0[k],error=quad(lambda y: (20.*(1./3.-y))/dx,x[k]-dx/2.,x[k]+dx/2.)
            elif np.abs(x[k]-41./50.)<=1./20.:
                rho0[k],error=quad(lambda y: (-20.*np.abs(y-41./50.))/dx,x[k]-dx/2.,x[k]+dx/2.)
            else:
                rho0[k]=-1.0+0.00001    
                     
#        dt=40.96*dx**2 # Time step
        dt=0.0001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps      
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.0001 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
        # No modification to Hc1, and use first order
        
        
    elif choice==35: # Example of 2 bumps from paper of Local Discontinuous Galerkin
                     # Double-well and mob=1
         
        n=80 # Number of cells
        x= np.linspace(0,1,n) # Spatial mesh  
        dx=x[2]-x[1] # Width of cells         
        
        
        epsilon=0.03162277660168379
        # Parameter epsilon   epsilon**2=1e-3
        rho0=np.zeros(n)
        
        for k in np.arange(n):
            if x[k]>=0 and x[k]<=1./3.-1./20.:
                rho0[k]=1.-0.00001
            elif np.abs(x[k]-1./3.)<=1./20.:
                rho0[k],error=quad(lambda y: (20.*(1./3.-y))/dx,x[k]-dx/2.,x[k]+dx/2.)
            elif np.abs(x[k]-41./50.)<=1./20.:
                rho0[k],error=quad(lambda y: (-20.*np.abs(y-41./50.))/dx,x[k]-dx/2.,x[k]+dx/2.)
            else:
                rho0[k]=-1.0+0.00001    
                     
#        dt=40.96*dx**2 # Time step
        dt=0.0001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps      
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.0001 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=1 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
        # No modification to Hc1, and use first order
        
    elif choice==36: # Example of 2 bumps from paper of Local Discontinuous Galerkin
                     # Double-well and mob=1-rho**2
         
        n=80 # Number of cells
        x= np.linspace(0,1,n) # Spatial mesh  
        dx=x[2]-x[1] # Width of cells         
        
        
        epsilon=0.03162277660168379
        # Parameter epsilon   epsilon**2=1e-3
        rho0=np.zeros(n)
        
        for k in np.arange(n):
            if x[k]>=0 and x[k]<=1./3.-1./20.:
                rho0[k]=1.
            elif np.abs(x[k]-1./3.)<=1./20.:
                rho0[k],error=quad(lambda y: (20.*(1./3.-y))/dx,x[k]-dx/2.,x[k]+dx/2.)
            elif np.abs(x[k]-41./50.)<=1./20.:
                rho0[k],error=quad(lambda y: (-20.*np.abs(y-41./50.))/dx,x[k]-dx/2.,x[k]+dx/2.)
            else:
                rho0[k]=-1.0 
                     
#        dt=40.96*dx**2 # Time step
        dt=0.0001 # Time step
        tmax=0.2
        ntimes=int(tmax/dt)# 4000 # Number of time steps      
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.0001 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2
        
    
    return n, dx, x, rho0, epsilon, dt, tmax, ntimes, pot, theta, thetac, choicemob
        
        