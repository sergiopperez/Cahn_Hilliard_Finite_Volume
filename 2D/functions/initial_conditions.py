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
from scipy.integrate import quad

##################
# DEFINE FUNCTION: INITIAL CONDITIONS
##################
    
def initial_conditions(choice):
    
    #############################################
    # 1) Order-of-convergence test with initial conditions as in the 1D case
    #############################################
    
    # Radial symmetry
    if choice==11: #Random initial configuration
        nx=10 # Number of cells per row
        ny=nx # Number of cells per column
        x= np.linspace(0,1,nx) # Spatial mesh 1D
        y= np.linspace(0,1,ny)
        dx=x[1]-x[0]
        rho0=-np.ones((nx,ny))
        epsilon=0.1 # Parameter epsilon 
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if np.abs(np.sqrt((x[i]-0.5)**2+(y[j]-0.5)**2))<=np.pi*epsilon/2.:
                    rho0[i,j]=np.cos(np.sqrt((x[i]-0.5)**2+(y[j]-0.5)**2)/epsilon)-1
       
        import matplotlib.pyplot as plt  
        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
  
#        dt=0.00002 # Time step for nx=160
        dt=0.0001 # Time step for nx<160
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.0 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=1 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=np.pi/2 # Angle for wetting boundary condition
        
    # x-symmetry
    elif choice==12: #Random initial configuration
        nx=160 # Number of cells per row
        ny=nx # Number of cells per column
        x= np.linspace(0,1,nx) # Spatial mesh 1D
        y= np.linspace(0,1,ny)
        dx=x[1]-x[0]
        rho0=-np.ones((nx,ny))
        epsilon=0.1 # Parameter epsilon 
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if np.abs((x[i]-0.5))<=np.pi*epsilon/2.:
                    rho0[i,j]=np.cos((x[i]-0.5)/epsilon)-1
       
        import matplotlib.pyplot as plt  
        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
  
#        dt=0.000005 # Time step for nx=160
        dt=0.0001 # Time step for nx<160
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.0 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=1 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=np.pi/2 # Angle for wetting boundary condition
        
        
        
    #############################################
    # 2) Random field with parameters from Barrett 1999
    #############################################
        
    # Constant mobility
    elif choice==21: #Random initial configuration
        n=256 # Number of cells per row
        nx=n
        ny=n
        dx=1/256. # Width of cells
        x= np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
        y=x
        rho0=np.fromfile('data/rho0_22.dat')
        epsilon=np.sqrt(0.00032) # Parameter epsilon   
        dt=0.0016 # Time step
        tmax=1.0
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=1 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=1 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=0 # Angle for wetting boundary condition
        
    # Degenerate mobility
    elif choice==22: #Random initial configuration
        n=256 # Number of cells per row
        nx=n
        ny=n
        dx=1/256. # Width of cells
        x= np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
        y=x
        rho0=np.fromfile('data/rho0_22.dat')
#        rho0=(0.05*np.random.random_sample([n*n])-0.025)-0.4
#        rho0.tofile('data/rho0_22.dat')
        epsilon=np.sqrt(0.00032) # Parameter epsilon   
        dt=0.0016 # Time step
        tmax=1.0
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=1 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=0 # Angle for wetting boundary condition
    
    # Constant mobility and deep-quench (Didn't work!)
    elif choice==23: #Random initial configuration
        n=256 # Number of cells per row
        nx=n
        ny=n
        dx=1/256. # Width of cells
        x= np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
        y=x
        rho0=np.fromfile('data/rho0_22.dat')
        epsilon=np.sqrt(0.00032) # Parameter epsilon   
        dt=0.0016 # Time step
        tmax=0.64
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.001 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=1 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=1 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=0 # Angle for wetting boundary condition
    
    # Degenerate mobility and deep-quench (Didn't work!)
    elif choice==24: #Random initial configuration
        n=256 # Number of cells per row
        nx=n
        ny=n
        dx=1/256. # Width of cells
        x= np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
        rho0=np.fromfile('data/rho0_22.dat')
        epsilon=np.sqrt(0.00032) # Parameter epsilon   
        dt=0.0016 # Time step
        tmax=0.64
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=2 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=1 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=0 # Angle for wetting boundary condition
        
    #############################################
    # 3) Wetting phenomena from Aymard 2019
    #############################################
        
    # One bubble with angle=pi/2
    elif choice==31: # Symmetric in columns! (Uncomment the symmetric function in columns)
        nx=256 # Number of cells per row
        ny=100 # Number of cells per column
        dx=1/256. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        radius=0.25
        index_circle=[]
        index_no_circle=[]
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if x[i]**2+y[j]**2<radius**2:
                    index_circle.append([i,j])
                    
                else:
                    index_no_circle.append([i,j])                      
        rho0=np.ones((nx,ny))-0.03
        for i,j in index_no_circle:
            rho0[i,j]=-0.97        
#        import matplotlib.pyplot as plt  
#        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
        epsilon=0.005 # Parameter epsilon   
        dt=0.001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=2 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=np.pi/2 # Angle for wetting boundary condition
        
        
        # One bubble with angle=pi/3
    elif choice==32: # Symmetric in columns! (Uncomment the symmetric function in columns)
        nx=256 # Number of cells per row
        ny=100 # Number of cells per column
        dx=1/256. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        radius=0.25
        index_circle=[]
        index_no_circle=[]
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if x[i]**2+y[j]**2<radius**2:
                    index_circle.append([i,j])
                    
                else:
                    index_no_circle.append([i,j])                      
        rho0=np.ones((nx,ny))-0.03
        for i,j in index_no_circle:
            rho0[i,j]=-0.97        
#        import matplotlib.pyplot as plt  
#        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
        epsilon=0.005 # Parameter epsilon   
        dt=0.001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=2 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=np.pi/3 # Angle for wetting boundary condition
        
        # One bubble with angle=2*pi/3
    elif choice==33: # Symmetric in columns! (Uncomment the symmetric function in columns)
        nx=256 # Number of cells per row
        ny=100 # Number of cells per column
        dx=1/256. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        radius=0.25
        index_circle=[]
        index_no_circle=[]
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if x[i]**2+y[j]**2<radius**2:
                    index_circle.append([i,j])
                    
                else:
                    index_no_circle.append([i,j])                      
        rho0=np.ones((nx,ny))-0.03
        for i,j in index_no_circle:
            rho0[i,j]=-0.97        
#        import matplotlib.pyplot as plt  
#        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
        epsilon=0.005 # Parameter epsilon   
        dt=0.001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=2 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=2*np.pi/3 # Angle for wetting boundary condition
        
      # One bubble with angle=5*pi/12
    elif choice==34: # Symmetric in columns! (Uncomment the symmetric function in columns)
        nx=256 # Number of cells per row
        ny=100 # Number of cells per column
        dx=1/256. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        radius=0.25
        index_circle=[]
        index_no_circle=[]
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if x[i]**2+y[j]**2<radius**2:
                    index_circle.append([i,j])
                    
                else:
                    index_no_circle.append([i,j])                      
        rho0=np.ones((nx,ny))-0.03
        for i,j in index_no_circle:
            rho0[i,j]=-0.97        
#        import matplotlib.pyplot as plt  
#        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
        epsilon=0.005 # Parameter epsilon   
        dt=0.001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=2 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=5*np.pi/12 # Angle for wetting boundary condition
        
        
      # One bubble with angle=7*pi/12
    elif choice==35: # Symmetric in columns! (Uncomment the symmetric function in columns)
        nx=256 # Number of cells per row
        ny=100 # Number of cells per column
        dx=1/256. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        radius=0.25
        index_circle=[]
        index_no_circle=[]
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if x[i]**2+y[j]**2<radius**2:
                    index_circle.append([i,j])
                    
                else:
                    index_no_circle.append([i,j])                      
        rho0=np.ones((nx,ny))-0.03
        for i,j in index_no_circle:
            rho0[i,j]=-0.97        
#        import matplotlib.pyplot as plt  
#        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
        epsilon=0.005 # Parameter epsilon   
        dt=0.001 # Time step
        tmax=0.1
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=2 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=7*np.pi/12 # Angle for wetting boundary condition
        
        
        
    #############################################
    # 4) Wetting phenomena from Aymard 2019 with two bubbles
    #############################################
        
    # Two bubble with angle=pi/4
    elif choice==41: # Symmetric in columns! (Uncomment the symmetric function in columns)
        nx=256 # Number of cells per row
        ny=64 # Number of cells per column
        dx=1/128. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        radius=0.30
        index_circle=[]
        index_no_circle=[]
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if ((x[i]-0.35)**2+y[j]**2<radius**2) or ((x[i]+0.35)**2+y[j]**2<radius**2):
                    index_circle.append([i,j])
                    
                else:
                    index_no_circle.append([i,j])                      
        rho0=np.ones((nx,ny))-0.03
        for i,j in index_no_circle:
            rho0[i,j]=-0.97      
        import matplotlib.pyplot as plt  
        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
        epsilon=0.012 # Parameter epsilon   
#        dt=0.001 # Time step
        dt=0.0005 # TAlternative time step
        tmax=15.

        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=2 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=np.pi/4 # Angle for wetting boundary condition
        
        
    # Two bubble with angle=3*pi/4
    elif choice==42: # Symmetric in columns! (Uncomment the symmetric function in columns)
        nx=256 # Number of cells per row
        ny=64 # Number of cells per column
        dx=1/128. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        radius=0.30
        index_circle=[]
        index_no_circle=[]
        # As initial condition we take a semi-circle centred at 0 and radius 0.25
        for i in range(nx):
            for j in range(ny):
                if ((x[i]-0.35)**2+y[j]**2<radius**2) or ((x[i]+0.35)**2+y[j]**2<radius**2):
                    index_circle.append([i,j])
                    
                else:
                    index_no_circle.append([i,j])                      
        rho0=np.ones((nx,ny))-0.03
        for i,j in index_no_circle:
            rho0[i,j]=-0.97      
        import matplotlib.pyplot as plt  
        plt.imshow(rho0)
        rho0=np.reshape(rho0,nx*ny)
        epsilon=0.012 # Parameter epsilon   
#        dt=0.001 # Time step
        dt=0.0005
#        tmax=15
        tmax=15
        ntimes=int(tmax/dt)# 4000 # Number of time steps        
        pot=1 # Choice of the potential: 1 is double-well, 2 is logarithmic
        theta=0.2 # Absolut temperature for the logarithmic potential
        thetac=1. # Critical temperature for the logarithmic potential
        choicemob=2 # Choice of mobility. 1 is constant mobility, 2 is 1-rho**2   
        bc=2 # Choice of boundary condition. 1 is no flux, 2 is wetting
        angle=3*np.pi/4 # Angle for wetting boundary condition
        
    
    
    return nx, ny, dx, x, y, rho0, epsilon, dt, tmax, ntimes, pot, theta, thetac, choicemob, bc, angle
        
        