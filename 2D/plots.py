
###############################################################################
###############################################################################
#
# CODE TO SOLVE THE 2D CAHN-HILLIARD EQUATION
#
# AUTHOR OF THE CODE: SERGIO P. PEREZ
#
# FINITE VOLUMES, RUNGE-KUTTA TVD 3rd ORDER
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
sys.path.append(dirpath+'\\data')
from plots4times import *
from initial_conditions import *
from Euler_implicit_row import *
from Euler_implicit_column import *
from createplots import *
plt.rc('text',usetex=True)
plt.rc('font',family='serif')



##################
# DEFINE FUNCTION: MAIN
##################
    
def main():
    
    ##################
    # Select the configuration of the problem (see initial_conditions.py)
    
    choiceinitialcond=3
    ##################
    
    
    
    if choiceinitialcond==11:
    
        rho10=np.fromfile('data/rho_10_11.dat')
        rho10=np.reshape(rho10,(10,10,1001))
        rho10 = rho10[:,:,-1]
        rho20=np.fromfile('data/rho_20_11.dat')
        rho20=np.reshape(rho20,(20,20,1001))
        rho20 = rho20[:,:,-1]
        rho40=np.fromfile('data/rho_40_11.dat')
        rho40=np.reshape(rho40,(40,40,1001))
        rho40 = rho40[:,:,-1]
        rho80=np.fromfile('data/rho_80_11.dat')
        rho80=np.reshape(rho80,(80,80,1001))
        rho80 = rho80[:,:,-1]
#        rho160 = np.fromfile('data/rho_160_11.dat')
#        rho160=np.reshape(rho160,(160,160,20001))
#        rho160 = rho160[:,:,-1]
        
        x10= np.linspace(0,1,10)
        x20= np.linspace(0,1,20)
        x40= np.linspace(0,1,40)
        x80= np.linspace(0,1,80)
#        x160= np.linspace(0,1,160)
        
        exactrho80=np.zeros(80)
        for i in np.arange(80):
            if np.abs(np.sqrt((x80[i]-0.5)**2+(x80[40]-0.5)**2))<=np.pi*epsilon:
                exactrho80[i]=-1.+((np.cos(np.sqrt((x80[i]-0.5)**2+(x80[40]-0.5)**2)/epsilon)+1.)/np.pi)    
            else:
                exactrho80[i]=-1
        

    
        # Plot temporal evolution of the free energy    
        figure1 = plt.figure(figsize=(10, 6.5))
        plt.clf()
        plt.plot(x10,rho10[5,:],label=r'$10$');
        plt.plot(x20,rho20[10,:],label=r'$20$');
        plt.plot(x40,rho40[20,:],label=r'$40$');
        plt.plot(x80,rho80[40,:],label=r'$80$');
#        plt.plot(x80,exactrho80,label=r'$exact$');
        plt.plot(x80,-0.8407877686478619+0.13238269*np.cos((x80-0.5)/epsilon),label=r'$exact$');
#        plt.plot(x160,rho160[80,:],label=r'$160$');
        plt.ylabel(r'$\phi$', fontsize=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.legend(fontsize=20,framealpha=1, edgecolor='black')
        #plt.show()
        plt.tight_layout()
#        figure1.savefig('figures/freeenergy_%s.pdf'%choiceinitialcond, bbox_inches='tight')
        
    if choiceinitialcond==12:
    
        rho80=np.fromfile('data/rho_80_12.dat')
        rho80=np.reshape(rho80,(80,80,1001))

        x80= np.linspace(0,1,80)
#        x160= np.linspace(0,1,160)
        
        exactrho80=np.zeros(80)
        for i in np.arange(80):
            if np.abs((x80[i]-0.5))<=np.pi*epsilon:
                exactrho80[i]=-1.+((np.cos((x80[i]-0.5)/epsilon)+1.)/np.pi)    
            else:
                exactrho80[i]=-1
        

    
        # Plot temporal evolution of the free energy    
        figure1 = plt.figure(figsize=(10, 6.5))
        plt.clf()
        plt.plot(x80,rho80[:,40,-1],label=r'$80$');
        plt.plot(x80,exactrho80,label=r'$exact$');
#        plt.plot(x160,rho160[80,:],label=r'$160$');
        plt.ylabel(r'$\phi$', fontsize=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.legend(fontsize=20,framealpha=1, edgecolor='black')
        #plt.show()
        plt.tight_layout()
#        figure1.savefig('figures/freeenergy_%s.pdf'%choiceinitialcond, bbox_inches='tight')
        
        
        # Plot of initial conditions
          
        figure2 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho80[:,:,0].T, cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x80[0],x80[-1],x80[0],x80[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([0,0.5,1],fontsize=14)
        plt.yticks([0,0.5,1],fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        #plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-IC-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
        
        # Phase field at t=0.1
        
        figure3 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho80[:,:,-1].T, cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x80[0],x80[-1],x80[0],x80[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([0,0.5,1],fontsize=14)
        plt.yticks([0,0.5,1],fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        #plt.show()
        plt.tight_layout()
        figure3.savefig('figures/density-01-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
        
    elif choiceinitialcond==2:
        
        n=256 # Number of cells per row
        dx=1/256. # Width of cells
        x= np.linspace(-n*dx/2,n*dx/2,n) # Spatial mesh
        
        rho21=np.fromfile('data/rho_21.dat')
        rho22=np.fromfile('data/rho_22.dat')
        
        F21=np.fromfile('data/F_21.dat')
        F22=np.fromfile('data/F_22.dat')
        
        t21=np.fromfile('data/t_21.dat')
        t22=np.fromfile('data/t_22.dat')

        rho21=np.reshape(rho21,(n,n,len(t21)))
        rho22=np.reshape(rho22,(n,n,len(t22)))
        
#        # Plot temporal evolution of the free energy  
        figure1 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.grid(True,linewidth=0.2)      
        ax = figure1.gca()
        ax.tick_params(axis="x", direction="in")        
        plt.plot(t21,F21,linewidth=3,linestyle='-',color=(0,0.5647,0.6196),label=r'$M=1$');
        plt.plot(t22,F22,linewidth=3,linestyle='--',color='darkorange',label=r'$M=(1-\phi)(1+\phi)$');
        plt.ylabel(r'$\mathcal{F}[\phi]$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$t$', fontsize=23)
        plt.xticks([0,0.5,  1],fontsize=14)
        plt.yticks([15,40],fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(3, 52)
        plt.legend(fontsize=20,framealpha=1, edgecolor='black')
        #plt.show()
        plt.tight_layout()
        figure1.savefig('figures/freeenergy-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
    
    
        
        time_plot = 0.1
        
        index21 = np.argmin(np.abs(np.array(t21)-time_plot))
        index22 = np.argmin(np.abs(np.array(t22)-time_plot))
          
        figure2 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho21[:,:,index21], cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],x[0],x[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([-0.5,0,0.5],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        #plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-cons-01-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')        
    
        figure3 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho22[:,:,index22], cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],x[0],x[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([-0.5,0,0.5],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        #plt.show()
        plt.tight_layout()
        figure3.savefig('figures/density-deg-01-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
        
        
        time_plot = 0.4
        
        index21 = np.argmin(np.abs(np.array(t21)-time_plot))
        index22 = np.argmin(np.abs(np.array(t22)-time_plot))
          
        figure2 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho21[:,:,index21], cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],x[0],x[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([-0.5,0,0.5],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        #plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-cons-04-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')        
    
        figure3 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho22[:,:,index22], cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],x[0],x[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([-0.5,0,0.5],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        #plt.show()
        plt.tight_layout()
        figure3.savefig('figures/density-deg-04-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
        
        
        
        time_plot = 1
        
        index21 = np.argmin(np.abs(np.array(t21)-time_plot))
        index22 = np.argmin(np.abs(np.array(t22)-time_plot))
          
        figure2 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho21[:,:,index21], cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],x[0],x[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([-0.5,0,0.5],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        #plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-cons-1-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')        
    
        figure3 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.imshow(rho22[:,:,index22], cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],x[0],x[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([-0.5,0,0.5],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        #plt.show()
        plt.tight_layout()
        figure3.savefig('figures/density-deg-1-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
        
        
        
    elif choiceinitialcond==3:
        
        ny=100 # Number of cells per row
        nx=256
        dx=1/256. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        
        rho31=np.fromfile('data/rho_31.dat')     
        F31=np.fromfile('data/F_31.dat')       
        t31=np.fromfile('data/t_31.dat')
        
        rho32=np.fromfile('data/rho_32.dat')     
        F32=np.fromfile('data/F_32.dat')       
        t32=np.fromfile('data/t_32.dat')
        
        rho33=np.fromfile('data/rho_33.dat')     
        F33=np.fromfile('data/F_33.dat')       
        t33=np.fromfile('data/t_33.dat')
        
        rho34=np.fromfile('data/rho_34.dat')     
        F34=np.fromfile('data/F_34.dat')       
        t34=np.fromfile('data/t_34.dat')
        
        rho35=np.fromfile('data/rho_35.dat')     
        F35=np.fromfile('data/F_35.dat')       
        t35=np.fromfile('data/t_35.dat')

        rho31=np.reshape(rho31,(nx,ny,len(t31)))
        rho32=np.reshape(rho32,(nx,ny,len(t32)))
        rho33=np.reshape(rho33,(nx,ny,len(t33)))
        rho34=np.reshape(rho34,(nx,ny,len(t34)))
        rho35=np.reshape(rho35,(nx,ny,len(t35)))
        
        
        time_plot = 0.1
        
        index31 = np.argmin(np.abs(np.array(t31)-time_plot))
          
        figure1 = plt.figure(figsize=(5.12, 2))
        plt.clf()
        plt.imshow(np.flipud(rho31[:,:,index31].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([0,0.4],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(0, 0.4)
        #plt.show()
        plt.tight_layout()
        figure1.savefig('figures/density-31-2D.pdf', bbox_inches='tight') 
        
        index32 = np.argmin(np.abs(np.array(t32)-time_plot))
        figure2 = plt.figure(figsize=(5.12, 2))
        plt.clf()
        plt.imshow(np.flipud(rho32[:,:,index32].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([0,0.4],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(0, 0.4)
        #plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-32-2D.pdf', bbox_inches='tight') 
        
        index33 = np.argmin(np.abs(np.array(t33)-time_plot))
        figure3 = plt.figure(figsize=(5.12, 2))
        plt.clf()
        plt.imshow(np.flipud(rho33[:,:,index33].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([0,0.4],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(0, 0.4)
        #plt.show()
        plt.tight_layout()
        figure3.savefig('figures/density-33-2D.pdf', bbox_inches='tight') 
        
        index34 = np.argmin(np.abs(np.array(t34)-time_plot))
        figure4 = plt.figure(figsize=(5.12, 2))
        plt.clf()
        plt.imshow(np.flipud(rho34[:,:,index34].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([0,0.4],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(0, 0.4)
        #plt.show()
        plt.tight_layout()
        figure4.savefig('figures/density-34-2D.pdf', bbox_inches='tight') 
        
        index35 = np.argmin(np.abs(np.array(t35)-time_plot))
        figure5 = plt.figure(figsize=(5.12, 2))
        plt.clf()
        plt.imshow(np.flipud(rho35[:,:,index35].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-0.5,0,0.5],fontsize=14)
        plt.yticks([0,0.4],fontsize=14)
        plt.xlim(-0.5, 0.5)
        plt.ylim(0, 0.4)
        #plt.show()
        plt.tight_layout()
        figure5.savefig('figures/density-35-2D.pdf', bbox_inches='tight') 
 

        nx=256 # Number of cells per row
        ny=100 # Number of cells per column
        dt=0.001
        tmax=0.1
        epsilon=0.005
        ntimes=int(tmax/dt)# 4000 # Number of time steps 
        rho31=np.reshape(rho31,(nx*ny,ntimes+1))
        rho32=np.reshape(rho32,(nx*ny,ntimes+1))
        rho33=np.reshape(rho33,(nx*ny,ntimes+1))
        rho34=np.reshape(rho34,(nx*ny,ntimes+1))
        rho35=np.reshape(rho35,(nx*ny,ntimes+1))


        F31=dx**2*((0.25*(rho31[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho31,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho31,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho31,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho31,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        
        F32=dx**2*((0.25*(rho32[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho32,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho32,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho32,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho32,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        F33=dx**2*((0.25*(rho33[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho33,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho33,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho33,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho33,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        F34=dx**2*((0.25*(rho34[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho34,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho34,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho34,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho34,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        F35=dx**2*((0.25*(rho35[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho35,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho35,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho35,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho35,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        rho31=np.reshape(rho31,(nx,ny,ntimes+1))
        rho32=np.reshape(rho32,(nx,ny,ntimes+1))
        rho33=np.reshape(rho33,(nx,ny,ntimes+1))
        rho34=np.reshape(rho34,(nx,ny,ntimes+1))
        rho35=np.reshape(rho35,(nx,ny,ntimes+1))

        angle=np.pi/2
        F31=F31+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho31[:,0,:]**3/3-rho31[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho31[:,-1,:]**3/3-rho31[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho31[0,:,:]**3/3-rho31[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho31[-1,:,:]**3/3-rho31[-1,:,:])).sum(axis=0))
                
        angle=np.pi/3
        F32=F32+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho32[:,0,:]**3/3-rho32[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho32[:,-1,:]**3/3-rho32[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho32[0,:,:]**3/3-rho32[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho32[-1,:,:]**3/3-rho32[-1,:,:])).sum(axis=0))
                
        angle=2*np.pi/3
        F33=F33+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho33[:,0,:]**3/3-rho33[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho33[:,-1,:]**3/3-rho33[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho33[0,:,:]**3/3-rho33[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho33[-1,:,:]**3/3-rho33[-1,:,:])).sum(axis=0))
                
        angle=5*np.pi/12
        F34=F34+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho34[:,0,:]**3/3-rho34[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho34[:,-1,:]**3/3-rho34[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho34[0,:,:]**3/3-rho34[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho34[-1,:,:]**3/3-rho34[-1,:,:])).sum(axis=0))
                
        angle=7*np.pi/12
        F35=F35+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho35[:,0,:]**3/3-rho35[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho35[:,-1,:]**3/3-rho35[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho35[0,:,:]**3/3-rho35[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho35[-1,:,:]**3/3-rho35[-1,:,:])).sum(axis=0))
        
        

         # Plot temporal evolution of the free energy  
        figure6 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.grid(True,linewidth=0.2)      
        ax = figure6.gca()
        ax.tick_params(axis="x", direction="in")        
        plt.semilogx(t32[:-2],F32[2:],linewidth=3,linestyle='-',color=(0,0.5647,0.6196),label=r'$\beta=\pi/3$');
        plt.semilogx(t34[:-2],F34[2:],linewidth=3,linestyle='--',color='darkorange',label=r'$\beta=5\pi/12$');
        plt.semilogx(t31[:-1],F31[1:],linewidth=3,linestyle='-.',color='forestgreen',label=r'$\beta=\pi/2$');
        plt.semilogx(t35[:-2],F35[2:],linewidth=3,linestyle=(0, (1, 3)),color='rebeccapurple',label=r'$\beta=7\pi/12$');
        plt.semilogx(t33[:-2],F33[2:],linewidth=3,linestyle=(0, (1, 1)),color='orangered',label=r'$\beta=2\pi/3$');
        plt.ylabel(r'$\mathcal{F}[\phi]$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$t$', fontsize=23)
        plt.xticks([0,0.001,0.01, 0.1],fontsize=14)
        plt.yticks([0.003,0.005,0.007],fontsize=14)
        plt.xlim(0, 0.1)
        plt.ylim(0.0015,0.0105)
        plt.legend(fontsize=20,framealpha=1, edgecolor='black',ncol=2)
        #plt.show()
        plt.tight_layout()
        figure6.savefig('figures/freeenergy-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
    
    elif choiceinitialcond==31: # These plots are created for the measurement of the angle
        
        ny=100 # Number of cells per row
        nx=256
        dx=1/256. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        
        thres = 0
        rho31=np.fromfile('data/rho_31.dat')  
        rho31[rho31>thres]=1
        rho31[rho31<=thres]=-1
        F31=np.fromfile('data/F_31.dat')       
        t31=np.fromfile('data/t_31.dat')
        
        rho32=np.fromfile('data/rho_32.dat') 
        rho32[rho32>thres]=1
        rho32[rho32<=thres]=-1
        F32=np.fromfile('data/F_32.dat')       
        t32=np.fromfile('data/t_32.dat')
        
        rho33=np.fromfile('data/rho_33.dat')  
        rho33[rho33>thres]=1
        rho33[rho33<=thres]=-1
        F33=np.fromfile('data/F_33.dat')       
        t33=np.fromfile('data/t_33.dat')
        
        rho34=np.fromfile('data/rho_34.dat')
        rho34[rho34>thres]=1
        rho34[rho34<=thres]=-1
        F34=np.fromfile('data/F_34.dat')       
        t34=np.fromfile('data/t_34.dat')
        
        rho35=np.fromfile('data/rho_35.dat')  
        rho35[rho35>thres]=1
        rho35[rho35<=thres]=-1
        F35=np.fromfile('data/F_35.dat')       
        t35=np.fromfile('data/t_35.dat')

        rho31=np.reshape(rho31,(nx,ny,len(t31)))
        rho32=np.reshape(rho32,(nx,ny,len(t32)))
        rho33=np.reshape(rho33,(nx,ny,len(t33)))
        rho34=np.reshape(rho34,(nx,ny,len(t34)))
        rho35=np.reshape(rho35,(nx,ny,len(t35)))
        
        
        time_plot = 0.1
        
        index31 = np.argmin(np.abs(np.array(t31)-time_plot))

        
        figure1 = plt.figure(figsize=(51.2, 20),frameon=False)
        plt.clf()
        plt.imshow(np.flipud(rho31[:,:,index31].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        plt.axis('off')
#        cb = plt.colorbar(ticks=[-1, 0, 1])
#        cb.ax.tick_params(labelsize=14)
#        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
#        plt.xlabel(r'$x$', fontsize=23)
#        plt.xticks([-0.5,0,0.5],fontsize=14)
#        plt.yticks([0,0.4],fontsize=14)
#        plt.xlim(-0.5, 0.5)
#        plt.ylim(0, 0.4)
        #plt.show()
#        plt.tight_layout()
        figure1.savefig('figures/density-31-2D.png',pad_inches=0) 
        
        index32 = np.argmin(np.abs(np.array(t32)-time_plot))
        figure2 = plt.figure(figsize=(51.2, 20),frameon=False)
        plt.clf()
        plt.imshow(np.flipud(rho32[:,:,index32].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        plt.axis('off')
#        cb = plt.colorbar(ticks=[-1, 0, 1])
#        cb.ax.tick_params(labelsize=14)
#        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
#        plt.xlabel(r'$x$', fontsize=23)
#        plt.xticks([-0.5,0,0.5],fontsize=14)
#        plt.yticks([0,0.4],fontsize=14)
#        plt.xlim(-0.5, 0.5)
#        plt.ylim(0, 0.4)
        #plt.show()
#        plt.tight_layout()
        figure2.savefig('figures/density-32-2D.png',pad_inches=0) 
        
        figure3 = plt.figure(figsize=(51.2, 20),frameon=False)
        plt.clf()
        plt.imshow(np.flipud(rho33[:,:,index31].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        plt.axis('off')
#        cb = plt.colorbar(ticks=[-1, 0, 1])
#        cb.ax.tick_params(labelsize=14)
#        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
#        plt.xlabel(r'$x$', fontsize=23)
#        plt.xticks([-0.5,0,0.5],fontsize=14)
#        plt.yticks([0,0.4],fontsize=14)
#        plt.xlim(-0.5, 0.5)
#        plt.ylim(0, 0.4)
        #plt.show()
#        plt.tight_layout()
        figure3.savefig('figures/density-33-2D.png',pad_inches=0) 
        
        figure4 = plt.figure(figsize=(51.2, 20),frameon=False)
        plt.clf()
        plt.imshow(np.flipud(rho34[:,:,index31].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        plt.axis('off')
#        cb = plt.colorbar(ticks=[-1, 0, 1])
#        cb.ax.tick_params(labelsize=14)
#        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
#        plt.xlabel(r'$x$', fontsize=23)
#        plt.xticks([-0.5,0,0.5],fontsize=14)
#        plt.yticks([0,0.4],fontsize=14)
#        plt.xlim(-0.5, 0.5)
#        plt.ylim(0, 0.4)
        #plt.show()
#        plt.tight_layout()
        figure4.savefig('figures/density-34-2D.png',pad_inches=0) 
        
        figure5 = plt.figure(figsize=(51.2, 20),frameon=False)
        plt.clf()

        plt.imshow(np.flipud(rho35[:,:,index31].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        plt.axis('off')
#        cb = plt.colorbar(ticks=[-1, 0, 1])
#        cb.ax.tick_params(labelsize=14)
#        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
#        plt.xlabel(r'$x$', fontsize=23)
#        plt.xticks([-0.5,0,0.5],fontsize=14)
#        plt.yticks([0,0.4],fontsize=14)
#        plt.xlim(-0.5, 0.5)
#        plt.ylim(0, 0.4)
        #plt.show()
#        plt.tight_layout()
        figure5.savefig('figures/density-35-2D.png',pad_inches=0) 
        
        

       
    
    elif choiceinitialcond==4:
        
        ny=64 # Number of cells per row
        nx=256
        dx=1/128. # Width of cells
        x= np.linspace(-nx*dx/2,nx*dx/2,nx) # Spatial mesh 1D
        y= np.linspace(0,ny*dx,ny)
        
        rho41=np.fromfile('data/rho_41.dat')     
        F41=np.fromfile('data/F_41.dat')       
        t41=np.fromfile('data/t_41.dat')
        
        rho42=np.fromfile('data/rho_42.dat')     
        F42=np.fromfile('data/F_42.dat')       
        t42=np.fromfile('data/t_42.dat')
        
        
        

        rho41=np.reshape(rho41,(nx,ny,len(t41)))
        rho42=np.reshape(rho42,(nx,ny,len(t42)))
        
        time_plot =0.1
        
        index41 = np.argmin(np.abs(np.array(t41)-time_plot))
        index42 = np.argmin(np.abs(np.array(t42)-time_plot))
          
        figure1 = plt.figure(figsize=(6, 3))
        plt.clf()
        plt.imshow(np.flipud(rho41[:,:,index41].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-1,0,1],fontsize=14)
        plt.yticks([0,0.5],fontsize=14)
        plt.xlim(-1, 1)
        plt.ylim(0, 0.5)
        #plt.show()
        plt.tight_layout()
        figure1.savefig('figures/density-pi4-01-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight') 
        
        figure2 = plt.figure(figsize=(6, 3))
        plt.clf()
        plt.imshow(np.flipud(rho42[:,:,index42].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-1,0,1],fontsize=14)
        plt.yticks([0,0.5],fontsize=14)
        plt.xlim(-1, 1)
        plt.ylim(0, 0.5)
#        plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-3pi4-01-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight') 
        
        
        time_plot =1
        
        index41 = np.argmin(np.abs(np.array(t41)-time_plot))
        index42 = np.argmin(np.abs(np.array(t42)-time_plot))
          
        figure1 = plt.figure(figsize=(6, 3))
        plt.clf()
        plt.imshow(np.flipud(rho41[:,:,index41].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-1,0,1],fontsize=14)
        plt.yticks([0,0.5],fontsize=14)
        plt.xlim(-1, 1)
        plt.ylim(0, 0.5)
        #plt.show()
        plt.tight_layout()
#        figure1.savefig('figures/density-pi4-1-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight') 
        
        figure2 = plt.figure(figsize=(6, 3))
        plt.clf()
        plt.imshow(np.flipud(rho42[:,:,index42].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-1,0,1],fontsize=14)
        plt.yticks([0,0.5],fontsize=14)
        plt.xlim(-1, 1)
        plt.ylim(0, 0.5)
#        plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-3pi4-1-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight') 
        
        
        time_plot =15
        
        index41 = np.argmin(np.abs(np.array(t41)-time_plot))
        index42 = np.argmin(np.abs(np.array(t42)-time_plot))
          
        figure1 = plt.figure(figsize=(6, 3))
        plt.clf()
        plt.imshow(np.flipud(rho41[:,:,index41].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-1,0,1],fontsize=14)
        plt.yticks([0,0.5],fontsize=14)
        plt.xlim(-1, 1)
        plt.ylim(0, 0.5)
        #plt.show()
        plt.tight_layout()
        figure1.savefig('figures/density-pi4-15-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight') 
        
        figure2 = plt.figure(figsize=(6, 3))
        plt.clf()
        plt.imshow(np.flipud(rho42[:,:,index42].T), cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[x[0],x[-1],y[0],y[-1]]);
        cb = plt.colorbar(ticks=[-1, 0, 1])
        cb.ax.tick_params(labelsize=14)
        plt.ylabel(r'$y$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$x$', fontsize=23)
        plt.xticks([-1,0,1],fontsize=14)
        plt.yticks([0,0.5],fontsize=14)
        plt.xlim(-1, 1)
        plt.ylim(0, 0.5)
#        plt.show()
        plt.tight_layout()
        figure2.savefig('figures/density-3pi4-15-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight') 
        

        nx=256 # Number of cells per row
        ny=64 # Number of cells per column
        dt=0.0005
        tmax=15
        epsilon=0.012
        ntimes=int(tmax/dt)# 4000 # Number of time steps 
        rho41=np.reshape(rho41,(nx*ny,ntimes+1))
        rho42=np.reshape(rho42,(nx*ny,ntimes+1))

        F41=dx**2*((0.25*(rho41[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho41,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho41,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho41,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho41,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        
        F42=dx**2*((0.25*(rho42[:,:]**2-1)**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho42,(nx,ny,ntimes+1))[1:,:,:]-np.reshape(rho42,(nx,ny,ntimes+1))[:-1,:,:])**2,((nx-1)*ny,ntimes+1))/dx**2).sum(axis=0)\
        +(epsilon**2/2*np.reshape((np.reshape(rho42,(nx,ny,ntimes+1))[:,1:,:]-np.reshape(rho42,(nx,ny,ntimes+1))[:,:-1,:])**2,((ny-1)*nx,ntimes+1))/dx**2).sum(axis=0))
        
        rho41=np.reshape(rho41,(nx,ny,ntimes+1))
        rho42=np.reshape(rho42,(nx,ny,ntimes+1))
        angle=np.pi/4
        F41=F41+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho41[:,0,:]**3/3-rho41[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho41[:,-1,:]**3/3-rho41[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho41[0,:,:]**3/3-rho41[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho41[-1,:,:]**3/3-rho41[-1,:,:])).sum(axis=0))
                
        angle=3*np.pi/4
        F42=F42+dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho42[:,0,:]**3/3-rho42[:,0,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho42[:,-1,:]**3/3-rho42[:,-1,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho42[0,:,:]**3/3-rho42[0,:,:])).sum(axis=0))\
                +dx*epsilon*((np.sqrt(2)/2*np.cos(angle)*(rho42[-1,:,:]**3/3-rho42[-1,:,:])).sum(axis=0))
        
        
        # Plot temporal evolution of the free energy  
        figure6 = plt.figure(figsize=(7, 5.2))
        plt.clf()
        plt.grid(True,linewidth=0.2)      
        ax = figure6.gca()
        ax.tick_params(axis="x", direction="in")        
        plt.semilogx(t41[:-1],F41[1:],linewidth=3,linestyle='-',color=(0,0.5647,0.6196),label=r'$\beta=\pi/4$');
        plt.semilogx(t42[:-1],F42[1:],linewidth=3,linestyle='--',color='darkorange',label=r'$\beta=3\pi/4$');
        plt.ylabel(r'$\mathcal{F}[\phi]$', fontsize=23,rotation=0,labelpad=23)
        plt.xlabel(r'$t$', fontsize=23)
        plt.xticks([0,0.001, 0.01, 0.1, 1, 10],fontsize=14)
        plt.yticks([0.01,0.04],fontsize=14)
        plt.xlim(0., 15)
        plt.ylim(0, 0.05)
        plt.legend(fontsize=20,framealpha=1, edgecolor='black')
        #plt.show()
        plt.tight_layout()
        figure6.savefig('figures/freeenergy-%s-2D.pdf'%choiceinitialcond, bbox_inches='tight')
        
        
        
        
    
    
    
    

    
    
    
    
if __name__ == '__main__':
    main()