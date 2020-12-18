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
# DEFINE FUNCTION: plot_imshow
##################

import matplotlib.pyplot as plt  
import numpy as np
 
def plot_imshow(x,rho,F,t,number):
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'


    #--------------------------------------------------------------------------
    # First plot: evolution of density in time
    #--------------------------------------------------------------------------
  
    figure1 = plt.figure(figsize=(7, 5.2))
    plt.clf()

    if number==11:
        plt.xlim(0, 30)
        plt.ylim(-40, 40)
        plt.xticks([0, 15, 30])
        plt.yticks([-30, 0, 30])
    elif number==12:
        plt.xlim(0, 30)
        plt.ylim(-40, 40)
        plt.xticks([0, 15, 30])
        plt.yticks([-30, 0, 30])
    
    plt.imshow(rho, cmap='Greys', vmin=-1, vmax=1,aspect='auto',extent=[t[0],t[-1],x[0],x[-1]]);
    cb = plt.colorbar(ticks=[-1, 0, 1])
    cb.ax.tick_params(labelsize=14)
    plt.ylabel(r'$x$', fontsize=23,rotation=0,labelpad=23)
    plt.xlabel(r'$t$', fontsize=23)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.show()
    plt.tight_layout()
    figure1.savefig('figures/density-'+str(number)+'.pdf', bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Second plot: free energy
    #--------------------------------------------------------------------------
    
    figure2 = plt.figure(figsize=(7, 5.2))
    plt.clf()
    plt.grid(True,linewidth=0.2)      
        
#    plt.xlim(-0.2, 3.2)
#    plt.ylim(1.5, 10.5)
#    plt.xticks([0.0,1.5,  3])
#    plt.yticks([4, 8,])    
    if number==11:
        plt.xlim(0, 30)
        plt.ylim(14, 21)
        plt.xticks([0,15,  30])
        plt.yticks([15,20])
    elif number==12:
        plt.xlim(0, 30)
        plt.ylim(25, 55)
        plt.xticks([0.0,15,  30])
        plt.yticks([30,50])
   
    ax = figure2.gca()
    ax.tick_params(axis="x", direction="in")
    
    plt.plot(t[10:],F[10:],linewidth=3,linestyle='-',color=(0,0.5647,0.6196));
    plt.ylabel(r'$\mathcal{F}[\phi]$', fontsize=23,rotation=0,labelpad=23)
    plt.xlabel(r'$t$', fontsize=23)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.show()
    plt.tight_layout()
    figure2.savefig('figures/freeenergy-'+str(number)+'.pdf', bbox_inches='tight')
    
    
    #--------------------------------------------------------------------------
    # Third plot: zoom in the density
    #--------------------------------------------------------------------------
    
    figure3 = plt.figure(figsize=(7, 5.2))
    plt.clf()
    plt.grid(True,linewidth=0.2)      
        

    if number==11:
        plt.xlim(-5, 5)
        plt.ylim(-1.1, 1.1)
        plt.xticks([-5,0,5])
        plt.yticks([-1,0,1])
    elif number==12:
        plt.xlim(-5, 5)
        plt.ylim(-1.1, 1.1)
        plt.xticks([-5,0,5])
        plt.yticks([-1,0,1])
   
    ax = figure3.gca()
    ax.tick_params(axis="x", direction="in")
    
    plt.plot(x[85:115],rho[85:115,-1],linewidth=3,linestyle='-',color=(0,0.5647,0.6196));
    plt.ylabel(r'$\phi$', fontsize=23,rotation=0,labelpad=23)
    plt.xlabel(r'$x$', fontsize=23)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.show()
    plt.tight_layout()
    figure3.savefig('figures/zoom-'+str(number)+'.pdf', bbox_inches='tight')
    
    
    