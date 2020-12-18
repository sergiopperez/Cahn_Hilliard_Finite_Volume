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
# DEFINE FUNCTION: plots4times
##################

import matplotlib.pyplot as plt  
import numpy as np
 
def plots4times(x,rho,F,dFdrho,t,times,t1,t2,t3,t4,number):
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    #--------------------------------------------------------------------------
    # First plot: density profile at four different times
    #--------------------------------------------------------------------------
    
    # Extract the density at four different times from the matrix rho
    Y1=rho[:,t1];
    Y2=rho[:,t2];
    Y3=rho[:,t3];
    Y4=rho[:,t4];
    
    
    figure1 = plt.figure(figsize=(7, 5.2))
    plt.clf()
    plt.grid(True,linewidth=0.2)
#    ax = figure1.gca()
#    ax.tick_params(axis="x", direction="in")
    
    if number==2:
        plt.xlim(-0.05, 1.05)
        plt.ylim(-1.05, 0.05)
        plt.xticks([0, 0.5, 1])
        plt.yticks([0, -0.5, -1])
    elif number==31 or number==32 or number==34 or number==35 or number==36:
        plt.xlim(-0.05, 1.05)
        plt.ylim(-1.1, 1.1)
        plt.xticks([0, 0.5, 1])
        plt.yticks([-1, 0, 1])

    
    
    plt.plot(x,Y1,label=times[0],linewidth=3,linestyle=':',color=(0,0.5647,0.6196));
    
    plt.plot(x,Y2,label=times[1],linewidth=3,linestyle='--',color=(0.4980,0.6902,0.0196));
    
    plt.plot(x,Y3,label=times[2],linewidth=3,linestyle='-.',color=(0.87058824300766,0.490196079015732,0));
    
    plt.plot(x,Y4,label=times[3],linewidth=3,linestyle='-',color=(0.4863,0.0784,0.3020));
    
    
    #plt.plot(x,exactrho);
    plt.ylabel(r'$\phi$', fontsize=23,rotation=0,labelpad=23)
    plt.xlabel(r'$x$', fontsize=23)
    plt.legend(fontsize=20,framealpha=1, edgecolor='black')
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
    
    if number==2:
#        plt.xlim(0, 0.1)
        plt.ylim(0.156, 0.204)
        plt.xticks([0.001, 0.01, 0.1])
        plt.yticks([0.2, 0.18, 0.16])
        t=t[1:]-(t[1]-t[0])
        F=F[1:]
    elif number==31:
        plt.xlim(-0.005, 0.105)
        plt.ylim(0.025, 0.075)
        plt.xticks([0.0,0.05,  0.1])
        plt.yticks([0.04, 0.06])
    elif number==32:
        plt.xlim(-0.1, 2.1)
        plt.ylim(0.025, 0.075)
        plt.xticks([0.0,1.0,  2.0])
        plt.yticks([0.04, 0.06])
    elif number==34:
        plt.xlim(-0.0005, 0.0105)
        plt.ylim(0.088, 0.104)
        plt.xticks([0.0,0.005,  0.01])
        plt.yticks([0.092, 0.1])
    elif number==35:
        plt.xlim(-0.005, 0.105)
        plt.ylim(0.03, 0.065)
        plt.xticks([0.0,0.05,  0.1])
        plt.yticks([0.04, 0.055])
    elif number==36:
        plt.xlim(-0.01, 0.21)
        plt.ylim(0.03, 0.065)
        plt.xticks([0.0,0.1,  0.2])
        plt.yticks([0.04, 0.055])
        
        
        
    ax = figure2.gca()
    ax.tick_params(axis="x", direction="in")
    
    if number==2:
        plt.semilogx(t,F,linewidth=3,linestyle='-',color=(0,0.5647,0.6196));
    else:
        plt.plot(t,F,linewidth=3,linestyle='-',color=(0,0.5647,0.6196));
    plt.ylabel(r'$\mathcal{F}[\phi]$', fontsize=23,rotation=0,labelpad=23)
    plt.xlabel(r'$t$', fontsize=23)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.show()
    plt.tight_layout()
    figure2.savefig('figures/freeenergy-'+str(number)+'.pdf', bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Third plot: variation of free energy
    #--------------------------------------------------------------------------
    
    figure3 = plt.figure(figsize=(7, 5.2))
    plt.clf()
    plt.grid(True,linewidth=0.2)
    
    if number==2:
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)
        plt.xticks([0, 0.5, 1])
        plt.yticks([0, 0.5, 1])
    elif number==31 or number==32:
        plt.xlim(-0.05, 1.05)
        plt.ylim(-1, 1)
        plt.xticks([0, 0.5, 1])
        plt.yticks([-0.75, 0, 0.75])        
    elif number==34:
        plt.xlim(-0.05, 1.05)
        plt.ylim(-1.1, 1.1)
        plt.xticks([0, 0.5, 1])
        plt.yticks([-1, 0, 1])
    elif number==35:
        plt.xlim(-0.05, 1.05)
        plt.ylim(-1.7, 1.7)
        plt.xticks([0, 0.5, 1])
        plt.yticks([-1.5, 0, 1.5])
    elif number==36:
        plt.xlim(-0.05, 1.05)
        plt.ylim(-1.7, 1.7)
        plt.xticks([0, 0.5, 1])
        plt.yticks([-1.5, 0, 1.5])
  
    
    # Extract the density at four different times from the matrix rho
    Y1=dFdrho[:,t1];
    Y2=dFdrho[:,t2];
    Y3=dFdrho[:,t3];
    Y4=dFdrho[:,t4];

    plt.plot(x,Y1,label=times[0],linewidth=3,linestyle=':',color=(0,0.5647,0.6196));
    
    plt.plot(x,Y2,label=times[1],linewidth=3,linestyle='--',color=(0.4980,0.6902,0.0196));
    
    plt.plot(x,Y3,label=times[2],linewidth=3,linestyle='-.',color=(0.87058824300766,0.490196079015732,0));
    
    plt.plot(x,Y4,label=times[3],linewidth=3,linestyle='-',color=(0.4863,0.0784,0.3020));
    
    
    plt.ylabel(r'$\frac{\delta \mathcal{F}}{\delta \rho}$', fontsize=27,rotation=0,labelpad=23)
    plt.xlabel(r'$x$', fontsize=23)
    if number==32 or number==34:
        plt.legend(fontsize=20,framealpha=1, edgecolor='black',loc='lower right')
    else:
        plt.legend(fontsize=20,framealpha=1, edgecolor='black')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.show()
    plt.tight_layout()
    figure3.savefig('figures/varfreeenergy-'+str(number)+'.pdf', bbox_inches='tight')
    