import matplotlib.pyplot as plt  
    
def plots4times(x,rho,F,dFdrho,t,times,t1,t2,t3,t4,number):
    
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
    plt.grid(True,linewidth=0.7)
    
    #plt.xlim(0, 0.2)
    #plt.ylim(0, 0.2)
    #plt.xticks(np.arange(0, 1, step=0.2))
    #plt.yticks(np.arange(0, 1, step=0.2))
    
    
    plt.plot(x,Y1,label=times[0],linewidth=3,linestyle=':',color=(0,0.5647,0.6196));
    
    plt.plot(x,Y2,label=times[1],linewidth=3,linestyle='--',color=(0.4980,0.6902,0.0196));
    
    plt.plot(x,Y3,label=times[2],linewidth=3,linestyle='-.',color=(0.87058824300766,0.490196079015732,0));
    
    plt.plot(x,Y4,label=times[3],linewidth=3,linestyle='-',color=(0.4863,0.0784,0.3020));
    
    
    #plt.plot(x,exactrho);
    plt.ylabel(r'$\rho$', fontsize=23,rotation=0,labelpad=23)
    plt.xlabel(r'$x$', fontsize=23)
    plt.legend(fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.show()
    plt.tight_layout()
    figure1.savefig('figures/densityprofile-'+str(number)+'.pdf', bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Second plot: free energy
    #--------------------------------------------------------------------------
    
    figure2 = plt.figure(figsize=(7, 5.2))
    plt.clf()
    plt.grid(True,linewidth=0.7)
    
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
    plt.grid(True,linewidth=0.7)
    
    # Extract the density at four different times from the matrix rho
    Y1=dFdrho[:,t1];
    Y2=dFdrho[:,t2];
    Y3=dFdrho[:,t3];
    Y4=dFdrho[:,t4];
    
    plt.clf()
    plt.grid(True,linewidth=0.7)
    
    plt.plot(x,Y1,label=times[0],linewidth=3,linestyle=':',color=(0,0.5647,0.6196));
    
    plt.plot(x,Y2,label=times[1],linewidth=3,linestyle='--',color=(0.4980,0.6902,0.0196));
    
    plt.plot(x,Y3,label=times[2],linewidth=3,linestyle='-.',color=(0.87058824300766,0.490196079015732,0));
    
    plt.plot(x,Y4,label=times[3],linewidth=3,linestyle='-',color=(0.4863,0.0784,0.3020));
    
    
    plt.ylabel(r'$\frac{\delta \mathcal{F}}{\delta \rho}$', fontsize=23,rotation=0,labelpad=23)
    plt.xlabel(r'$t$', fontsize=23)
    plt.legend(fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.show()
    plt.tight_layout()
    figure3.savefig('figures/varfreeenergy-'+str(number)+'.pdf', bbox_inches='tight')
    