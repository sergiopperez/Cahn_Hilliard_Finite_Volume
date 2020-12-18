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
# DEFINE FUNCTION: createplots
##################

import numpy as np
from plots4times import *
from plot_imshow import *
from scipy import interpolate

def createplots(x,rho,F,dFdrho,t,choiceinitialcond,theta,choicemob):

    if choiceinitialcond==11 or choiceinitialcond==12:
        
       plot_imshow(x,rho,F,t,choiceinitialcond)
       
       
    
    elif choiceinitialcond==2:
            
        index1=np.argmin(np.abs(t[:]-0.01))
        index2=np.argmin(np.abs(t[:]-0.02))
        index3=np.argmin(np.abs(t[:]-0.1))
        timespos=[0,index1,index2,index3]
        times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
        tck = interpolate.splrep(x, dFdrho[:,0], s=0.000001)
        dFdrho[:,0]=interpolate.splev(x, tck, der=0)
        plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
            
    elif choiceinitialcond==31: #theta==0.3 and choicemob==1
    
        index1=np.argmin(np.abs(t[:]-0.01))
        index2=np.argmin(np.abs(t[:]-0.02))
        index3=np.argmin(np.abs(t[:]-0.1))
        timespos=[0,index1,index2,index3]
        times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
        plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
            
    elif choiceinitialcond==32: #theta==0.3 and choicemob==2
        
        index1=np.argmin(np.abs(t[:]-0.02))
        index2=np.argmin(np.abs(t[:]-0.1))
        index3=np.argmin(np.abs(t[:]-2))
        timespos=[0,index1,index2,index3]
        times=['$$t=0$$   ', '$$t=0.02$$ ' ,'$$t=0.1$$ ', '$$t=2$$']
        plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
        
    elif choiceinitialcond==33: #theta==0.0 and choicemob==1            
        
        index1=np.argmin(np.abs(t[:]-0.01))
        index2=np.argmin(np.abs(t[:]-0.02))
        index3=np.argmin(np.abs(t[:]-2))
        timespos=[0,index1,index2,index3]
        times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=2$$ ']
        plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
        
    elif choiceinitialcond==34: #theta==0.0 and choicemob==2            
        
        index1=np.argmin(np.abs(t[:]-0.01))
        index2=np.argmin(np.abs(t[:]-0.02))
        index3=np.argmin(np.abs(t[:]-0.1))
        timespos=[0,index1,index2,index3]
        times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
        dFdrho[np.isnan(dFdrho)==1]=dFdrho[0,index3]
        #dFdrho[indexnan]=dFdrho[0,index3]
        
        plots4times(x,rho,F[:index1],dFdrho,t[:index1],times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
    
    elif choiceinitialcond==35: # Double-well and choicemob==1            
        
        index1=np.argmin(np.abs(t[:]-0.01))
        index2=np.argmin(np.abs(t[:]-0.02))
        index3=np.argmin(np.abs(t[:]-0.1))
        timespos=[0,index1,index2,index3]
        times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
              
        plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
    
    elif choiceinitialcond==36: # Double-well and choicemob==2            
        
        index1=np.argmin(np.abs(t[:]-0.01))
        index2=np.argmin(np.abs(t[:]-0.02))
        index3=np.argmin(np.abs(t[:]-0.1))
        timespos=[0,index1,index2,index3]
        times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
       
        plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
    

