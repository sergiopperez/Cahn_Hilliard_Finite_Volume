# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:00:56 2019

@author: sp3215
"""
import numpy as np
from plots4times import *

def createplots(x,rho,F,dFdrho,t,choiceinitialcond,theta,choicemob):

    if choiceinitialcond==2:
            
            index1=np.argmin(np.abs(t[:]-0.01))
            index2=np.argmin(np.abs(t[:]-0.02))
            index3=np.argmin(np.abs(t[:]-0.1))
            timespos=[0,index1,index2,index3]
            times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
            plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
            
    elif choiceinitialcond==3:
        
        if theta==0.3 and choicemob==1:
        
            index1=np.argmin(np.abs(t[:]-0.01))
            index2=np.argmin(np.abs(t[:]-0.02))
            index3=np.argmin(np.abs(t[:]-0.1))
            timespos=[0,index1,index2,index3]
            times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
            plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
            
        elif theta==0.3 and choicemob==2:
        
            index1=np.argmin(np.abs(t[:]-0.01))
            index2=np.argmin(np.abs(t[:]-0.02))
            index3=np.argmin(np.abs(t[:]-0.1))
            timespos=[0,index1,index2,index3]
            times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
            plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
        
        else:
            
        
            index1=np.argmin(np.abs(t[:]-0.01))
            index2=np.argmin(np.abs(t[:]-0.02))
            index3=np.argmin(np.abs(t[:]-0.1))
            timespos=[0,index1,index2,index3]
            times=['$$t=0$$   ', '$$t=0.01$$ ' ,'$$t=0.02$$ ', '$$t=0.1$$']
            plots4times(x,rho,F,dFdrho,t,times,timespos[0],timespos[1],timespos[2],timespos[3],choiceinitialcond)
        


