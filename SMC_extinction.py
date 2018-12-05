# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 17:54:46 2018

@author: mugdhapolimera
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
smc = pd.read_csv('smcbar_ext.csv')
C1 = -4.959
C2 = 2.264
C3 = 0.389
C4 = 0.461 
x0 = 4.600
gamma = 1.000 

def F(x):
    if x >= 5.9:    
        return 0.5392*((x - 5.9)**2) + 0.05644*((x - 5.9)**3) 
    else:
        return 0
def D(x, gamma, x0):
    return x**2/((x**2 - x0**2)**2 + x**2 * gamma**2)

#wavelength = np.arange(1,0.125,0.001)
def smcextinction(x):        
    Rv = 2.75    
    if (x < 3.3):
        p = np.poly1d(np.polyfit(smc['#x'][smc['#x'] < 3.4],smc['Al/Av'][smc['#x'] < 3.4],1))
        return p(x)*Rv
    #elif( (x > 3.3) & (x < 3.7)):
     #   return (p(x) + (Exv_Ebv/Rv + 1))/2 #Ax_Av = 
    else:
        Exv_Ebv = C1 + C2*x + C3*D(x, gamma, x0) + C4*F(x)
        return ((Exv_Ebv/Rv + 1)*Rv)
#xline = np.arange(1,4,0.01)
'''
xlim = np.arange(1,8,0.01)
#x = 1.818
#print extinction(x)
#plt.figure()

for x in xlim:
   plt.plot(x,smcextinction(x)/2.75,'b.')
plt.plot(smc['#x'],smc['Al/Av'],'ro')
wavelength = np.arange(0.3,0.8001,0.0003)
Al_Ebv = np.array([smcextinction(1./wave) for wave in wavelength])
plt.plot(1/wavelength, Al_Ebv/2.75,'go')
#df = pd.DataFrame(np.column_stack((1/wavelength, 10000*wavelength, Al_Ebv)), columns = ['x', 'Wavelength', 'A_lambda/E(B-V)'])
#print df
#df.to_csv('smcbar_extinction.csv')
'''