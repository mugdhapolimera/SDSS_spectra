# -*- coding: utf-8 -*-
"""
Created on Mon Dec 03 11:40:26 2018

@author: mugdhapolimera
"""
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 17:54:46 2018

@author: mugdhapolimera
"""
#import numpy as np
#import matplotlib.pyplot as plt
#import pandas as pd

def mwextinction(x):        
    Rv = 3.07   
    y = x - 1.82
    a = 1.0 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 - 1.718*y**5 \
    - 0.827*y**6 + 1.647*y**7 - 0.505*y**8        
    b = 1.952*y + 2.908*(y**2) - 3.989*(y**3) - 7.985*(y**4) + 11.102*(y**5) \
    + 5.491*(y**6) - 10.805*(y**7) + 3.347*(y**8)        
                
    Al_Av = a + b/Rv
    return (Al_Av*Rv)
#xline = np.arange(1,4,0.01)
#xlim = np.arange(1,8,0.01)
#x = 1.818
# Cardelli, Clayton & Mathis (1989) with A_V = 1 and R_V = 3.1
#plt.plot(mw['Wavelength'],mw['Al_Ebv'],'ro')
#wavelength = np.arange(3000,8000,3, dtype = float)
#Al_Ebv = np.array([mwextinction(float(1e4/wave)) for wave in wavelength])
#plt.plot(wavelength, Al_Ebv,'go')
#plt.plot(wavelength,ext*3.1, 'bo')
#df = pd.DataFrame(np.column_stack((1/wavelength, 10000*wavelength, Al_Ebv)), columns = ['x', 'Wavelength', 'A_lambda/E(B-V)'])
#print df
#df.to_csv('smcbar_extinction.csv')
