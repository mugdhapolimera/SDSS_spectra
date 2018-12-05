# -*- coding: utf-8 -*-
"""
Created on Tue Dec 04 18:40:53 2018

@author: mugdhapolimera
"""

def F(x):
    if x >= 5.9:    
        return 0.5392*((x - 5.9)**2) + 0.05644*((x - 5.9)**3) 
    else:
        return 0
def D(x, gamma, x0):
    return x**2/((x**2 - x0**2)**2 + x**2 * gamma**2)

#wavelength = np.arange(1,0.125,0.001)
def smcextinction(x):        
    C1 = -4.959
    C2 = 2.264
    C3 = 0.389
    C4 = 0.461 
    x0 = 4.600
    gamma = 1.000 
    Rv = 2.75    

    if (x < 3.3):
        Al_Av = 0.69383845*(x**2) - 0.29124909*x
    #elif( (x > 3.3) & (x < 3.7)):
     #   return (p(x) + (Exv_Ebv/Rv + 1))/2 #Ax_Av = 
    else:
        Ex_Ebv = C1 + C2*x + C3*D(x, gamma, x0) + C4*F(x)
        Al_Av = Ex_Ebv/Rv + 1
    return (Al_Av*Rv)
        
def mwextinction(x):        
    Rv = 3.07   
    y = x - 1.82
    a = 1.0 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 - 1.718*y**5 \
    - 0.827*y**6 + 1.647*y**7 - 0.505*y**8        
    b = 1.952*y + 2.908*(y**2) - 3.989*(y**3) - 7.985*(y**4) + 11.102*(y**5) \
    + 5.491*(y**6) - 10.805*(y**7) + 3.347*(y**8)        
                
    Al_Av = a + b/Rv
    return (Al_Av*Rv)
