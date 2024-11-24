# imports

import numpy
import scipy
from matplotlib import pyplot

# key variables
#geometric
Rc = 33e-3
Rt = 13e-3
Re = 34e-3
R1 = 62e-3
R2 = 20.9e-3
Lcyl = 66e-3
Lc = 120e-3
Le = 54e-3
L = Lc + Le
#material
J=200
W=100

Thta = [] # empty state
# R3, T3, & Tn not used for approximation (might include later)


g1 = lambda y : Rc
g2 = lambda y : Rc - R1 + (R1**2 - (y - Lcyl)**2)**(1/2)
g3 = lambda y : Rt + R2 - (R2**2 - (y - Lc)**2)**(1/2)
g4 = lambda y : (Re - Rt)/Le * (y - Lc) + Rt


# function

def mesh(n, L_total, state_vector): # n: number of elements -> n + 1 points
    y = []
    R = []
    
    for i in range(n+1):
        state_vector.append(300) # room temp
        y.append(i/n * L)
        if y[i] < Lcyl:
            R.append(g1(y[i]))
        elif y[i] < 106.386: # this value needs automation
            R.append(g2(y[i]))
        elif y[i] < Lc:
            R.append(g2(y[i]))
        elif y[i] < L:
            R.append(g2(y[i]))
            
    
    
    

