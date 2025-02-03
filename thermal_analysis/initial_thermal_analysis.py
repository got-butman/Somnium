import numpy as np
import matplotlib.pyplot as plt

# Constants
N_layers = 30 #number of layers

gamma_ox= 7/5 #adiabatic index for diatomic gas
gamma_ipa = 4/3 #adiabatic index for non linear poliatomic gas

ox_ipa_ration=1.86 #ratio between oxygen and IPA

gamma_mix = (ox_ipa_ration*gamma_ox+gamma_ipa)/(ox_ipa_ration+1) #adiabatic index for the mixture

#initial conditions
T_initial = 2000 #K         DONT KNOW!!
R1_Rt_ratio=1.5
R2_Rt_ratio=0.382
R1=0.5 #m  DONT KNOW!! 
Rt = R1/R1_Rt_ratio
contr_angle=30 #degrees betweeen chamber and throat
exp_angle=15 #degrees between throat and exit       DONT KNOW!!

def chambertothroat_radial_profile(R1, Rt, contr_angle): #calculate radii of different layers

    #initial radial profile
    R = np.linspace(R1, Rt, N_layers)          #WILL CHANGE ONCE I HAVE MORE INFO ON NOZZLE GEOMTRY
    length= R1*(1-1/R1_Rt_ratio)/np.tan(np.radians(contr_angle))
    return R, length

def throattoexit_radial_profile(Rt, R2_Rt_ratio, exp_angle):
    R2 = Rt*R2_Rt_ratio
    length_exit= np.abs(Rt*(1-R2_Rt_ratio)/np.tan(np.radians(exp_angle)))
    R= np.linspace(Rt, R2, N_layers)        #WILL CHANGE ONCE I HAVE MORE INFO ON NOZZLE GEOMTRY
    return R, length_exit

def Volume(R, length):
    x_frustrums= length/N_layers    #length of each frustrum (truncated cone)          WILL CHANGE ONCE I HAVE MORE INFO ON NOZZLE GEOMTRY

    #calculate different volumes of each layer
    V= []
    for i in range(N_layers-1): #          WILL CHANGE ONCE I HAVE MORE INFO ON NOZZLE GEOMTRY
        V_x= 1/3*np.pi*(R[i]**2+R[i]*R[i+1]+R[i+1]**2)*x_frustrums #volume of frustrum
        V.append(V_x)
    return V

def Layer_volums(R1, R1_Rt_ratio, contr_angle, Rt, R2_Rt_ratio, exp_angle):
    global R_total
    global V_total
    R, length = chambertothroat_radial_profile(R1, R1_Rt_ratio, contr_angle)
    V = Volume(R, length)
    R_exit, length_exit = throattoexit_radial_profile(Rt, R2_Rt_ratio, exp_angle)
    V_exit = Volume(R_exit, length_exit)
    V_total= np.array(V + V_exit)
    R_total = np.array(R + R_exit)
    #make R_total global so it can be used in other functions
    return V_total

V_total = np.array(Layer_volums(R1, R1_Rt_ratio, contr_angle, Rt, R2_Rt_ratio, exp_angle))

#calculate temperature of each layer knowing temperature of first layer and using adiaabatic relations

def adiabatic_temperature(T_initial, V, gamma_mix):
    T= [T_initial]
    for i in range(N_layers-1):
        T_x= T[i]*(V[i]/V[i+1])**(gamma_mix-1)
        T.append(T_x)

    return T

T= adiabatic_temperature(T_initial, V_total, gamma_mix)

#use Bartz equation to calculate heat flow at each layer

def Bartz(Tg, Tw, h):
    q=h*(Tg-Tw)
    return q
def bartz_coeff(D, mu,cp,Pr,pc,c_char,rc,At,A,sigma):
    h=0.026/D**0.2*(mu**0.2)*(cp)/(Pr**0.6)*(pc**0.8)/(c_char**0.8)*(D**0.1)/(rc**0.1)*(At**0.9)/(A**0.9)*(sigma)
    return h

def sigma(Tw, Tc, gamma,M,w):
    sigma= 1/((0.5*Tw/Tc*(1+(gamma-1)/2*M**2)+0.5)**(0.8-w/5)*(1+(gamma-1)/2*M**2)**(w/5))
    return sigma


#constants:
D=Rt*2
mu=1.8e-5       #DON'T KNOW!!!
cp=1000         #DON'T KNOW!!!
Pr=0.7          #DON'T KNOW!!!
pc=101325       #DON'T KNOW!!!
c_char=1210     #DON'T KNOW!!!
rc=287         #DON'T KNOW!!!
At=np.pi*Rt**2
A=np.pi*R_total**2            #DON'T KNOW!!!
w=0.4           #DON'T KNOW!!!
M=R_total*1.5           #DON'T KNOW!!!

#what is Tw???
Tw= 2000 #       #DON'T KNOW!!!????
sigma= sigma(Tw, T, gamma_mix,M,w)
h= bartz_coeff(D, mu,cp,Pr,pc,c_char,rc,At,A,sigma)