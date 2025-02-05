import numpy as np
import matplotlib.pyplot as plt

# Constants
N_layers = 1000 #number of layers

gamma_ox= 7/5 #adiabatic index for diatomic gas
gamma_ipa = 4/3 #adiabatic index for non linear poliatomic gas

ox_ipa_ration=1.86 #ratio between oxygen and IPA

gamma_mix = (ox_ipa_ration*gamma_ox+gamma_ipa)/(ox_ipa_ration+1) #adiabatic index for the mixture

#initial conditions/parameters
T_initial = 2000 #K         DONT KNOW!!

Rc=63.23/2 #mm 
Rt = 31.62/2 #mm
Re=57.07/2 #mm
Lcyl= 73.53 #mm
Lc= 120 #mm
Le= 60.68 #mm

R2=47.14 #mm
R1= 23.71 #mm

contr_angle=np.radians(30)  #angle betweeen chamber and throat
exit_angle=np.radians(8) #angle at exit of parabolic nozzle      


def cylindrical_radial_profile(R=Rc):
    return R

def beginning_chamber_radial_profile(x, R=R2, xc=Lcyl, yc=Rc-R2):
    return np.sqrt(R**2-(x-xc)**2)+yc

def middle_chamber_radial_profile(x, La, Rla, contr_angle=contr_angle, Lcyl=Lcyl):
    m=-np.tan(contr_angle)
    q= Rla-m*(La + Lcyl -t)
    
    return m*x+q

def end_chamber_radial_profile(x, R=R1, xc=Lc, yc=Rt+R1):
    return -np.sqrt(R**2-(x-xc)**2)+yc

def nozzle_radial_profile(x, Re=Re, exit_angle=exit_angle, Rt=Rt, Lc=Lc, Le=Le):
    alpha= np.tan(exit_angle)
    Lt=Lc+Le
    e=Lc/Rt
    f=Lt/Re
    g=2*Re*alpha

    A=np.array([[Rt, 1, 1/Rt], [Re, 1, 1/Re], [g, alpha, 0]])
    b=np.array([e, f, 1])
    a,b,c=np.linalg.solve(A,b)

#check!!! seems correct
    ba=b/a
    y=(-ba+np.sqrt(ba**2-4*(c-x)/a))/2

    return y

def Volume(R):
    x_frustrums= t    #length of each frustrum (truncated cone)

    #calculate different volumes of each layer
    V= []
    for i in range(N_layers): #          WILL CHANGE ONCE I HAVE MORE INFO ON NOZZLE GEOMTRY
        V_x= 1/3*np.pi*(R[i]**2+R[i]*R[i+1]+R[i+1]**2)*x_frustrums #volume of frustrum
        V.append(V_x)
    return V

def Layer_volums(N_layers=N_layers):
    global R
    global V
    global t
    global x
    R=np.ones(N_layers+1) #initialize R
    t= (Lc+Le)/N_layers #thickness of each layer
    x= t*np.arange(N_layers+1) #distance from chamber to each layer

    La=R2*np.sin(contr_angle) #length of beginning of chamber
    Lb=108 #mm #length before beginning of end of chamber


    #calculate radial profile of each layer
    Rla=100000
    for i in range(len(x)):
        if x[i]<Lcyl:
            R[i]= cylindrical_radial_profile()
        elif x[i]<(Lcyl + La) and x[i]>Lcyl:
            R[i]= beginning_chamber_radial_profile(x[i])
            if Rla> R[i]:
                Rla= R[i]
        elif x[i]<Lb and x[i]>(Lcyl + La):
            R[i]= middle_chamber_radial_profile(x[i], La, Rla)
        elif x[i]<Lc and x[i]>(Lb):
            R[i]= end_chamber_radial_profile(x[i])
        else:
            R[i]= nozzle_radial_profile(x[i])

    V = Volume(R)
    #make R_total global so it can be used in other functions
    return V

V = np.array(Layer_volums())

#calculate temperature of each layer knowing temperature of first layer and using adiaabatic relations

def adiabatic_temperature(T_initial, V, gamma_mix):
    T= [T_initial]
    for i in range(N_layers-1):
        T_x= T[i]*(V[i]/V[i+1])**(gamma_mix-1)
        T.append(T_x)

    return T

T= adiabatic_temperature(T_initial, V, gamma_mix)
T= np.array(T)

# Create mesh grid for plotting
X, Y = np.meshgrid(x, np.linspace(-np.max(R), np.max(R), 100))

# Initialize temperature grid with NaNs
T_2D = np.full_like(X, np.nan, dtype=float)  # Fill with NaN initially

# Assign temperatures inside the nozzle
for i in range(len(R) - 1):
    mask = (np.abs(Y) <= R[i])  # Keep values inside the nozzle
    
    T_2D[:, i][mask[:, i]] = T[i]  # Assign corresponding temperature

# Set color limits to match only valid temperature values
vmin, vmax = np.min(T), np.max(T)

# Plot the temperature distribution
c = plt.pcolormesh(X, Y, T_2D, cmap='hot', shading='auto', vmin=vmin, vmax=vmax)  # Fix color scaling
plt.plot(x, R, 'k', linewidth=2)  # Upper nozzle boundary
plt.plot(x, -R, 'k', linewidth=2)  # Lower nozzle boundary
plt.xlabel('Axial Position')
plt.ylabel('Radial Position')
plt.title('Nozzle Chamber Temperature Distribution')
plt.colorbar(c, label='Temperature (K)')
plt.show()




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
D=Rt*2  #throat diameter
mu=1.8e-5       #DON'T KNOW!!! # dynamic viscosity
cp=1000         #DON'T KNOW!!!  #specific heat
Pr=0.7          #DON'T KNOW!!!  #Prandtl number
pc=30e5       #chamber pressure
c_char=1210     #DON'T KNOW!!!  #characteristic velocity
rc=287         #DON'T KNOW!!!    #radius of curvature of throat (?)
At=np.pi*Rt**2  #throat cross sectional area
A=np.pi*R**2            #cross sectional area at each layer
w=0.4           #DON'T KNOW!!!  #exponent of viscosity-temperature power law for the combustion gas
M=R*1.5           #DON'T KNOW!!!  #Mach number at each layer

#what is Tw???
Tw= 2000 #       #DON'T KNOW!!!???? #wall temperature at each layer
sigma= sigma(Tw, T, gamma_mix,M,w)
h= bartz_coeff(D, mu,cp,Pr,pc,c_char,rc,At,A,sigma)