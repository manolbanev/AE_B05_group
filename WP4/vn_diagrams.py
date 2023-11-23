import numpy as np
import matplotlib.pyplot as plt

g = 9.81 #m/s
R = 287 #J/kg
C_Lmax1 = 2.28 #flaps retracted
C_Lmax0 = 2.99 #Landing
S = 287.9 #m^2

def ISA(h):
    #For h<11000
    T = 288.15-0.0065*h
    p = 101325 * (T/288.15)**(-g/(R*-0.0065))
    rho = 1.225 * (T/288.15)**((-g/(R*-0.0065))-1)
    return T, p, rho

def EAS(V_TAS, h):
    V_EAS = V_TAS * (ISA(h)[2]/ISA(0)[2])**(1/2)
    return V_EAS

def mach(V, h):
    M = V / (1.4*R*ISA(h)[0])**(1/2)
    return M

def msTOknt(x_ms):
    x_knt = 1.9438*x_ms
    return x_knt

def kgTOlb(x_kg):
    x_lb = 0.454*x_kg
    return x_lb

#Altitudes, cruise altitude = 31,000 ft = 9448.8 m
h_0 = 0#m
h_c = 9448.8 #m
heights = np.array([h_0, h_c])

#Densities
densities = ISA(heights)[2]

#Weights
MTOW = 140000 #kg
OEW = 62000 #kg
W_f = 70000 #kg
masses = np.array([OEW, MTOW, MTOW-W_f]) #kg

#Velocities
#velocities = np.array([V_s0, V_s1, V_a, V_c, V_d, V_f])
V_c = 257.6 #m/s
V_d0 = 1.25*V_c #m/s - sea level
V_dc = 301.8 #m/s due to sound barrier (so M_d = 1) - cruise altitude
V_s1 = np.sqrt((2*g*masses[1])/(densities[0]*S*C_Lmax1))
V_s0 = np.sqrt((2*g*masses[1])/(densities[0]*S*C_Lmax0))

#loads
n_max = 2.1 + 24000/(kgTOlb(masses)+10000)
n_max[n_max > 3.8] = 3.8
n_max[n_max < 2.5] = 2.5
n_min = -1

V_EAS = np.linspace(0, 700, 1000)

plt.plot(V_EAS, (V_EAS/V_s1)**2, color = "black")
plt.plot(V_EAS, (V_EAS/V_s0)**2, color = "black")

plt.hlines(y=n_max[1], xmin = 0, xmax = 350, colors = "black")
plt.axhline(y=n_min, color = "black")
plt.axhline(y=2, color = "black")
plt.axvline(x = V_d0, color = "black")
#plt.axhline(y=n_max[1], xmin = V_s1*n_max[1]**(1/2), xmax=V_d)












# Add horizontal axis line at y=0
plt.axhline(0, color='black',linewidth=0.5)

# Add vertical axis line at x=0
plt.axvline(0, color='black',linewidth=0.5)

plt.ylim(-2, 4)
plt.xlim(0, 350)
plt.show()