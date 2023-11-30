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
Masses = np.array(["OEW", "MTOW", "MTOW-W_f"]) #kg

Heights = np.array(["standard sea level altitude", "cruise altitude"])

#Velocities
V_c = 257.6 #m/s
V_d0 = 1.25*V_c #m/s - sea level
V_dc = 301.8 #m/s due to sound barrier (so M_d = 1) - cruise altitude
V_d = np.array([V_d0, V_dc])
V_s1 = np.sqrt((2*g*masses[1])/(densities[0]*S*C_Lmax1))
V_s0 = np.sqrt((2*g*masses[1])/(densities[0]*S*C_Lmax0))
V_f = 1.8*V_s0

#loads
n_max = 2.1 + 24000/(kgTOlb(masses)+10000)
n_max[n_max > 3.8] = 3.8
n_max[n_max < 2.5] = 2.5
n_min = -1


V_EAS = np.linspace(0, 700, 1000)

for j in range(len(masses)):
    for i in range(len(heights)):
    
        V_C = msTOknt(EAS(V_c, heights[i]))
        V_D = msTOknt(EAS(V_d[i], heights[i]))
        
        V_S1 = msTOknt(EAS(np.sqrt((2*g*masses[j])/(densities[i]*S*C_Lmax1)), heights[i]))
        V_S0 = msTOknt(EAS(np.sqrt((2*g*masses[j])/(densities[i]*S*C_Lmax0)), heights[i]))
        V_F = msTOknt(EAS(1.8*V_s0, heights[i]))
        
        
        plt.plot(V_EAS[:np.sum(V_EAS<V_S1*n_max[j]**(1/2))], ((V_EAS[:np.sum(V_EAS<V_S1*n_max[j]**(1/2))]/V_S1)**2), color = "black")
        plt.plot(V_EAS[:np.sum(V_EAS<V_S0*2**(1/2))], (V_EAS[:np.sum(V_EAS<V_S0*2**(1/2))]/V_S0)**2, color = "black")
        plt.plot(V_EAS[:np.sum(V_EAS<V_S1)], -(V_EAS[:np.sum(V_EAS<V_S1)]/V_S1)**2, color = "black")
        plt.plot(V_EAS[np.sum(V_EAS<V_C):np.sum(V_EAS<V_D)], ((1)/(V_D-V_C))*(V_EAS[np.sum(V_EAS<V_C):np.sum(V_EAS<V_D)]-V_D), color = "black")
        
        print(f"V_C = {V_C}, V_S1 = {V_S1}, V_S0 = {V_S0}, V_F = {V_F}, V_D = {V_D}, h = {heights[i]}")
        
        plt.hlines(y=n_max[j], xmin = V_S1*n_max[j]**(1/2), xmax = V_D, colors = "black")
        plt.hlines(y=n_min, xmin = V_S1, xmax = V_C, colors = "black")
        plt.hlines(y=2, xmin = V_S0*2**(1/2), xmax = V_S1*2**(1/2), colors = "black")
        plt.vlines(x = V_D, ymin = 0, ymax = n_max[j], colors = "black")
        
        plt.vlines(x=V_S1, ymin = 0, ymax = 1, colors = "#404040", linestyle="--", lw = 0.8)
        plt.vlines(x=V_S1*n_max[j]**(1/2), ymin = 0, ymax = n_max[j], colors = "#404040", linestyle="--", lw = 0.8)
        plt.vlines(x=V_C, ymin = -1, ymax = n_max[j], colors = "#404040", linestyle="--", lw = 0.8)
        plt.hlines(y=1, xmin = 0, xmax = V_D+V_D/30, colors = "#404040", linestyles = "dashed", lw = 0.8)
        plt.hlines(y=-1, xmin = 0, xmax = V_F, colors = "#404040", linestyles = "dashed", lw = 0.8)
        plt.hlines(y=0, xmin = 0, xmax = V_D*1.3, colors = "Black", lw = 1)

        #plt.axhline(y=n_max[1], xmin = V_s1*n_max[1]**(1/2), xmax=V_d)
        
        
        #plt.text(V_S1+V_D/100, 0.1, f'V_S', fontsize=8, ha = 'left')
        #plt.text(V_S1*n_max[j]**(1/2)+V_D/100, 0.1, f"V_S*(n_max)^(1/2)", fontsize = 8, ha = "left")
        #plt.text(V_C+V_D/100, 0.1, f"V_C", fontsize = 8, ha = "left")
        #plt.text(V_D+V_D/100, 0.1, f"V_D", fontsize = 8, ha = "left")
        #plt.text(V_C*0.7, n_max[j]+0.12, f"n_max", fontsize = 8, ha = "center")
        #plt.text(V_C*0.6, -1.2, f"n_min", fontsize = 8, ha = "center")

        
        plt.title(f"V-n diagram for M = {Masses[j]} at {Heights[i]}")
        
        # Add horizontal axis line at y=0
        #plt.axhline(0, color='black',h=0.5)
        
        # Add vertical axis line at x=0
        plt.axvline(0, color='black',linewidth=0.5)
        
        plt.grid(True)
        plt.ylim(-1.5, n_max[1]+0.5)
        #plt.xlim(0, V_D*1.1)
        plt.xlim(0, 750)
        plt.xlabel("V_EAS [KTS]")
        plt.ylabel("n", rotation = 0)

        plt.show()
