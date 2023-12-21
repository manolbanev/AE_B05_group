import matplotlib.pyplot as plt
import numpy as np
import math 
from USE_THIS_loads_main import get_integrated  

wing_length = 22.445
c_root = 10.10
c_tip = 2.73
E = 68.9 *10**9 #Pa
v = 0.33 #poisson's ratio
t = 0.009 # thickness of the web 
ks =  9.5 # clamped edges
Vy = 541111 #N
T = 38300  #Nm
Vx = 0 #N
margin = 1.5


def cord_function(y):# c(y) gives the cord lenght as function of span location y
    return c_root - ((c_root-c_tip)/wing_length)*y

def short_side_of_the_spar_function(y):
    return 0.0942*cord_function(y)

def compute_tau_ave(V,t,h):
    return V/(2*h*t)

def compute_tau_torque(T,t,c,h):
    return T/(2*t*c*h)
def compute_tau_total(tau_ave,tau_torque):
    k_v = 1.5
    return tau_ave*k_v + tau_torque

def compute_tau_critical(h):
    return -((np.pi)**(2)*ks*E/(12*(1-v**2)))*(t/h)**2


def main():
    b = np.linspace(0,wing_length,50)
    V = get_integrated()[0]
    T = get_integrated()[1]
    tau_total_arr = []
    tau_critical = []
    for i in range(50):
        c = cord_function(b[i])
        h = short_side_of_the_spar_function(b[i])
        tau_ave = compute_tau_ave(V[i],t,h)
        tau_torque = compute_tau_torque(T[i],t,c,h)
        tau_tot = compute_tau_total(tau_ave,tau_torque)
        tau_total_arr.append(margin*tau_tot/10**6)
        tau_critical.append(compute_tau_critical(h)/10**6)
    plt.plot(b,tau_total_arr,label= 'total applied shear stress')
    plt.plot(b,tau_critical,label = 'Critical shear stress')
    plt.xlabel('Span position [m]')
    plt.ylabel('$\u03C4_{total}$ [MPa]')
    plt.title('Total applied and critical shear stress for tickness = '+ str(t*1000)+ 'mm')
    plt.legend()
    plt.show()
    
main()