# Assuming that the spars is an I beam. 
import matplotlib.pyplot as plt
import numpy as np
from USE_THIS_loads_main import get_integrated  

wing_length = 22.445
c_root = 10.10
c_tip = 2.73
E = 68.9 *10**9 #Pa
v = 0.33 #poisson's ratio
t = 0.02 # thickness of the web 
ks =  9.5 # clamped edges

def cord_function(y):# c(y) gives the cord lenght as function of span location y
    return c_tip+((c_root-c_tip)/wing_length)*y

def short_side_of_the_spar_function(y):
    return 0.0942*cord_function(y)

def calculate_critical_shear_stress(y):
    return ((np.pi)**(2)*ks*E/(12*(1-v**2)))*(t/short_side_of_the_spar_function(y))**2

b = np.linspace(0,wing_length,50)
plt.plot(b,calculate_critical_shear_stress(b),color = 'blue')
plt.plot(b,get_integrated()[1]/(2*t*short_side_of_the_spar_function(b)),color = 'red')
plt.xlabel("Span [m]")
plt.ylabel("Critical shear stress in the spar web [Pa]")
plt.show()

