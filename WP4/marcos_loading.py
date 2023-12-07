import numpy as np
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt

c_root = 10.10  # m
c_tip = 2.73    # m
x = np.linspace(0, 22.445, 1000)

def chord_distribution(x, c_root, c_tip, wingspan):
    return c_root - ((c_root - c_tip) / wingspan) * x

#plt.plot(span, chord_distribution(span, c_root, c_tip, 22.445))

def A(x):
    chord_at_x = chord_distribution(x, c_root, c_tip, wingspan=22.445)
    height = 0.0942 * chord_at_x
    length = 0.5 * chord_at_x
    return 0.95*height * length


g = sp.interpolate.interp1d(x,A(x),kind='cubic',fill_value="extrapolate")

def Wfps(x):
    W_fperspan = (9.81*A(x)*800) [:800]
    W_fperspan = np.concatenate((W_fperspan, np.zeros(200)))
    return W_fperspan

def WW(x):
    WA = A(x)*9.81*173.7434
    return(WA)

g = sp.interpolate.interp1d(x,WW(x),kind='cubic',fill_value="extrapolate")
integrate.quad(g, 0, 22.445)

plt.plot(x, WW(x), label = "Wing structural weight")

f = sp.interpolate.interp1d(x,Wfps(x),kind='cubic',fill_value="extrapolate")

plt.plot(x, f(x), label = "Fuel weight")

print("The fuel weight in N for one wing is ", integrate.quad(f, 0, 22.445)[0])
print("The fuel weight in N for one wing is ", integrate.quad(g, 0, 22.445)[0])

plt.xlabel('Span (m)')
plt.ylabel('Weight per unit span N/m')
plt.title('Weight Distribution Along the Span')
plt.legend()
plt.show()
