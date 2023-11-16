import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt

N_stringers_top = [4]
N_stringers_bottom = [4]
Stringer_change = [0]
S_stringers = 0.5 #m^2

spar_height = 0.0942 #m
width = 0.5 #m 
t = 0.001 #m
c_root = 10.10
c_tip = 2.73
span = 44.89

#chord
def chord(a) :
    c = c_root - (c_root - c_tip) * 2 * a / span
    return(c)

def d(n) :
    distance = 0
    if n == 1 :
        return(spar_height/2)
    elif n % 2 == 0 :
        for k in range(int(n/2)) :
            distance = distance + ((spar_height/2)**2 + (width/2 - width/n * k )**2)**0.5
        distance = 2 * distance
    else :
        for k in range(int((n-1)/2)) :
            distance = distance + ((spar_height/2)**2 + (width/2 - width/n * k )**2)**0.5
        distance = distance * 2 + spar_height/2
    return(distance)

#Number of stringers :
nstringertop = sp.interpolate.interp1d(Stringer_change,N_stringers_top,kind="previous",fill_value="extrapolate")
nstringerbot = sp.interpolate.interp1d(Stringer_change,N_stringers_bottom,kind="previous",fill_value="extrapolate")

#Ixx
def I_xx(y) :  
    I = S_stringers * (width * spar_height * chord(y))**2 * (nstringertop(y) + nstringerbot(y)) + 1/3 * width * chord(y) * (spar_height * chord(y))**2 * t
    return(I)

def J_z(y) :
    J = 1/3 * spar_height * chord(y) * width * chord(y) * t * (spar_height * chord(y) + width * chord(y)) + S_stringers *chord(y) * (d(nstringertop(y)) + d(nstringerbot(y)))
    return(J)


x = [0]
I = [I_xx(0)]
J = [J_z(0)]
j = 0
while x[-1] < span/2 :
    j = j + 0.01
    x.append(j)
    I.append(I_xx(x[-1]))
    J.append(J_z(x[-1]))

plt.subplot(121)
plt.plot(x, I)
plt.title("Moment of inertia diagram")
plt.xlabel("Spanwise position")
plt.ylabel("Moment of inertia")

plt.subplot(122)
plt.plot(x, J)
plt.title("Torsional stiffness")
plt.xlabel("Spanwise position")
plt.ylabel("Polar moment of inertia")
print(d(2))
plt.show()

