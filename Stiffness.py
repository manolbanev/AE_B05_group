import scipy as sp
import matplotlib.pyplot as plt
import math

#given
spar_height = 0.0942 #m
width = 0.5 #m 
c_root = 10.10 #m
c_tip = 2.73 #m
span = 44.89 #m
E = 68.9 * 10**9 #m
G = 26 * 10**9 #m
rho = 2700
nmax = 2.729
nmin = -1
nsafety = 1.5
#Variables
N_stringers_top = [4]
N_stringers_bottom = [4]
Stringer_change = [0]
S_stringers = 0.001 #m^2
t1 = 0.001 #m
t2 = 0.002

#chord and distance
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


#Stiffness
def I_xx(y) :  
    w = width * chord(y)
    h = spar_height * chord(y)
    I = S_stringers * (h/2)**2 * (nstringertop(y) + nstringerbot(y)) + 1/12 * t2 * h**3 + 2 * w * t1 * (h/2)**2
    return(I)

def J_z(y) :
    w = width * chord(y)
    h = spar_height * chord(y)
    J = S_stringers * ((chord(y))**2 * (d(nstringertop(y))**2 + d(nstringerbot(y))**2)) + 2 * w * t1 * (w**2 + t1**2)/12 + 2 * h * t2 * (h**2 + t2**2)/12 + 2 * w * t1 * (h/2)**2 + 2 * h * t2 *(w/2)**2
    return(J)


#Intigrals 
def v(y) :
    second = - 1/(I_xx(y) * E) #change 1 to bending
    return(second)
def intigral(y) :
    intigral, error0 = sp.integrate.quad(v, 0, y)
    return(intigral)

def theta(y) :
    tha = 1/(G * J_z(y)) #change 1 to torsion
    return(tha)

"""
#plotting graph
x = [0]
I = [I_xx(0)]
J = [J_z(0)]
j = 0
Deflection = [0]
Angle = [0]
while x[-1] < span/2 :
    j = j + 0.01
    x.append(j)
    I.append(I_xx(x[-1]))
    J.append(J_z(x[-1]))
    estimate1, error1 = sp.integrate.quad(intigral, 0, j)
    estimate2, error2 = sp.integrate.quad(theta, 0, j)
    Deflection.append(estimate1)
    Angle.append(estimate2 * 180/math.pi)



plt.subplot(221)
plt.plot(x, I)
plt.title("Moment of inertia diagram")
plt.ylabel("Moment of inertia [m^4]")

plt.subplot(223)
plt.plot(x, J)
plt.title("Torsional stiffness")
plt.xlabel("Spanwise position [m]")
plt.ylabel("Polar moment of inertia [m^4]")

plt.subplot(222)
plt.plot(x, Deflection)
plt.title("Deflection diagram")
plt.ylabel("Deflection[m]")

plt.subplot(224)    
plt.plot(x, Angle)
plt.title("Angle Diagram")
plt.xlabel("Spanwise position [m]")
plt.ylabel("Angle [deg]")

plt.show()
"""
#wingbox validity and weight
V_total = 2 * t1 * (chord(0) + chord(span/2))/2 * width * span/2 + 2 * t2 * (chord(0) + chord(span/2))/2 * spar_height * span/2 
l = 0
for i in range(int((len(Stringer_change)-1))) :
    l = l + 1
    V_total = V_total + (nstringerbot(i) + nstringertop(i)) * (Stringer_change[i+1] - Stringer_change[i]) * S_stringers
V_total = V_total + (nstringerbot(l) + nstringerbot(l)) * (span/2 - Stringer_change[-1] ) * S_stringers

estimate1, error1 = sp.integrate.quad(intigral, 0, span/2)
estimate2, error2 = sp.integrate.quad(theta, 0, span/2)
print("the deflection is between", nsafety * nmax * estimate1, "and", nsafety * nmin * estimate1, "m")
print("the twist angle is between", nsafety * nmax * estimate2 * 180/math.pi, "and", nsafety * nmin * estimate2 * 180/math.pi, "°")
print("The weight is ", 2 * V_total * rho, "kg")