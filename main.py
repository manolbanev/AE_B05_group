import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

#const
Cd_0 = 0.05
A = 7#9
e = 0.447 #0.9
V = 256
ro_cruise = 0.441
w_s=0
ro=1.225
Vs=75
CLs_land= [2.6,2.8,3]
CLs_takeoff= [2.0,2.2,2.4]
CLs_clean= [1.6,1.8,2]
color_takeoff = ['salmon','red','darkred']
color_clean = ['darkblue','blue','deepskyblue']
color_land = ['bisque','gold','darkorange']
f = 0.65
k = 240
S_land = 1981.2
climb_grad = 0.04
convert_value = 47.88

# plot database points
file = pd.read_excel('AC database WP1.xlsx',sheet_name='Sheet2',header=None)
df = pd.DataFrame(file)
y_axis = file[0]
x_axis = file[1]
plt.figure(figsize = (10, 5))

plt.ylabel("T/W")
plt.xlabel("W/S")
plt.ylim(0,0.8)  
plt.grid(True)
plt.xlim(0,10000)  

#sizing for stall speed
def size_for_stall(Cl_max):
    w_s = ro*Vs*Vs*Cl_max/2
    return w_s
n= 0 
for cl in CLs_clean:
    plt.axvline(x = size_for_stall(cl), color = '{}'.format(color_clean[n]), label = 'CLmax = {} Stall'.format(cl))
    n = n+1
#sizing for landing
def sizing_for_landing(Cl_max):
        w_s2 = (Cl_max*ro*S_land/(0.5847))/(2*f)
        return w_s2
i= 0 
for cl in CLs_land:
    plt.axvline(x = sizing_for_landing(cl), color = '{}'.format(color_land[i]), label = 'CLmax = {} Land'.format(cl))
    i = i+1


#used for the filling the space between the functions
x1 = size_for_stall(1.8)

#sizing for criuse speed 
x=[]
y=[]
y1=[]
x = np.arange(1,10000)
def sizing_for_cruise(x):
    for j in range(9999):
        y.append(x[j]/(math.pi*A*e*0.5*ro_cruise*V**2)+(Cd_0*0.5*ro_cruise*V**2)/x[j])
    return y 
#used for the filling the space between the functions
for m in range(6196):
    y1.append(x[m]/(math.pi*A*e*0.5*ro_cruise*V**2)+(Cd_0*0.5*ro_cruise*V**2)/x[m])

plt.plot(x,sizing_for_cruise(x),color = 'purple',label = 'Vcr=0.85 Mach')

# sizing for take off 
def sizing_for_take_off(Cl_max):
    Cl_to = Cl_max/(1.1**2)
    return Cl_to
Take_off_x = [0,10000]
Take_off_y = []
p=0
for cl in CLs_takeoff:
    Take_off_y = []
    for i in range(0,2):
        Take_off_y.append(Take_off_x[i]/(k*sizing_for_take_off(cl)*convert_value))
    plt.plot(Take_off_x,Take_off_y,color = '{}'.format(color_takeoff[p]),label = 'CLmax = {} Take off'.format(cl))
    p= p+1


# size for climb gradient

y = climb_grad+2*(Cd_0/(math.pi*A*e))**(0.5)
x = [0,10000]
arr_y = [y,y]
# draw design point
plt.plot([6196],[0.30], 'ro')

plt.plot(x,arr_y,color = 'green',label= 'Climb gradient')
x = np.linspace(0, x1, 6196)
plt.fill_between(x,y1,1,color='lime')
plt.fill_between([0,10000],[0,0.4786271233639654],color = 'white')
#show database points
plt.scatter(x_axis, y_axis,color = 'blue')
plt.legend(loc = 'upper right')
plt.show()