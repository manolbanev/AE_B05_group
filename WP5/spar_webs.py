# Assuming that the spars is an I beam. (from the manual) In order to reduce the workload, you may assume that only the spar webs carry any shear flow due
# to the shear force, and you may use your engineering judgement to determine the maximum shear
# stress due to the shear force in the webs (rather than going through the hassle of determining the
# shear flow distribution for a wing box with stringers);
import matplotlib.pyplot as plt
import numpy as np
# from USE_THIS_loads_main import get_integrated  
import pandas as pd
import math

wing_length = 22.445
c_root = 10.10
c_tip = 2.73
E = 68.9 *10**9 #Pa
v = 0.33 #poisson's ratio
t = 0.02 # thickness of the web 
ks =  9.5 # clamped edges
Vy = 541111 #N
T = 38300  #Nm
Vx = 0 #N

def cord_function(y):# c(y) gives the cord lenght as function of span location y
    return c_root+((c_root-c_tip)/wing_length)*y

def short_side_of_the_spar_function(y):
    return 0.0942*cord_function(y)

def calculate_critical_shear_stress(y):
    return ((np.pi)**(2)*ks*E/(12*(1-v**2)))*(t/short_side_of_the_spar_function(y))**2

def position_of_booms():
    return [[0.25,0.044],[0.25,-0.044],[-0.25,-0.044],[-0.25,0.044]] #[x,y] in m

def booms_area():
    return [1,1,1,1] #m^2

def calculate_Ixx(position_of_booms,booms_area):
    Ixx = 0
    for i in range(4):
        Ixx = Ixx + booms_area[i]*position_of_booms[i][1]**2
    return Ixx
        

def calculate_Iyy(position_of_booms,booms_area):
    Iyy = 0
    for i in range(4):
        Iyy = Iyy + booms_area[i]*(position_of_booms[i][0]-0.35)**2
    return Iyy

def calculate_total_shear(Vy,T,h,t,c):
    total_shear = []
    for i in range(50):
        total_shear.append((0.75*Vy[i])/(h[i]*t) + (T[i])/(2*c[i]*h[i]*t))


# def plot_wing_box(position_of_booms):

#     x_booms = []
#     y_booms = []
#     plt.xlim(-3, 3)
#     plt.ylim(-3, 3)
#     plt.axis('equal')
#     X, Y = [], []
#     for line in open('NASA_sc_20412.dat', 'r'):
#         values = [float(s) for s in line.split()]
#         X.append(values[0])
#         Y.append(values[1])
#     for i in range(4):
#         x_booms.append(position_of_booms[i][0])
#         y_booms.append(position_of_booms[i][1])
#     plt.plot(X, Y,color='red')
#     plt.plot(x_booms[:2],y_booms[:2],marker='o', linestyle='-', color='blue')
#     plt.plot(x_booms[-2:],y_booms[-2:],marker='o', linestyle='-', color='blue')
#     plt.show()

def compute_delta_q_b(Vy,Vx,Ixx,Iyy,position,A):
    delta_q_b = []
    for i in range(4):
        delta_q_b.append(-(Vy/Ixx)*A[i]*position[i][1]-(Vx/Iyy)*A[i]*position[i][0])
    return delta_q_b

# def compute_q_s0():
#      return  -np.sum(compute_delta_q_b(Vy,Vx,calculate_Ixx(position_of_booms(),booms_area()),calculate_Iyy(position_of_booms(),booms_area()),position_of_booms(),booms_area()))

def compute_q_critical():
    return compute_delta_q_b(Vy,Vx,calculate_Ixx(position_of_booms(),booms_area()),calculate_Iyy(position_of_booms(),booms_area()),position_of_booms(),booms_area())[0]+compute_q_due_to_torque(T,position_of_booms())


def compute_q_due_to_torque(T,position_of_booms):
    enclosed_area = (-position_of_booms[0][0] + position_of_booms[2][0])*(abs(position_of_booms[0][1])+ abs(position_of_booms[1][1]))
    q_t  = -T/(2*enclosed_area)
    return q_t



def main():
    # print(calculate_Ixx(position_of_booms(),booms_area()))
    # print(calculate_Iyy(position_of_booms(),booms_area()))
    #plot_wing_box(position_of_booms())
    # print(compute_delta_q_b(Vy,Vx,calculate_Ixx(position_of_booms(),booms_area()),calculate_Iyy(position_of_booms(),booms_area()),position_of_booms(),booms_area()))
    # print('q_cr',compute_q_critical(),compute_q_due_to_torque(T,position_of_booms()))
    # print(calculate_t(0.044*c_root,0.088*0.5*c_root*c_root))
    return 1
main()

 

# b = np.linspace(0,wing_length,50)
# plt.plot(b,calculate_critical_shear_stress(b),color = 'blue')
# plt.plot(b,get_integrated()[1]/(2*t*short_side_of_the_spar_function(b)),color = 'red')
# plt.xlabel("Span [m]")
# plt.ylabel("Critical shear stress in the spar web [Pa]")
# plt.show()
