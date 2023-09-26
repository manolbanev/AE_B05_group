"""Runs an XFOIL analysis for a given airfoil and flow conditions"""
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt


# airfoil_name = "NASA_sc_20710"
names_of_airfoils = ['NACA2415','NACA_22112','NASA_sc_20610','NASA_sc_20710','NACA_23012','BACXXX','NASA_sc_20410','boing_106','boing_737_midspan']
file_name = "polar_file.txt"
alpha_i = 0
alpha_f = 3
alpha_step = 0.1
Re = 29610844.47
n_iter = 100
M = 0.67

def run_xfoil(airfoil_name):
    if os.path.exists("polar_file_{0}.txt".format(airfoil_name)):
        os.remove("polar_file_{0}.txt".format(airfoil_name))

    input_file = open("input_file.in", 'w')
    input_file.write("LOAD {0}.dat\n".format(airfoil_name))
    input_file.write(airfoil_name + '\n')
    input_file.write("PANE\n")
    input_file.write("OPER\n")
    input_file.write("Visc {0}\n".format(Re))
    input_file.write("M {0}\n".format(M))
    input_file.write("PACC\n")
    input_file.write("polar_file_{0}.txt\n\n".format(airfoil_name))
    input_file.write("ITER {0}\n".format(n_iter))
    input_file.write("ASeq {0} {1} {2}\n".format(alpha_i, alpha_f,
                                                alpha_step))
    input_file.write("\n\n")
    input_file.write("quit\n")
    input_file.close()

    subprocess.call("xfoil.exe < input_file.in",shell=True)

    # print(polar_data)
run_xfoil(names_of_airfoils[-2])

# figure, axis = plt.subplots(3, 2, figsize=(10,7))
# figure.tight_layout(pad=2.0)
# m=0
# for i in range(3):
#     for n in range(2):
#         data = polar_data = np.loadtxt("polar_file_{0}.txt".format(names_of_airfoils[m]), skiprows=12)
#         cd=[]
#         cl=[]
#         alpha = []
#         for d in data:
#             cd.append(d[2])
#             cl.append(d[1])
#             alpha.append(d[0])
#         axis[0,0].plot(alpha, cl,'o')
#         axis[0,0].set_xlabel('Alpha')
#         axis[0,0].set_ylabel('Cl')
#         axis[0,0].set_title('Cl/Alpha')
#         m=m+1
# plt.show()

