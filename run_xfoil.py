import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt


names_of_airfoils = ['NACA2415','NACA_22112','NASA_sc_20710','NACA_23012','BACXXX','boeing_106','boeing_737_midspan','NASA_sc_20610','NASA_sc_20410','NASA_sc_20412','NASA_sc_20414']
file_name = "polar_file.txt"
names_of_airfoils_2 = ['NASA_sc_20410','NASA_sc_20412','NASA_sc_20414']
names_of_airfoils_3 = ['NASA_sc_20710','NASA_sc_20610']
alpha_i = -5
alpha_f = 30
alpha_step = 0.1
Re = 29610844.47
Re_2 = 22300000
M_2 = 0.2
n_iter = 100
M = 0.67
cl_needed = 0.656
run_simulation = False
cd_alpha = False

def run_xfoil(airfoil_name):
    if os.path.exists("polar_file_{0}_take_off.txt".format(airfoil_name)):
        os.remove("polar_file_{0}_take_off.txt".format(airfoil_name))

    input_file = open("input_file.in", 'w')
    input_file.write("LOAD {0}.dat\n".format(airfoil_name))
    input_file.write(airfoil_name + '\n')
    input_file.write("PANE\n")
    input_file.write("OPER\n")
    input_file.write("Visc {0}\n".format(Re_2))
    input_file.write("M {0}\n".format(M_2))
    input_file.write("PACC\n")
    input_file.write("polar_file_{0}_take_off.txt\n\n".format(airfoil_name))
    input_file.write("ITER {0}\n".format(n_iter))
    input_file.write("ASeq {0} {1} {2}\n".format(alpha_i, alpha_f,
                                                alpha_step))
    input_file.write("\n\n")
    input_file.write("quit\n")
    input_file.close()

    subprocess.call("xfoil.exe < input_file.in",shell=True)
if run_simulation:
    for i in range(3):
        run_xfoil(names_of_airfoils_2[i])
else:
    if not cd_alpha:
        figure, axis = plt.subplots(1, 2, figsize=(10,7))
        figure.tight_layout(pad=2.0)
        for i in range(3):
            data = np.loadtxt("polar_file_{0}_take_off.txt".format(names_of_airfoils_2[i]), skiprows=12)
            # data = np.loadtxt('polar_9m_m0_1408.txt',skiprows=12)
            cd=[]
            cl=[]
            alpha = []
            for d in data:
                cd.append(d[2])
                cl.append(d[1])
                alpha.append(d[0])
            axis[0].plot(alpha,cl,label = names_of_airfoils_2[i]+" slope: "+  str(np.format_float_positional(np.polyfit(alpha,cl,1)[0],2)),linewidth=3)
            axis[0].plot(alpha,cl,label = 'R = 22.3M , M = 0.2, XFoil data'+ names_of_airfoils_2[i],linewidth=3)
            axis[1].plot(cd,cl,'o',label = 'R = 22.3M , M = 0.2, XFoil data'+ names_of_airfoils_2[i],linewidth=3)
        axis[0].set_xlabel(u'\u0251'+u'\u00B0',fontsize=16)
        axis[0].set_ylabel('Cl',fontsize = 16)
        axis[0].legend(fontsize = 14, loc="lower right")
        # axis[0].axhline(y = cl_needed,color = 'red')
        axis[0].grid()
        axis[1].set_xlabel('Cd',fontsize=16)
        axis[1].set_ylabel('Cl',fontsize = 16)
        axis[1].legend(fontsize = 14, loc="lower right")
        axis[1].grid()
        plt.show()

    # cd-alpha and cl/cd - alpha
    else:
        figure, axis = plt.subplots(1, 2, figsize=(10,7))
        figure.tight_layout(pad=2.0)
        for i in range(3):
            data = np.loadtxt("polar_file_{0}_take_off.txt".format(names_of_airfoils_2[i]), skiprows=12)
            cd=[]
            cl=[]
            cl_cd = []
            alpha = []
            for d in data:
                cd.append(d[2])
                cl.append(d[1])
                cl_cd.append(d[1]/d[2])
                alpha.append(d[0])
            axis[0].plot(alpha,cd,'o',label = names_of_airfoils_2[i],linewidth=3)
            axis[1].plot(alpha,cl_cd,'o',label = names_of_airfoils_2[i],linewidth=3)
        axis[0].set_xlabel(u'\u0251'+u'\u00B0',fontsize=16)
        axis[0].set_ylabel('Cd',fontsize = 16)
        axis[0].legend(fontsize = 14, loc="upper left")
        axis[0].grid()
        axis[1].set_xlabel(u'\u0251'+u'\u00B0',fontsize=16)
        axis[1].set_ylabel('Cl/Cd',fontsize = 16)
        axis[1].legend(fontsize = 14, loc="lower right")
        axis[1].grid()
        plt.show()



