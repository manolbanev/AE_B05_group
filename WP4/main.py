import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import interpolate 
import random

q = 61.25
c_y = [9.84,2.45]
y = [0,21.665]
C_y_func = sp.interpolate.interp1d(y,c_y,kind='linear',fill_value="extrapolate")# c(y) = 9.84 - 0.3411y

file_names = ['MainWing_a=0.00_v=10.00ms.txt','MainWing_a=10.00_v=10.00ms.txt']

def load_file(filename):
    data = [] # [y-span,chord,ai,cl,icd,cm]
    used_cols = [0,1,2,3,5,7]
    positive_indeces = []
    for col in range(len(used_cols)):
        data.append(np.genfromtxt(fname=filename,skip_header=21,skip_footer=1029,usecols=used_cols[col]))
    Cllst = []
    Cdlst = []
    Cmlst = []
    ylst = [x for x in data[0] if x >= 0]
    for m in range(len(ylst)):
        for i in range(len(data[0])):
            if ylst[m] == data[0][i]:
                positive_indeces.append(i)
    for n in range(len(ylst)):
        Cllst.append(data[3][positive_indeces[n]])
        Cdlst.append(data[4][positive_indeces[n]])
        Cmlst.append(data[5][positive_indeces[n]])
    Cl_func = sp.interpolate.interp1d(ylst,Cllst,kind='cubic',fill_value="extrapolate")
    Cd_func = sp.interpolate.interp1d(ylst,Cdlst,kind='cubic',fill_value="extrapolate")
    Cm_func = sp.interpolate.interp1d(ylst,Cmlst,kind='cubic',fill_value="extrapolate")

    L_per_unit_span_func = q*Cl_func(np.linspace(0,21.665,50))*C_y_func(np.linspace(0,21.665,50))
    D_per_unit_span_func = q*Cd_func(np.linspace(0,21.665,50))*C_y_func(np.linspace(0,21.665,50))
    M_per_unit_span_func = q*Cm_func(np.linspace(0,21.665,50))*C_y_func(np.linspace(0,21.665,50))**2

    # for t in range(len(ylst)):
    #     ylst[t] = ylst[t]+ random.uniform(0, 1) #testing how good is the interpolation
    # plt.plot(data[0],data[5],'o')
    # plt.plot(ylst,Cm_func(ylst),'red')
    # plt.show()
    plt.plot(np.linspace(0,21.665,50),L_per_unit_span_func,label = 'Lift per unit span [N/m]')
    plt.plot(np.linspace(0,21.665,50),D_per_unit_span_func,label= 'Drag per unit span [N/m]')
    plt.plot(np.linspace(0,21.665,50),M_per_unit_span_func,label= 'Moment per unit span [N]')
    plt.legend()
    plt.xlabel('Half wing span [m]')
    plt.ylabel('Aerodynamic loading per unit span [N/m]')
    plt.show()
    

load_file(file_names[1])
load_file(file_names[0])

