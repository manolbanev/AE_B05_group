import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import interpolate
import math
import matplotlib.animation as animation

q = 61.25
c_y = [9.84, 2.45]
y = [0, 21.665]
half_span = 21.665

alpha_step = 0

CL_0 = 0.148880
CL_10 = 0.901779
Cd_0 = 0.001000
Cd_10 = 0.036043
Cm_0 = -0.24373
Cm_10 = -1.17556

C_y_func = sp.interpolate.interp1d(y, c_y, kind='linear', fill_value="extrapolate")  # c(y) = 9.84 - 0.3411y

file_names = ['MainWing_a=0.00_v=10.00ms.txt', 'MainWing_a=10.00_v=10.00ms.txt']


def load_file(filename):
    data = []  # [y-span,chord,ai,cl,icd,cm]
    used_cols = [0, 1, 2, 3, 5, 7]
    positive_indeces = []
    for col in range(len(used_cols)):
        data.append(np.genfromtxt(fname=filename, skip_header=21, skip_footer=1029, usecols=used_cols[col]))
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
    return Cllst, Cdlst, Cmlst, ylst


def interpolate(Cllst, Cdlst, Cmlst, ylst):
    Cl_func = sp.interpolate.interp1d(ylst, Cllst, kind='cubic', fill_value='extrapolate')
    Cd_func = sp.interpolate.interp1d(ylst, Cdlst, kind='cubic', fill_value='extrapolate')
    Cm_func = sp.interpolate.interp1d(ylst, Cmlst, kind='cubic', fill_value='extrapolate')
    return Cl_func, Cd_func, Cm_func


def L_prime_func(y, Cl_func):
    return Cl_func * q * C_y_func(y)


def D_prime_func(y, Cd_func):
    return Cd_func * q * C_y_func(y)


def M_prime_func(y, Cm_func):
    return Cm_func * q * C_y_func(y) ** 2


def find_cl_d(alpha):
    coef = math.sin(math.radians(alpha)) / math.sin(math.radians(10))
    Cl_d = coef * (CL_10 - CL_0) + CL_0
    return Cl_d


def find_Cl_alpha(y, alpha):
    Cl_d = find_cl_d(alpha)
    Cl_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                            load_file(file_names[0])[3])[0]
    Cl_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                             load_file(file_names[1])[3])[0]
    return Cl_func_0(y) + ((Cl_d - CL_0) / (CL_10 - CL_0)) * (Cl_func_10(y) - Cl_func_0(y))


def find_cd_d(alpha):
    return Cd_0 + ((Cd_10 - Cd_0) / (CL_10 ** 2 - CL_0 ** 2) * (find_cl_d(alpha) ** 2 - CL_0 ** 2))


def find_Cd_alpha(y, alpha):
    Cd_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                            load_file(file_names[0])[3])[1]
    Cd_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                             load_file(file_names[1])[3])[1]
    return Cd_func_0(y) + ((Cd_func_10(y) - Cd_func_0(y)) / (Cd_10 - Cd_0)) * (find_cd_d(alpha) - Cd_0)


def find_cm_d(alpha):
    return math.sin(math.radians(alpha)) / math.sin(math.radians(10)) * (Cm_10 - Cm_0) + Cm_0


def find_Cm_alpha(y, alpha):
    Cm_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                            load_file(file_names[0])[3])[2]
    Cm_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                             load_file(file_names[1])[3])[2]
    return Cm_func_0(y) + ((find_cm_d(alpha) - Cm_0) / (Cm_10 - Cm_0)) * (Cm_func_10(y) - Cm_func_0(y))


# animation
alpha = 0


def animate(i):
    line1.set_ydata(L_prime_func(np.linspace(0, 21.665, 50), find_Cl_alpha(np.linspace(0, 21.665, 50), alpha + i / 10)))
    line2.set_ydata(D_prime_func(np.linspace(0, 21.665, 50), find_Cd_alpha(np.linspace(0, 21.665, 50), alpha + i / 10)))
    line3.set_ydata(M_prime_func(np.linspace(0, 21.665, 50), find_Cm_alpha(np.linspace(0, 21.665, 50), alpha + i / 10)))
    ax.set_title('Alpha ' + str(round(alpha + i / 10, 2)))
    return line1, line2, line3


fig, ax = plt.subplots()
ax.set_xlabel('Half wing span [m]')
ax.set_ylabel('Lift per unit span [N/m]')
ax.set_ylim(-600, 600)
ax.set_xlim(0, 22)
x = np.linspace(0, 21.665, 50)
line1, = ax.plot(x, L_prime_func(x, find_Cl_alpha(x, alpha)))
line2, = ax.plot(x, D_prime_func(x, find_Cd_alpha(x, alpha)))
line3, = ax.plot(x, M_prime_func(x, find_Cm_alpha(x, alpha)))
ani = animation.FuncAnimation(
    fig, animate, interval=0.200, blit=False, save_count=50, frames=120)
plt.show()

# plot for particular alpha
# def main(alpha):
# plt.plot(np.linspace(0,21.665,50),L_prime_func(np.linspace(0,21.665,50),find_Cl_alpha(np.linspace(0,21.665,50),alpha)),label = 'Lift per unit span [N/m]')
# plt.plot(np.linspace(0,21.665,50),D_prime_func(np.linspace(0,21.665,50),find_Cd_alpha(np.linspace(0,21.665,50),alpha)),label= 'Drag per unit span [N/m]')
# plt.plot(np.linspace(0,21.665,50),M_prime_func(np.linspace(0,21.665,50),find_Cm_alpha(np.linspace(0,21.665,50),alpha)),label= 'Moment per unit span [N]')
# plt.legend()
# plt.xlabel('Half wing span [m]')
# plt.ylabel('Aerodynamic loading per unit span [N/m]')
# plt.title('Aerodynamic loading for '+ str(alpha) +' deg angle of attack')
#plt.show()


print(L_prime_func(np.linspace(0, 21.665, 50), find_Cl_alpha(np.linspace(0, 21.665, 50),0)))
