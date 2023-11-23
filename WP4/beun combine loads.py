import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy import integrate


span = np.linspace(0, 22.445, 100)


def get_aerodynamic(x):

    q = 61.25
    c_y = [9.84, 2.45]
    y = [0, 22.445]
    half_span = 22.445

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

    aerodynamic_output = L_prime_func(x, find_Cl_alpha(x,0))
    return aerodynamic_output


def get_inertial(point):

    c_root = 10.10  # m
    c_tip = 2.73  # m

    def chord_distribution(x, c_root, c_tip, wingspan):
        return c_root - ((c_root - c_tip) / wingspan) * x

    # plt.plot(span, chord_distribution(span, c_root, c_tip, 22.445))

    def A(x):
        chord_at_x = chord_distribution(x, c_root, c_tip, wingspan=22.445)
        height = 0.0942 * chord_at_x
        length = 0.5 * chord_at_x
        return 0.95 * height * length


    def Wfps(x):
        W_fperspan = (9.81 * A(x) * 800)[:800]
        #W_fperspan = np.concatenate((W_fperspan, np.zeros(200)))
        return W_fperspan


    def WW(x):
        WA = A(x) * 9.81 * 173.7434
        return (WA)

    inertial_loading = WW(point) + Wfps(point)

    return inertial_loading


def get_moment(distribution, limit):
    moment_distribution, moment_error = sp.integrate.quad(distribution, 0, limit)
    return moment_distribution


def get_integrated_idiot():
    values = []
    for point in span:
        values.append(get_moment(get_aerodynamic, point))
        print('value: ', values, 'point: ', point)

    return values


def engine_weight():
    values = []
    for i in span:
        if abs(7.5 - i) < 000.1:
            values.append(32480)
        else:
            values.append(0)
    return values


def combined_load(x):
    combined = 625 * get_aerodynamic(x) - get_inertial(x) - engine_weight()
    return combined


fig, ax = plt.subplots()
ax.set_xlabel('Half wing span [m]')
ax.set_ylabel('Load [N]')
x = span
line1 = ax.plot(x, combined_load(x))
#line2 = ax.plot(get_integrated_idiot())
plt.show()




