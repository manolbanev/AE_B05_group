import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy import integrate

# constants
n_step = 50
wing_length = 22.445
span = np.linspace(0, wing_length, n_step)
pylon = span[int(0.32 * n_step):int(0.36 * n_step)]
engine_load = np.zeros_like(span)
engine_torques = np.zeros_like(span)
q_scale = 625
G = 26 * 10 ** 9
q = 61.25
cL_0 = 0.149649
cL_10 = 0.907325
Cm_0 = -0.24373
Cm_10 = -1.17556
Cd_0 = 0.001000
Cd_10 = 0.036043
alpha = 0
c_root = 10.10
c_tip = 2.73
spar_height = 0.0942
width = 0.5
E = 68.9 * 10 ** 9
rho = 2700
nmax = 2.729
nmin = -1
nsafety = 1.5
S_stringers = 0.0016
t1 = 0.001
t2 = 0.001
stringertop = 4
stringerbot = 4


# obtain aerodynamic loads distributions
def get_aerodynamic(x):
    c_y_func = sp.interpolate.interp1d([0, wing_length], [c_root, c_tip], kind='linear', fill_value="extrapolate")
    file_names = ['MainWing_a=0.00_v=10.00ms.txt', 'MainWing_a=10.00_v=10.00ms.txt']
    def load_file(filename): # read aerodynamic data from file
        data = []
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

    def interpolate(Cllst, Cdlst, Cmlst, ylst): # interpolate aerodynamic data to obtain a function
        Cl_func = sp.interpolate.interp1d(ylst, Cllst, kind='cubic', fill_value='extrapolate')
        Cd_func = sp.interpolate.interp1d(ylst, Cdlst, kind='cubic', fill_value='extrapolate')
        Cm_func = sp.interpolate.interp1d(ylst, Cmlst, kind='cubic', fill_value='extrapolate')
        return Cl_func, Cd_func, Cm_func

    # aerodynamic load distributions
    def L_prime_func(y, Cl_func):
        return Cl_func * q * c_y_func(y)

    def M_prime_func(y, Cm_func):
        return Cm_func * q * c_y_func(y) ** 2

    def D_prime_func(y, Cd_func):
        return Cd_func * q * c_y_func(y)

    def find_cl_d(alpha): # lift coefficient slope
        coef = math.sin(math.radians(alpha)) / math.sin(math.radians(10))
        Cl_d = coef * (cL_10 - cL_0) + cL_0
        return Cl_d

    def find_Cl_alpha(y, alpha): # lift coefficient vs angle of attack function
        Cl_d = find_cl_d(alpha)
        Cl_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                                load_file(file_names[0])[3])[0]
        Cl_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                                 load_file(file_names[1])[3])[0]
        return Cl_func_0(y) + ((Cl_d - cL_0) / (cL_10 - cL_0)) * (Cl_func_10(y) - Cl_func_0(y))

    def find_cm_d(alpha): # moment coefficient slope
        return math.sin(math.radians(alpha)) / math.sin(math.radians(10)) * (Cm_10 - Cm_0) + Cm_0

    def find_Cm_alpha(y, alpha): # moment coefficient vs angle of attack function
        Cm_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                                load_file(file_names[0])[3])[2]
        Cm_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                                 load_file(file_names[1])[3])[2]
        return Cm_func_0(y) + ((find_cm_d(alpha) - Cm_0) / (Cm_10 - Cm_0)) * (Cm_func_10(y) - Cm_func_0(y))

    def find_cd_d(alpha): # drag coefficient slope
        return Cd_0 + ((Cd_10 - Cd_0) / (cL_10 ** 2 - cL_0 ** 2) * (find_cl_d(alpha) ** 2 - cL_0 ** 2))

    def find_Cd_alpha(y, alpha): # drag coefficient vs angle of attack function
        Cd_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                                load_file(file_names[0])[3])[1]
        Cd_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                                 load_file(file_names[1])[3])[1]
        return Cd_func_0(y) + ((Cd_func_10(y) - Cd_func_0(y)) / (Cd_10 - Cd_0)) * (find_cd_d(alpha) - Cd_0)

    aerodynamic_output = L_prime_func(x, find_Cl_alpha(x, alpha))  # edit angle of attack
    torque_output = M_prime_func(x, find_Cm_alpha(x, alpha))
    drag_output = D_prime_func(x, find_Cd_alpha(x, alpha))

    return aerodynamic_output, torque_output, drag_output


# obtain inertial loads distributions
def get_inertial(point):
    def chord_distribution(x, c_root, c_tip, wingspan):
        return c_root - ((c_root - c_tip) / wingspan) * x

    def A(x): # obtain area, height and length of the wingbox
        chord_at_x = chord_distribution(x, c_root, c_tip, wingspan=22.445)
        height = 0.0942 * chord_at_x
        length = 0.5 * chord_at_x
        return 0.95 * height * length, height, length

    def Wfps(x): # obtain fuel per span
        W_fperspan = (9.81 * A(x)[0] * 800)[:int(0.8 * n_step)]
        W_fperspan = np.concatenate((W_fperspan, np.zeros(int(0.2 * n_step))))
        return W_fperspan

    def WW(x): # obtain wing weight per span
        WA = A(x)[0] * 9.81 * 173.7434
        return (WA)
    return WW(point) + Wfps(point), A(point)


def get_shear(distribution, limit):
    shear_distribution, shear_error = sp.integrate.quad(distribution, 0, limit)
    return shear_distribution


def get_moment(distribution, limit):
    moment_distribution, moment_error = sp.integrate.quad(distribution, 0, limit)
    return moment_distribution


def get_reaction_moment(distribution, span):
    reaction_moment, moment_error = sp.integrate.quad(distribution, 0, span)
    return reaction_moment


def get_reaction_force(distribution, span):
    reactio_force, moment_error = sp.integrate.quad(distribution, 0, span)
    return reactio_force


# integrate combined load distribution to obtain shear and moment distributions
def get_integrated():
    values_shear = []
    values_moment = []
    combined_load_distribution = sp.interpolate.interp1d(span, combined_load(span), kind='quadratic',fill_value='extrapolate')
    # until pylon
    for i in span[:int(0.32 * n_step)]:
        values_shear.append(get_shear(combined_load_distribution, i))
    # pylon
    for i in span[int(0.32 * n_step):int(0.36 * n_step)]:
        values_shear.append(get_shear(combined_load_distribution, i))
    # until end of fuel tank
    for i in span[int(0.36 * n_step):int(0.8 * n_step)]:
        values_shear.append(get_shear(combined_load_distribution, i))
    # after end of fuel tank
    for i in span[int(0.8 * n_step):]:
        values_shear.append(get_shear(combined_load_distribution, i))
    get_reaction_force(combined_load_distribution, span[-1])

    new_values_shear = [i * -1 + get_reaction_force(combined_load_distribution, span[-1]) for i in values_shear]
    shear_function = sp.interpolate.interp1d(span, new_values_shear, kind='quadratic', fill_value='extrapolate')

    # until pylon
    for i in span[:int(0.32 * n_step)]:
        values_moment.append(get_moment(shear_function, i))
    # pylon
    for i in span[int(0.32 * n_step):int(0.36 * n_step)]:
        values_moment.append(get_moment(shear_function, i))
    # until end of fuel tank
    for i in span[int(0.36 * n_step):int(0.8 * n_step)]:
        values_moment.append(get_moment(shear_function, i))
    # after end of fuel tank
    for i in span[int(0.8 * n_step):]:
        values_moment.append(get_moment(shear_function, i))

    new_values_moment = [i - get_reaction_moment(shear_function, span[-1]) for i in values_moment]
    return new_values_shear, new_values_moment


# get engine weight distribution
def engine_weight():
    a = 0
    for i in span:
        if i in pylon:
            engine_load[a] = 32480 / 2
            a += 1
        else:
            engine_load[a] = 0
            a += 1
    return engine_load


def engine_torque():
    trust = 191270
    pylon_length = 0.725
    a = 0
    for i in span:
        if i in pylon:
            engine_torques[a] = trust * pylon_length
            a += 1
        else:
            engine_torques[a] = 0
            a += 1
    return engine_torques


# get wingbox stiffness and deflection
def get_stiffness():
    def chord(a):
        c = c_root - (c_root - c_tip) * 2 * a / 2 * wing_length
        return (c)

    def d(n):
        distance = 0
        if n == 1:
            return (spar_height / 2) ** 2
        elif n % 2 == 0:
            for k in range(int(n / 2)):
                distance = distance + ((spar_height / 2) ** 2 + (width / 2 - width / n * k) ** 2)
            distance = 2 * distance
        else:
            for k in range(int((n - 1) / 2)):
                distance = distance + ((spar_height / 2) ** 2 + (width / 2 - width / n * k) ** 2)
            distance = distance * 2 + (spar_height / 2) ** 2
        return distance

    def stringers_top(shape, value):
        n = shape.shape
        string_top = np.full(n, value)
        return string_top

    def stringers_bottom(shape, value):
        n = shape.shape
        string_bot = np.full(n, value)
        return string_bot

    # moment of inertia and polar moment of inertia
    def I_xx(y):
        w = width * chord(y)
        h = spar_height * chord(y)
        I = S_stringers * (h / 2) ** 2 * (
                stringers_top(span, stringertop) + stringers_bottom(span,
                                                                    stringerbot)) + 1 / 12 * t2 * h ** 3 + 2 * w * t1 * (
                    h / 2) ** 2
        return (I)

    def J_z(y):
        w = width * chord(y)
        h = spar_height * chord(y)
        J = S_stringers * ((chord(y)) ** 2 * (d(stringertop) + d(stringerbot))) + 4 * (w * h) ** 2 / (
                    (2 * w / t1) + (2 * h / t2))
        return J

    def v(dist):  # get second derivative of deflection
        second = []
        for i in get_integrated()[1]:
            a = 0
            second.append((-1 * i) / (I_xx(dist)[a] * E))
            a += 1
        return second

    def integral(function, limit):
        deflection_distribution, deflection_error = sp.integrate.quad(function, 0, limit)
        return deflection_distribution

    # get deflection function
    def get_deflection():
        values_deflection = []
        values_deflection_actual = []
        deflection_function = sp.interpolate.interp1d(span, v(span), kind='quadratic', fill_value='extrapolate')
        # until pylon
        for i in span[:int(0.32 * n_step)]:
            values_deflection.append(integral(deflection_function, i))
        # pylon
        for i in span[int(0.32 * n_step):int(0.36 * n_step)]:
            values_deflection.append(integral(deflection_function, i))
        # until end of fuel tank
        for i in span[int(0.36 * n_step):int(0.8 * n_step)]:
            values_deflection.append(integral(deflection_function, i))
        # after end of fuel tank
        for i in span[int(0.8 * n_step):]:
            values_deflection.append(integral(deflection_function, i))

        deflection_prime = sp.interpolate.interp1d(span, values_deflection, kind='quadratic', fill_value='extrapolate')

        # until pylon
        for i in span[:int(0.32 * n_step)]:
            values_deflection_actual.append(integral(deflection_prime, i))
        # pylon
        for i in span[int(0.32 * n_step):int(0.36 * n_step)]:
            values_deflection_actual.append(integral(deflection_prime, i))
        # until end of fuel tank
        for i in span[int(0.36 * n_step):int(0.8 * n_step)]:
            values_deflection_actual.append(integral(deflection_prime, i))
        # after end of fuel tank
        for i in span[int(0.8 * n_step):]:
            values_deflection_actual.append(integral(deflection_prime, i))

        return values_deflection_actual

    # wingbox validity and weight
    V_total = 2 * t1 * (chord(0) + chord(2 * wing_length / 2)) / 2 * width * 2 * wing_length / 2 + 2 * t2 * (
            chord(0) + chord(2 * wing_length / 2)) / 2 * spar_height * 2 * wing_length / 2 + (
                          stringertop + stringerbot) * 2 * wing_length / 2 * S_stringers
    mass = V_total * 2 * rho
    return J_z(span), get_deflection(), mass


def combined_load(dist):  # get combined load distribution
    combined = q_scale * get_aerodynamic(dist)[0] - get_inertial(dist)[0] - engine_weight()
    return combined


def combined_torque(dist):  # get combined torque distribution
    combined_torque = q_scale * get_aerodynamic(dist)[1] - engine_torque()
    return combined_torque


# get twist and twist angle distributions
def get_twist(dist):
    torque_function = combined_torque(dist)
    twist_values = []
    J_z = get_stiffness()[0]
    d_theta = []

    def get_integral(func, lim):
        integral, err_tw = sp.integrate.quad(func, 0, lim)
        return integral

    for i in range(0, n_step):
        val = (torque_function[i] / (J_z[i] * G))
        d_theta.append(val)

    d_theta_function = sp.interpolate.interp1d(span, d_theta, kind='quadratic', fill_value='extrapolate')

    # until pylon
    for i in span[:int(0.32 * n_step)]:
        twist_values.append(get_integral(d_theta_function, i))
    # pylon
    for i in span[int(0.32 * n_step):int(0.36 * n_step)]:
        twist_values.append(get_integral(d_theta_function, i))
    # until end of fuel tank
    for i in span[int(0.36 * n_step):int(0.8 * n_step)]:
        twist_values.append(get_integral(d_theta_function, i))
    # after end of fuel tank
    for i in span[int(0.8 * n_step):]:
        twist_values.append(get_integral(d_theta_function, i))

    return twist_values


# plotting loads distribution, shear, moment, torque and twist angle
fig, axs = plt.subplots(3, 2)
plt.tight_layout()
x = span
axs[0, 0].plot(x, combined_load(x), color='red')
axs[0, 0].set_xlabel('Spanwise location [m]')
axs[0, 0].set_ylabel('Load distribution [N/m]')
axs[0, 0].set_title('Load distribution')
axs[1, 0].plot(x, get_integrated()[0], color='lime')
axs[1, 0].set_title('Shear')
axs[1, 0].set_xlabel('Spanwise location [m]')
axs[1, 0].set_ylabel('Shear [N]')
axs[2, 0].plot(x, get_integrated()[1], color='blue')
axs[2, 0].set_xlabel('Spanwise location [m]')
axs[2, 0].set_ylabel('Moment [Nm]')
axs[2, 0].set_title('Moment')
axs[0, 1].plot(x, combined_torque(x), color='cyan')
axs[0, 1].set_xlabel('Spanwise location [m]')
axs[0, 1].set_ylabel('Torque [Nm]')
axs[0, 1].set_title('Torque')
axs[1, 1].plot(x, get_twist(x), color='pink')
axs[1, 1].set_xlabel('Spanwise location [m]')
axs[1, 1].set_ylabel('Twist Angle [rad]')
axs[1, 1].set_title('Twist Angle')
axs[-1, -1].axis('off')
plt.show()

# obtaining maximum angle of twist, deflection and weight of the wingbox
result1 = get_twist(span)
result2 = get_stiffness()
print(2.729 * 1.5 * result1[-1] * 180 / math.pi, "°")
print(2.729 * 1.5 * result2[1][-1], "m")
print(result2[2], "kg")
load_case_1 = 2.729 * 1.5
load_case_2 = -1.5

# plotting twist angle and deflection for load cases
fig, axs = plt.subplots(2, 2)
plt.tight_layout()
x = span
axs[0, 0].plot(x, np.multiply(get_twist(x), load_case_1 * 180 / math.pi), color='red')
axs[0, 0].set_title('Twist Angle')
axs[0, 0].set_xlabel('Spanwise location [m]')
axs[0, 0].set_ylabel('Twist Angle [°]')
axs[1, 0].plot(x, np.multiply(get_twist(x), load_case_2 * 180 / math.pi), color='red')
axs[1, 0].set_title('Twist Angle')
axs[1, 0].set_xlabel('Spanwise location [m]')
axs[1, 0].set_ylabel('Twist Angle [°]')
axs[0, 1].plot(x, np.multiply(get_stiffness()[1], load_case_1), color='lime')
axs[0, 1].set_title('Deflection')
axs[0, 1].set_xlabel('Spanwise location [m]')
axs[0, 1].set_ylabel('Deflection [m]')
axs[1, 1].plot(x, np.multiply(get_stiffness()[1], load_case_2), color='lime')
axs[1, 1].set_title('Deflection')
axs[1, 1].set_xlabel('Spanwise location [m]')
axs[1, 1].set_ylabel('Deflection [m]')
plt.show()
