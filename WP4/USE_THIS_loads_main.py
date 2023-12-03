import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy import integrate


# constants
span = np.linspace(0, 22.445, 50)
pylon = span[16:18]
q_scale = 625
G = 26 * 10 ** 9


# obtaining aerodynamics load from simulations
def get_aerodynamic(x):
    # constants for obtaining aerodynamic load
    q = 61.25
    c_y = [9.84, 2.45]
    y = [0, 22.445]
    cL_0 = 0.149649
    cL_10 = 0.907325
    Cm_0 = -0.24373
    Cm_10 = -1.17556
    Cd_0 = 0.001000
    Cd_10 = 0.036043
    alpha = 0

    c_y_func = sp.interpolate.interp1d(y, c_y, kind='linear', fill_value="extrapolate")

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
        return Cl_func * q * c_y_func(y)

    def M_prime_func(y, Cm_func):
        return Cm_func * q * c_y_func(y) ** 2

    def D_prime_func(y, Cd_func):
        return Cd_func * q * c_y_func(y)

    def find_cl_d(alpha):
        coef = math.sin(math.radians(alpha)) / math.sin(math.radians(10))
        Cl_d = coef * (cL_10 - cL_0) + cL_0
        return Cl_d

    def find_Cl_alpha(y, alpha):
        Cl_d = find_cl_d(alpha)
        Cl_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                                load_file(file_names[0])[3])[0]
        Cl_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                                 load_file(file_names[1])[3])[0]
        return Cl_func_0(y) + ((Cl_d - cL_0) / (cL_10 - cL_0)) * (Cl_func_10(y) - Cl_func_0(y))

    def find_cm_d(alpha):
        return math.sin(math.radians(alpha)) / math.sin(math.radians(10)) * (Cm_10 - Cm_0) + Cm_0

    def find_Cm_alpha(y, alpha):
        Cm_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                                load_file(file_names[0])[3])[2]
        Cm_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                                 load_file(file_names[1])[3])[2]
        return Cm_func_0(y) + ((find_cm_d(alpha) - Cm_0) / (Cm_10 - Cm_0)) * (Cm_func_10(y) - Cm_func_0(y))

    def find_cd_d(alpha):
        return Cd_0 + ((Cd_10 - Cd_0) / (cL_10 ** 2 - cL_0 ** 2) * (find_cl_d(alpha) ** 2 - cL_0 ** 2))

    def find_Cd_alpha(y, alpha):
        Cd_func_0 = interpolate(load_file(file_names[0])[0], load_file(file_names[0])[1], load_file(file_names[0])[2],
                                load_file(file_names[0])[3])[1]
        Cd_func_10 = interpolate(load_file(file_names[1])[0], load_file(file_names[1])[1], load_file(file_names[1])[2],
                                 load_file(file_names[1])[3])[1]
        return Cd_func_0(y) + ((Cd_func_10(y) - Cd_func_0(y)) / (Cd_10 - Cd_0)) * (find_cd_d(alpha) - Cd_0)

    aerodynamic_output = L_prime_func(x, find_Cl_alpha(x, alpha))  # edit angle of attack
    torque_output = M_prime_func(x, find_Cm_alpha(x, alpha))
    drag_output = D_prime_func(x, find_Cd_alpha(x, alpha))

    return aerodynamic_output, torque_output, drag_output


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
        return 0.95 * height * length, height, length

    def Wfps(x):
        W_fperspan = (9.81 * A(x)[0] * 800)[:40]
        W_fperspan = np.concatenate((W_fperspan, np.zeros(10)))
        return W_fperspan

    def WW(x):
        WA = A(x)[0] * 9.81 * 173.7434
        return (WA)

    inertial_loading = WW(point) + Wfps(point)

    return inertial_loading, A(point)


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


def get_integrated_idiot():
    values_shear = []
    values_moment = []
    combined_load_distribution = sp.interpolate.interp1d(span, combined_load(span), kind='quadratic',
                                                         fill_value='extrapolate')
    # 0 - 16
    for i in span[:16]:
        values_shear.append(get_shear(combined_load_distribution, i))
    # 16 - 18
    for i in span[16:18]:
        values_shear.append(get_shear(combined_load_distribution, i))
    # 18 - 40
    for i in span[18:40]:
        values_shear.append(get_shear(combined_load_distribution, i))
    # 40 - 50
    for i in span[40:50]:
        values_shear.append(get_shear(combined_load_distribution, i))
    get_reaction_force(combined_load_distribution, span[-1])

    new_values_shear = [i * -1 + get_reaction_force(combined_load_distribution, span[-1]) for i in values_shear]
    shear_function = sp.interpolate.interp1d(span, new_values_shear, kind='quadratic', fill_value='extrapolate')
    for i in span[:16]:
        values_moment.append(get_moment(shear_function, i))
        # 16 - 18
    for i in span[16:18]:
        values_moment.append(get_moment(shear_function, i))
        # 18 - 40
    for i in span[18:40]:
        values_moment.append(get_moment(shear_function, i))
        # 40 - 50
    for i in span[40:50]:
        values_moment.append(get_moment(shear_function, i))

    new_values_moment = [i - get_reaction_moment(shear_function, span[-1]) for i in values_moment]
    return new_values_shear, new_values_moment


def engine_weight():
    values = np.zeros_like(span)
    a = 0
    for i in span:
        if i in pylon:
            values[a] = 32480 / 2
            a += 1
        else:
            values[a] = 0
            a += 1
    return values


def engine_torque():
    values = np.zeros_like(span)
    trust = 191270
    pylon_lenght = 0.725
    a = 0
    for i in span:
        if i in pylon:
            values[a] = trust * pylon_lenght
            a += 1
        else:
            values[a] = 0
            a += 1
    return values


def get_stiffness():
    # given
    spar_height = 0.0942  # m
    width = 0.5  # m
    c_root = 10.10  # m
    c_tip = 2.73  # m
    span_elliot = 44.89  # m
    E = 68.9 * 10 ** 9  # m
    rho = 2700
    nmax = 2.729
    nmin = -1
    nsafety = 1.5
    # Variables
    N_stringers_top = [4]
    N_stringers_bottom = [4]
    Stringer_change = [0]
    S_stringers = 0.001  # m^2
    t1 = 0.001  # m
    t2 = 0.002

    # chord and distance
    def chord(a):
        c = c_root - (c_root - c_tip) * 2 * a / span_elliot
        return (c)

    """""
    def d(function):
        distance = 0
        for i in span:
            n = function
            if n[i] == 1:
                return (spar_height / 2)
            elif n[i] % 2 == 0:
                for k in range(int(n[i] / 2)):
                    distance = distance + ((spar_height / 2) ** 2 + (width / 2 - width / n[i] * k) ** 2) ** 0.5
                distance = 2 * distance
            else:
                for k in range(int((n - 1) / 2)):
                    distance = distance + ((spar_height / 2) ** 2 + (width / 2 - width / n[i] * k) ** 2) ** 0.5
                distance = distance * 2 + spar_height / 2
        return distance

    # Number of stringers :
    nstringertop = sp.interpolate.interp1d(Stringer_change, N_stringers_top, kind="previous", fill_value="extrapolate")
    nstringerbot = sp.interpolate.interp1d(Stringer_change, N_stringers_bottom, kind="previous",
                                           fill_value="extrapolate")
    """""

    def stringers_top(shape, value):
        n = shape.shape
        string_top = np.full(n, value)
        return string_top

    def stringers_bottom(shape, value):
        n = shape.shape
        string_bot = np.full(n, value)
        return string_bot

    # Stiffness
    def I_xx(y):
        w = width * chord(y)
        h = spar_height * chord(y)
        I = S_stringers * (h / 2) ** 2 * (
                    stringers_top(span, 4) + stringers_bottom(span, 4)) + 1 / 12 * t2 * h ** 3 + 2 * w * t1 * (
                    h / 2) ** 2
        return (I)

    def J_z(y):
        w = width * chord(y)
        h = spar_height * chord(y)
        J = S_stringers * (
                    (chord(y)) ** 2 * (stringers_top(span, 4) ** 2 + stringers_bottom(span, 4) ** 2)) + 2 * w * t1 * (
                    w ** 2 + t1 ** 2) / 12 + 2 * h * t2 * (h ** 2 + t2 ** 2) / 12 + 2 * w * t1 * (
                    h / 2) ** 2 + 2 * h * t2 * (w / 2) ** 2
        return J

    # Intigrals
    def v(y):
        second = []
        for i in get_integrated_idiot()[1]:
            a = 0
            second.append((-1 * i) / (I_xx(y)[a] * E))
            a += 1
        return second

    def intigral(function, limit):
        deflection_distribution, deflection_error = sp.integrate.quad(function, 0, limit)
        return deflection_distribution

    def get_deflection():
        values_deflection = []
        values_deflection_actual = []
        deflection_function = sp.interpolate.interp1d(span, v(span), kind='quadratic', fill_value='extrapolate')
        # 0 - 16
        for i in span[:16]:
            values_deflection.append(intigral(deflection_function, i))
        # 16 - 18
        for i in span[16:18]:
            values_deflection.append(intigral(deflection_function, i))
        # 18 - 40
        for i in span[18:40]:
            values_deflection.append(intigral(deflection_function, i))
        # 40 - 50
        for i in span[40:50]:
            values_deflection.append(intigral(deflection_function, i))

        deflection_prime = sp.interpolate.interp1d(span, values_deflection, kind='quadratic', fill_value='extrapolate')
        # 0 - 16
        for i in span[:16]:
            values_deflection_actual.append(intigral(deflection_prime, i))
        # 16 - 18
        for i in span[16:18]:
            values_deflection_actual.append(intigral(deflection_prime, i))
        # 18 - 40
        for i in span[18:40]:
            values_deflection_actual.append(intigral(deflection_prime, i))
        # 40 - 50
        for i in span[40:50]:
            values_deflection_actual.append(intigral(deflection_prime, i))

        return values_deflection_actual

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
        """

    # wingbox validity and weight
    V_total = 2 * t1 * (chord(0) + chord(span_elliot / 2)) / 2 * width * span_elliot / 2 + 2 * t2 * (
            chord(0) + chord(span_elliot / 2)) / 2 * spar_height * span_elliot / 2
    l = 0
    for i in range(int((len(Stringer_change) - 1))):
        l = l + 1
        V_total = V_total + (stringers_bottom(span, 4) + stringers_top(span, 4)) * (
                Stringer_change[i + 1] - Stringer_change[i]) * S_stringers
    V_total = V_total + (stringers_top(span, 4) + stringers_bottom(span, 4)) * (
                span_elliot / 2 - Stringer_change[-1]) * S_stringers



    return J_z(span)


def combined_load(x):
    combined = q_scale * get_aerodynamic(x)[0] - get_inertial(x)[0] - engine_weight()
    return combined


def combined_torque(x):
    combined_torque = q_scale * get_aerodynamic(x)[1] - engine_torque()
    return combined_torque


def get_twist(x):
    torque_function = combined_torque(x)
    twist_values = []
    thickness = 0.002
    J_z = get_stiffness()
    area = get_inertial(span)[1][0]
    height = get_inertial(span)[1][1]
    length = get_inertial(span)[1][2]
    d_theta = []

    def get_integral(func, lim):
        integral, err_tw = sp.integrate.quad(func, 0, lim)
        return integral

    for i in range(0, 50):
        val = (torque_function[i] / (4 * (area[i] ** 2) * G)) * (2 * (length[i] + height[i]) / thickness)
        d_theta.append(val)

    d_theta_function = sp.interpolate.interp1d(span, d_theta, kind='quadratic', fill_value='extrapolate')

    for i in span[:16]:
        twist_values.append(get_integral(d_theta_function, i))
    # 16 - 18
    for i in span[16:18]:
        twist_values.append(get_integral(d_theta_function, i))
    # 18 - 40
    for i in span[18:40]:
        twist_values.append(get_integral(d_theta_function, i))
    # 40 - 50
    for i in span[40:50]:
        twist_values.append(get_integral(d_theta_function, i))

    return twist_values


# plotting the graphs
fig, axs = plt.subplots(3, 2)
plt.tight_layout()
x = span
axs[0, 0].plot(x, combined_load(x), color='red')
axs[0, 0].set_title('Load distribution')
axs[1, 0].plot(x, get_integrated_idiot()[0], color='lime')
axs[1, 0].set_title('Shear')
axs[2, 0].plot(x, get_integrated_idiot()[1], color='blue')
axs[2, 0].set_title('Moment')
axs[0, 1].plot(x, combined_torque(x), color='cyan')
axs[0, 1].set_title('Torque')
axs[1, 1].plot(x, get_twist(x), color='pink')
axs[1, 1].set_title('Twist Angle')
axs[-1, -1].axis('off')
plt.show()



