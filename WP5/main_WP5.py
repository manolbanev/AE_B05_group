from WP4.USE_THIS_loads_main import get_integrated
from WP4.USE_THIS_loads_main import combined_torque
from WP5.wing import wing1
from WP5.wingbox import wingbox1
from WP5.ribs import rib1
import matplotlib.pyplot as plt
import scipy as sp

sigma_yield = 240 * 1e6


def get_shear():
    return get_integrated()[0]


def get_moment():
    return get_integrated()[1]


def get_torque():
    return combined_torque(wing1.span)


def get_sigma():
    sigma_list = []
    moments = get_moment()
    a = 0
    for i in wing1.span:
        sigma = (moments[a] * (wingbox1.get_height(i) / 2)) / wingbox1.get_moment_of_inertia(i)
        sigma_list.append(sigma)
        a += 1
    return sigma_list


def get_sigma_absolute():
    sigma_list = []
    sigma = get_sigma()
    for i in sigma:
        sigma_list.append(abs(i))
    return sigma_list


def get_shear_flow():
    torque = get_torque()
    shear_flow = []
    a = 0
    for i in wing1.span:
        shear_flow.append(torque[a] / 2 * wingbox1.get_area(i))
        a += 1
    return shear_flow


def get_shear_stress():
    tau_list = []
    shear = get_shear()
    shear_flow = get_shear_flow()
    a = 0
    for i in wing1.span:
        tau_list.append(wingbox1.get_tau(i, shear[a]) + shear_flow[a] * wingbox1.spar_rear.thickness)
        a += 1
    return tau_list


def get_stringer_buckling():
    stringer_buckling = []
    # fix how rib distance is passed here,
    for i in wing1.span:
        stringer_buckling.append(wingbox1.stringer.get_buckling(rib1.get_distance(i)))
    return stringer_buckling


def get_web_buckling():
    web_buckling = []
    for i in wing1.span:
        web_buckling.append(wingbox1.spar_rear.get_buckling(i))
    return web_buckling


def get_skin_buckling_function():
    skin_buckling = []
    for i in wing1.span:
        skin_buckling.append(wingbox1.get_skin_buckling(i))
    return skin_buckling


def check_stringer_failure():
    sigma = get_sigma_absolute()
    stringer_buckling = get_stringer_buckling()
    for i in range(len(wing1.span)):
        if sigma[i] > stringer_buckling[i]:
            print('Stringer buckling at spanwise location: ', wing1.span[i])
    print('Done checking stringer failure')


def check_spar_failure():
    tau = get_shear_stress()
    web_buckling = get_web_buckling()
    for i in range(len(wing1.span)):
        if tau[i] > web_buckling[i]:
            print('Web buckling at spanwise location: ', wing1.span[i])
    print('Done checking web failure')


def check_compressive_failure():
    sigma = get_sigma_absolute()
    for i in range(len(wing1.span)):
        if sigma[i] > sigma_yield:
            print('Compressive failure at spanwise location: ', wing1.span[i])
    print('Done checking compressive failure')


def check_skin_failure():
    sigma = []
    for i in get_sigma():
        sigma.append(abs(i))
    for i in range(len(wing1.span)):
        if sigma[i] > wingbox1.get_skin_buckling(wing1.span[i]):
            print('Skin buckling at spanwise location: ', wing1.span[i])
    print('Done checking skin failure')


check_stringer_failure()
check_spar_failure()
check_compressive_failure()
check_skin_failure()



'''
figs, axs = plt.subplots(5)
plt.tight_layout()
axs[0].plot(wing1.span, get_sigma(), color='magenta')
axs[0].set_title('Sigma')
axs[1].plot(wing1.span, get_stringer_buckling(), color='cyan')
axs[1].set_title('Stringer Buckling')
axs[2].plot(wing1.span, get_web_buckling(), color='lime')
axs[2].set_title('Web Buckling')
axs[3].plot(wing1.span, get_shear_stress(), color='yellow')
axs[3].set_title('Shear Stress')
axs[4].plot(wing1.span, get_skin_buckling_function(), color='purple')
axs[4].set_title('Skin Buckling')
plt.show()
'''


fig, ax = plt.subplots(2)
x = wing1.span
ax[0].plot(x, get_sigma_absolute(), color='magenta')
ax[0].plot(x, get_stringer_buckling(), color='cyan', linestyle='dashed')
ax[0].plot(x, get_skin_buckling_function(), color='purple', linestyle='dashed')
ax[1].plot(x, get_shear_stress(), color='yellow')
ax[1].plot(x, get_web_buckling(), color='lime', linestyle='dashed')
plt.show()


# implement (in)sanity checks here -> a/b ratio, size of the stringers, rib number vs rib thickness
def check_stringer_height():
    for i in wing1.span:
        if wingbox1.spar_rear.get_height(i) <= 2 * wingbox1.stringer.height:
            return True


def check_stringer_width():
    for i in wing1.span:
        if wingbox1.stringer.number * wingbox1.stringer.width / 2 >= wing1.get_chord(i):
            return True


def check_ribs():  # find and swap reference values for rib thickness and rib % in the wing
    if rib1.number * 0.01 > wing1.length / 4:
        return True


def insanity_check():
    if check_stringer_height():
        print('Warning, stringers too tall')
    if check_stringer_width():
        print('Warning, stringers too wide')
    if check_ribs():
        print('Warning, too many ribs')


def get_mass():
    density = 2700
    spars = (wingbox1.spar_rear.get_area() * wingbox1.spar_rear.thickness
             + wingbox1.spar_front.get_area() * wingbox1.spar_front.thickness) * density
    skin = wing1.surface * wingbox1.skin_thickness * density
    stringers = wingbox1.stringer.number * wingbox1.stringer.area * wing1.length * density
    ribs = rib1.number * 1.7 * 0.033 * density / 2
    total = spars + skin + stringers + ribs
    print('Total mass of the wing is: ', total)
    pass


def get_safety_margin():
    sm_stringer = []
    sm_skin = []
    sm_web =[]
    sigmas = get_sigma_absolute()
    str_buckling = get_stringer_buckling()
    skin = get_skin_buckling_function()
    shear = get_shear_stress()
    web_buckling = get_web_buckling()
    for i in range(len(wing1.span)-1):
        sm_stringer.append(str_buckling[i] / sigmas[i])
        sm_skin.append(skin[i] / sigmas[i])
        sm_web.append(web_buckling[i] / shear[i])
    span = wing1.span[:-1]
    sm_str_func = sp.interpolate.interp1d(span, sm_stringer, kind='linear', fill_value='extrapolate')
    sm_skin_func = sp.interpolate.interp1d(span, sm_skin, kind='linear', fill_value='extrapolate')
    sm_web_func = sp.interpolate.interp1d(span, sm_web, kind='linear', fill_value='extrapolate')
    return sm_str_func, sm_skin_func, sm_web_func


insanity_check()
get_mass()

# plot safety margins here
fig, ax = plt.subplots(3)
x = wing1.span
ax[0].plot(x, get_safety_margin()[0](x), color='magenta')
ax[0].set_title('Stringer Safety Margin')
ax[1].plot(x, get_safety_margin()[1](x), color='cyan')
ax[1].set_title('Skin Safety Margin')
ax[2].plot(x, get_safety_margin()[2](x), color='purple')
ax[2].set_title('Web Safety Margin')


ax[0].set_ylim(0, 10)
ax[1].set_ylim(0, 10)
ax[2].set_ylim(0, 10)
ax[0].axhline(y=1, color='black', linestyle='dashed')
ax[1].axhline(y=1, color='black', linestyle='dashed')
ax[2].axhline(y=1, color='black', linestyle='dashed')


plt.show()


