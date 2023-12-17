from WP4.USE_THIS_loads_main import get_integrated
from WP5.wing import wing1
from WP5.wingbox import wingbox1
from WP5.ribs import rib1
import matplotlib.pyplot as plt

sigma_yield = 240 * 1e6


def get_shear(): # change calculation of the stress in spars using appendix f or sth
    return get_integrated()[0]


def get_moment():
    return get_integrated()[1]


def get_sigma():
    sigma_list = []
    moments = get_moment()
    a = 0
    for i in wing1.span:
        sigma = (moments[a] * (wingbox1.get_height(i) / 2)) / wingbox1.get_moment_of_inertia(i)
        sigma_list.append(sigma)
        a += 1
    return sigma_list


def get_stringer_buckling():
    stringer_buckling = []
    for i in wing1.span:
        stringer_buckling.append(wingbox1.stringer.get_buckling(rib1.get_distance(i)))
    return stringer_buckling


def get_web_buckling():
    web_buckling = []
    for i in wing1.span:
        web_buckling.append(wingbox1.spar_rear.get_buckling(i))
    return web_buckling


def check_stringer_failure():
    sigma = []
    for i in get_sigma():
        sigma.append(abs(i))
    stringer_buckling = get_stringer_buckling()
    for i in range(len(wing1.span)):
        if sigma[i] > stringer_buckling[i]:
            print('Stringer buckling at spanwise location: ', wing1.span[i])
    print('Done checking stringer failure')


def check_spar_failure():
    tau = get_shear()
    web_buckling = get_web_buckling()
    for i in range(len(wing1.span)):
        if tau[i] > web_buckling[i]:
            print('Web buckling at spanwise location: ', wing1.span[i])
    print('Done checking web failure')


def check_compressive_failure():
    sigma = []
    for i in get_sigma():
        sigma.append(abs(i))
    for i in range(len(wing1.span)):
        if sigma[i] > sigma_yield:
            print('Compressive failure at spanwise location: ', wing1.span[i])
    print('Done checking compressive failure')


check_stringer_failure()
check_spar_failure()
check_compressive_failure()


figs, axs = plt.subplots(3)
plt.tight_layout()
axs[0].plot(wing1.span, get_sigma(), color='magenta')
axs[0].set_title('Sigma')
axs[1].plot(wing1.span, get_stringer_buckling(), color='cyan')
axs[1].set_title('Stringer Buckling')
axs[2].plot(wing1.span, get_web_buckling(), color='lime')
axs[2].set_title('Web Buckling')
plt.show()



# implement (in)sanity checks here -> a/b ratio, size of the stringers
