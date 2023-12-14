from WP4.USE_THIS_loads_main import get_integrated
from WP5.wing import wing1
from WP5.wingbox import Wingbox

wingbox1 = Wingbox(0.003)


def get_shear():
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


print(max(get_shear()))
print(get_sigma())
