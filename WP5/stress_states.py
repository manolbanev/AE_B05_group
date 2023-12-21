from WP4.USE_THIS_loads_main import get_integrated, combined_torque
from WP5.components import wing1, wingbox1, rib1
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
        sigma_list.append(sigma * 1e-6)
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
        tau_list.append((wingbox1.get_tau(i, shear[a]) + shear_flow[a] * wingbox1.spar_rear.thickness) * 1e-6)
        a += 1
    return tau_list


def get_stringer_buckling():
    stringer_buckling = []
    for i in wing1.span:
        stringer_buckling.append(wingbox1.stringer.get_buckling(rib1.get_distance(i)) * 1e-6)
    return stringer_buckling


def get_web_buckling():
    web_buckling = []
    for i in wing1.span:
        web_buckling.append(wingbox1.spar_rear.get_buckling(i) * 1e-6)
    return web_buckling


def get_skin_buckling_function():
    skin_buckling = []
    for i in wing1.span:
        skin_buckling.append(wingbox1.get_skin_buckling(i) * 1e-6)
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


def get_safety_margin():
    sm_stringer = []
    sm_skin = []
    sm_web = []
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


def check_failure():
    check_stringer_failure()
    check_spar_failure()
    check_compressive_failure()
    check_skin_failure()
