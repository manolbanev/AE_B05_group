from WP4.USE_THIS_loads_main import combined_load
import numpy as np
import scipy as sp

span = np.linspace(0, 22.445, 50)
stress = []


def calculate_stress(force: float, area: float) -> float:
    """
    This function returns the stress given the force and the area at a point
    :param force: applied force
    :param area: area of the cross-section
    :return: stress value
    """
    return force / area


def check_yield(stress: float, yield_stress: float) -> bool:
    """
    This function returns True if the structure fails at a point and False if it does not
    :param stress: applied stress
    :param yield_stress: yield stress of the material
    :return: True or False
    """
    if abs(stress) > yield_stress:
        return True
    elif abs(stress) < yield_stress:
        return False


def get_shear(dist: np.array) -> any:
    """
    This function returns the shear force as a function of distance from the root.
    :param dist: x-axis for the function
    :return: shear force function
    """
    shear_lst = combined_load(dist)
    shear_func = sp.interpolate.interp1d(dist, shear_lst, kind='cubic', fill_value='extrapolate')
    return shear_func


def stress_distribution():
    for i in span:
        stress.append(calculate_stress(get_shear(span)(i), 1))
    return stress
