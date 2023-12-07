import numpy as np
from WP4.USE_THIS_loads_main import combined_load
import scipy as sp


def get_shear(dist: np.array) -> any:
    """
    This function returns the shear force as a function of distance from the root.
    :param dist: x-axis for the function
    :return: shear force function
    """
    shear_lst = combined_load(dist)
    shear_func = sp.interpolate.interp1d(dist, shear_lst, kind='cubic', fill_value='extrapolate')
    return shear_func


def calculate_margin(dist: np.array, load_func: any, failure_func: any) -> any:
    """
    This function returns the margin of safety as a function of distance from the root.
    :param dist: x-axis for the function
    :param load_func: applied load function
    :param failure_func: function of values of load at which failure occurs
    :return: margin of safety function
    """
    margin_lst = []
    for i in dist:
        margin_lst.append(load_func(i) / failure_func(i))
    margin_func = sp.interpolate.interp1d(dist, margin_lst, kind='cubic', fill_value='extrapolate')
    return margin_func



