import numpy as np
import scipy as sp


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



