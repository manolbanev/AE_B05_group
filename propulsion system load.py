import numpy as np
import scipy as sp
from scipy import interpolate

# forces undder
forces = np.array([191270, 0, -32480.91])

position = np.array([-1.45, 8.53, 1.7])


def calc_distance(point):
    distance = np.subtract(point, position)
    return distance


def calc_moment(force, pos):
    mom = np.cross(force, calc_distance(pos))
    return mom






