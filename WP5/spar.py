from WP5.wing import wing1
import numpy as np
import scipy as sp


class Spar:
    def __init__(self, thickness, position, plate_ar_coefficient, young_modulus, poisson_ratio):
        self.thickness = thickness
        self.position = position
        self.plate_ar_coefficient = plate_ar_coefficient
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio

    def get_location(self, y):
        return self.position * wing1.get_chord(y)

    def get_height(self, y):
        return 0.0942 * wing1.get_chord(y)

    def get_moment_of_inertia(self, y):
        return self.thickness * self.get_height(y) ** 3 / 12

    def get_buckling(self, y):
        shear_crit = ((np.pi ** 2 * self.plate_ar_coefficient * self.young_modulus) /
                      (12 * (1 - self.poisson_ratio ** 2)) * (self.thickness / self.get_height(y)) ** 2)
        return shear_crit

    def get_area(self):
        spar_lst = []
        for i in wing1.span:
            spar_lst.append(self.get_height(i))
        spar_func = sp.interpolate.interp1d(wing1.span, spar_lst, kind='linear', fill_value='extrapolate')
        spar_area, spar_error = sp.integrate.quad(spar_func, 0, wing1.length)
        return spar_area

