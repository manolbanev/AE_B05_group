from WP5.stringer import Stringer
from WP5.spar import Spar
from WP5.wing import wing1
import numpy as np


class Wingbox:
    def __init__(self, skin_thickness, shear_factor=1.5, skin_k_factor=4):
        self.skin_thickness = skin_thickness
        self.shear_factor = shear_factor
        self.skin_k_factor = skin_k_factor
        self.spar_front = Spar(0.01, 0.1, 9, 68.9 * 1e9, 0.33)
        self.spar_rear = Spar(0.01, 0.6, 9, 68.9 * 1e9, 0.33)
        self.stringer = Stringer(0.007, 0.05, 0.08, 20, 0.25, 68.9 * 1e9)

    def get_area(self, y):
        return self.spar_rear.get_height(y) * (self.spar_rear.get_location(y) - self.spar_front.get_location(y))

    def get_height(self, y):
        return self.spar_rear.get_height(y)

    def get_moment_of_inertia(self, y):
        stringers_moi = self.stringer.get_moment_of_inertia()
        stringers_steiner = self.stringer.area * ((self.spar_rear.get_height(y) / 2) ** 2)
        spars_moi = self.spar_rear.get_moment_of_inertia(y) + self.spar_front.get_moment_of_inertia(y)
        skin_steiner = ((self.spar_rear.position - self.spar_front.position) * wing1.get_chord(y)
                        * self.skin_thickness * (self.spar_rear.get_height(y) / 2) ** 2)
        return self.stringer.number * (stringers_moi + stringers_steiner) + spars_moi + (2 * skin_steiner)

    def get_tau(self, point, shear):
        tau_av = (shear / (self.spar_rear.thickness * self.spar_rear.get_height(point) +
                           self.spar_front.thickness * self.spar_front.get_height(point)))

        tau_max = tau_av * self.shear_factor
        return tau_max

    def get_skin_buckling(self, y):
        sigma_crit = (np.pi ** 2 * self.skin_k_factor * self.spar_front.young_modulus
                      / (12 * (1 - self.spar_front.poisson_ratio ** 2)) *
                      (self.skin_thickness * (self.stringer.number + 1) / wing1.get_chord(y)) ** 2)
        return sigma_crit


wingbox1 = Wingbox(0.01)
