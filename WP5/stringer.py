import numpy as np


class Stringer:
    def __init__(self, thickness, height, width, number, clamp_factor, young_modulus):
        self.thickness = thickness
        self.height = height
        self.width = width
        self.number = number
        self.clamp_factor = clamp_factor
        self.young_modulus = young_modulus
        self.area = self.thickness * self.height + self.width * self.thickness

    def get_moment_of_inertia(self):  # change the method of calculating the moment of inertia
        moment_of_inertia = (self.thickness * (self.height ** 2) * ((self.height / 3) + self.width) / 4)
        return moment_of_inertia

    def get_buckling(self, distance_ribs):
        sigma_crit = ((self.clamp_factor * np.pi ** 2 * self.young_modulus * self.get_moment_of_inertia()) / (
                (distance_ribs ** 2) * self.area))
        return sigma_crit



