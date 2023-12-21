from WP4.USE_THIS_loads_main import get_integrated
import numpy as np
import scipy as sp


class Rib:
    def __init__(self, number):
        self.moment = get_integrated()[1]
        self.sigma = []
        self.positions = [0]
        self.number = number

    def get_spacing(self):
        distance = wing1.length / self.number
        return distance

    def get_positions(self):
        spacing = self.get_spacing()
        for i in range(len(wing1.span)):
            if wing1.span[i] - self.positions[-1] >= spacing:
                self.positions.append(self.positions[-1] + spacing)
                if self.positions[-1] >= wing1.length:
                    self.positions[-1] = wing1.length
                    break
        return self.positions

    def get_distance(self, point):      # works for an integer might be an issue with passing the list in
        self.positions = self.get_positions()
        if point <= self.positions[-2]:
            distance = self.get_spacing()
            return distance
        else:
            return 22.445 - self.positions[-2]


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


class Wing:
    def __init__(self, length, c_tip, c_root, surface):
        self.length = length
        self.c_tip = c_tip
        self.c_root = c_root
        self.surface = surface
        self.span = np.linspace(0, self.length, 50)

    def get_chord(self, position):
        return float(self.c_root - ((self.c_root - self.c_tip) / self.length) * position)


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


wing1 = Wing(22.445, 2.73, 10.10, 287.9)
wingbox1 = Wingbox(0.01)
rib1 = Rib(25)
