from WP4.USE_THIS_loads_main import get_integrated
import numpy as np
import scipy as sp

# this file contains the classes for the wingbox, wing, spar, stringer and rib


class Wing:  # class used to represent the wing
    def __init__(self, length: float, c_tip: float, c_root: float, surface: float):
        self.length = length
        self.c_tip = c_tip
        self.c_root = c_root
        self.surface = surface
        self.span = np.linspace(0, self.length, 50)

    def get_chord(self, position: float) -> float:  # get the chord length at a certain position along the span
        return float(self.c_root - ((self.c_root - self.c_tip) / self.length) * position)


class Spar:  # class used to represent the spars
    def __init__(self, thickness: float, position: float,
                 plate_ar_coefficient=9, young_modulus=68.9 * 1e9, poisson_ratio=0.33):
        self.thickness = thickness
        self.position = position * 0.01  # position of the spar along the chord, in percentage
        self.plate_ar_coefficient = plate_ar_coefficient  # the following 3 values are constant for all design options
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio

    def get_location(self, point: float) -> float:  # get the location of the spar along the span in meters
        return self.position * wing1.get_chord(point)

    def get_height(self, point: float) -> float:  # get the height of the spar at a certain point along the span
        return 0.0942 * wing1.get_chord(point)

    def get_moment_of_inertia(self, point):  # get the spars moment of inertia at a certain point along the span
        return self.thickness * self.get_height(point) ** 3 / 12

    def get_buckling(self, point: float) -> float:  # get the critical buckling stress of the spar at a certain point
        shear_crit = ((np.pi ** 2 * self.plate_ar_coefficient * self.young_modulus) /
                      (12 * (1 - self.poisson_ratio ** 2)) * (self.thickness / self.get_height(point)) ** 2)
        return shear_crit

    def get_area(self) -> float:  # get the total area of the spar along the wing
        spar_lst = []
        for i in wing1.span:
            spar_lst.append(self.get_height(i))
        spar_func = sp.interpolate.interp1d(wing1.span, spar_lst, kind='linear', fill_value='extrapolate')
        spar_area, spar_error = sp.integrate.quad(spar_func, 0, wing1.length)
        return spar_area


class Stringer:  # class used to represent the stringers
    def __init__(self, thickness: float, height: float,
                 width: float, number: int, clamp_factor=0.25, young_modulus=68.9 * 1e9):
        self.thickness = thickness
        self.height = height
        self.width = width
        self.number = number
        self.clamp_factor = clamp_factor  # the following 2 values are constant for all design options
        self.young_modulus = young_modulus
        self.area = self.thickness * self.height + self.width * self.thickness

    def get_moment_of_inertia(self) -> float:  # get the stringers moment of inertia
        moment_of_inertia = (self.thickness * (self.height ** 2) * ((self.height / 3) + self.width) / 4)
        return moment_of_inertia

    def get_buckling(self, distance_ribs: float) -> float:  # get the critical buckling stress of the stringers
        sigma_crit = ((self.clamp_factor * np.pi ** 2 * self.young_modulus * self.get_moment_of_inertia()) / (
                (distance_ribs ** 2) * self.area))
        return sigma_crit


class Wingbox:  # class used to represent the wingbox
    def __init__(self, skin_thickness, shear_factor=1.5, skin_k_factor=4):
        self.skin_thickness = skin_thickness
        self.shear_factor = shear_factor
        self.skin_k_factor = skin_k_factor
        self.spar_front = Spar(0.01, 10)  # front spar
        self.spar_rear = Spar(0.01, 60)  # rear spar
        self.stringer = Stringer(0.007, 0.05, 0.08, 20)  # stringers

    def get_area(self, point: float) -> float:  # get the area of the wingbox at a certain point along the span
        return (self.spar_rear.get_height(point) *
                (self.spar_rear.get_location(point) - self.spar_front.get_location(point)))

    def get_height(self, point: float) -> float:  # get height of the wingbox at a point along the span
        return self.spar_rear.get_height(point)

    def get_moment_of_inertia(self, point: float) -> float:  # get the moment of inertia of the wingbox at a point
        stringers_moi = self.stringer.get_moment_of_inertia()
        stringers_steiner = self.stringer.area * ((self.spar_rear.get_height(point) / 2) ** 2)
        spars_moi = self.spar_rear.get_moment_of_inertia(point) + self.spar_front.get_moment_of_inertia(point)
        skin_steiner = ((self.spar_rear.position - self.spar_front.position) * wing1.get_chord(point)
                        * self.skin_thickness * (self.spar_rear.get_height(point) / 2) ** 2)
        return self.stringer.number * (stringers_moi + stringers_steiner) + spars_moi + (2 * skin_steiner)

    def get_tau(self, point: float, shear: float) -> float:  # get the shear stress at a point along the span
        tau_av = (shear / (self.spar_rear.thickness * self.spar_rear.get_height(point) +
                           self.spar_front.thickness * self.spar_front.get_height(point)))

        tau_max = tau_av * self.shear_factor
        return tau_max

    def get_skin_buckling(self, point: float) -> float:  # get the critical buckling stress of the skin
        sigma_crit = (np.pi ** 2 * self.skin_k_factor * self.spar_front.young_modulus
                      / (12 * (1 - self.spar_front.poisson_ratio ** 2)) *
                      (self.skin_thickness * (self.stringer.number + 1) / wing1.get_chord(point)) ** 2)
        return sigma_crit


class Rib:  # class used to represent ribs
    def __init__(self, number: int):
        self.positions = [0]
        self.number = number  # desired number of ribs

    def get_spacing(self) -> float:  # get the required, even spacing for the desired number of ribs
        distance = wing1.length / self.number
        return distance

    def get_positions(self) -> list:  # get the positions of the ribs along the span
        spacing = self.get_spacing()
        for i in range(len(wing1.span)):
            if wing1.span[i] - self.positions[-1] >= spacing:
                self.positions.append(self.positions[-1] + spacing)
                if self.positions[-1] >= wing1.length:
                    self.positions[-1] = wing1.length
                    break
        return self.positions

    def get_distance(self, point) -> float:  # get the distance between two ribs at a certain point along the span
        self.positions = self.get_positions()
        if point <= self.positions[-2]:
            distance = self.get_spacing()
            return distance
        else:
            return 22.445 - self.positions[-2]


# create wingbox for desired design option
wing1 = Wing(22.445, 2.73, 10.10, 287.9)
wingbox1 = Wingbox(0.01)
rib1 = Rib(25)
