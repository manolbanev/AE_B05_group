import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy import integrate

# constants
r_integration_step = 50
r_span = np.linspace(0, 22.445, r_integration_step)
r_tip_chord = 2.73
r_root_chord = 10.10
r_span_length = 22.445

class Wing:
    def __init__(self, span, tip_chord, root_chord, span_length, integration_step):
        self.span = span
        self.tip_chord = tip_chord
        self.root_chord = root_chord
        self.span_length = span_length
        self.integration_step = integration_step

    def get_chord_distribution(self, x):
        return self.root_chord - ((self.root_chord - self.tip_chord) / self.span_length) * x


class Force:
    def __init__(self, values_list, interpolate_kind, span, tip_chord, root_chord, span_length):
        self.values_list = values_list
        self.interpolate_kind = interpolate_kind
        wing1 = Wing(span, tip_chord, root_chord, span_length)

    def get_function(self, x_values):
        return sp.interpolate.interp1d(x_values, self.values_list, kind=self.interpolate_kind, fill_value='extrapolate')

    def integrate(self, x_values):
        for i in x_values[:self.wing1.integration_step * 0.32]:
            self.values_list.append(sp.integrate.quad(self.get_function(self, x_values), 0, i)[0])
        for i in x_values[self.wing1.integration_step * 0.32:self.wing1.integration_step * 0.36]:
            self.values_list.append(sp.integrate.quad(self.get_function(self, x_values), 0, i)[0])
        for i in x_values[self.wing1.integration_step * 0.36:self.wing1.integration_step * 0.8]:
            self.values_list.append(sp.integrate.quad(self.get_function(self, x_values), 0, i)[0])
        for i in x_values[self.wing1.integration_step * 0.8:self.wing1.integration_step]:
            self.values_list.append(sp.integrate.quad(self.get_function(self, x_values), 0, i)[0])
        return self.values_list

aerodynamic_list = []

Aerodynamic = Force(aerodynamic_list, 'cubic', r_span, r_tip_chord, r_root_chord, r_span_length)


