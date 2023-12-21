import numpy as np


class Wing:
    def __init__(self, length, c_tip, c_root, surface):
        self.length = length
        self.c_tip = c_tip
        self.c_root = c_root
        self.surface = surface
        self.span = np.linspace(0, self.length, 50)

    def get_chord(self, position):
        return float(self.c_root - ((self.c_root - self.c_tip) / self.length) * position)


wing1 = Wing(22.445, 2.73, 10.10, 287.9)
