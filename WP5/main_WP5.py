import numpy as np
from WP5 import compressive_strenght
from WP4.USE_THIS_loads_main import combined_load


def get_shear(dist):
    return combined_load(dist)


