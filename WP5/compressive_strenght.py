def calculate_stress(force: float, area: float) -> float:
    """
    This function returns the stress given the force and the area at a point
    :param force: applied force
    :param area: area of the cross-section
    :return: stress value
    """
    return force / area


def check_yield(stress: float, yield_stress: float) -> bool:
    """
    This function returns True if the structure fails at a point and False if it does not
    :param stress: applied stress
    :param yield_stress: yield stress of the material
    :return: True or False
    """
    if abs(stress) > yield_stress:
        return True
    elif abs(stress) < yield_stress:
        return False

