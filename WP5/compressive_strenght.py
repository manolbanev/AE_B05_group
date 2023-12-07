def calculate_stress(force: float, area: float) -> float:
    """Calculate stress given force and area"""
    return force / area


def check_yield(stress: float, yield_stress: float) -> bool:
    """Check if applied stress is greater than yield stress"""
    if abs(stress) > yield_stress:
        return True
    elif abs(stress) < yield_stress:
        return False

