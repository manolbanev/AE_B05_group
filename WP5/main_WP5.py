from WP5.components import wing1, wingbox1, rib1
from WP5.stress_states import check_failure
from WP5.plots import plot_safety_margin, plot_failures_combined, plot_failures


def check_stringer_height() -> bool:  # checks if stringers fit in the wingbox (height)
    for i in wing1.span:
        if wingbox1.spar_rear.get_height(i) <= 2 * wingbox1.stringer.height:
            return True


def check_stringer_width() -> bool:  # checks if stringers fit in the wingbox (width)
    for i in wing1.span:
        if wingbox1.stringer.number * wingbox1.stringer.width / 2 >= wing1.get_chord(i):
            return True


def check_ribs() -> bool:  # checks if ribs are too close to each other
    if rib1.number * 0.01 > wing1.length / 4:
        return True


def insanity_check():  # checks if the wingbox is feasible
    if check_stringer_height():
        print('Warning, stringers too tall')
    if check_stringer_width():
        print('Warning, stringers too wide')
    if check_ribs():
        print('Warning, too many ribs')
    else:
        print('Insanity check passed')


def get_mass() -> float:  # calculates the total mass of the wing
    density = 2700
    spars = (wingbox1.spar_rear.get_area() * wingbox1.spar_rear.thickness
             + wingbox1.spar_front.get_area() * wingbox1.spar_front.thickness) * density
    skin = wing1.surface * wingbox1.skin_thickness * density
    stringers = wingbox1.stringer.number * wingbox1.stringer.area * wing1.length * density
    ribs = rib1.number * 1.7 * 0.033 * density / 2
    total = spars + skin + stringers + ribs
    print('Total mass of the wing is: ', total)
    return total


insanity_check()
check_failure()
get_mass()
plot_failures_combined()
plot_safety_margin()
