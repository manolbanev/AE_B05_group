import numpy as np
import matplotlib.pyplot as plt

c_root = 10.10  # m
c_tip = 2.73    # m
span = np.linspace(0, 22.445, 1000)

def chord_distribution(span, c_root, c_tip, wingspan):
    return c_root - ((c_root - c_tip) / wingspan) * span

def wingbox_area_distribution(span):
    chord_at_x = chord_distribution(span, c_root, c_tip, wingspan=22.445)
    height = 0.0942 * chord_at_x
    length = 0.5 * chord_at_x
    return height * length

def wingbox_weight_distribution(span):
    area_distribution = wingbox_area_distribution(span)
    weight_distribution = area_distribution * 9.81 * 840
    return weight_distribution

def total_weight(span):
    total_weight = np.trapz(wingbox_weight_distribution(span), span)
    return total_weight

# Calculate the weight distribution
weight_distribution = wingbox_weight_distribution(span)

fuel_scaling = 51500/397925.7085
wing_scaling = 156381.21/397925.7085

# Add an arrow at a specific point (e.g., at span = 10)
arrow = plt.arrow(3.367, 0, 0, -32481, width=0.4, length_includes_head=True, head_width=0.9, head_length=8120, label = "Engine Weight")


# Plot the results
plt.plot(span, fuel_scaling*-weight_distribution, label='Fuel Weight Distribution', color = "red")
plt.plot(span, wing_scaling*-weight_distribution, label='Wing Structural Weight Distribution', color = "green")
plt.xlabel('Span (m)')
plt.ylabel('Weight Distribution (N/m)')
plt.title('Wingbox Weight Distribution Along the Span')
plt.legend()
plt.show()

# Calculate and print the total weight of the fuel
total_fuel = fuel_scaling*total_weight(span)
print(f"Total Weight of the Fuel: {total_fuel:.0f} Newtons")

total_wing = wing_scaling*total_weight(span)
print(f"Total Weight of the Fuel: {total_wing:.0f} Newtons")