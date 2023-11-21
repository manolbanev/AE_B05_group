import matplotlib.pyplot as plt
import numpy as np


def shear(pos):
    value_shear = pos * 2
    return value_shear


def moment(pos):
    value_moment = pos + 1
    return value_moment


def torque(pos):
    value_torque = pos
    return value_torque


x = np.linspace(0, 23)
y_shear = shear(x)
y_moment = moment(x)
y_torque = torque(x)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
plt.subplots_adjust(hspace=0.8)
ax1.plot(x, y_shear, color='red')
ax1.set_title('Shear')
ax2.plot(x, y_moment, color='lime')
ax2.set_title('Moment')
ax3.plot(x, y_torque, color='blue')
ax3.set_title('Torque')
plt.show()
