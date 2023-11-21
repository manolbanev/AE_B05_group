import matplotlib.pyplot as plt
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.axis('equal')
X, Y = [], []
x_1 = 0.0015
x_2 = 0.06
delta_y = 0.0312
chord = 4.8

for line in open('NASA_sc_20412.dat', 'r'):
  values = [float(s) for s in line.split()]
  X.append(values[0])
  Y.append(values[1])
x_smallest = min(X)
x_smallest_pos = X.index(x_smallest)
x_biggest = max(X)
x_biggest_pos = X.index(x_biggest)
y_max = max(Y)
y_min = min(Y)
thickness = chord*(y_max-y_min) #t
thickness_to_chord = thickness/chord #t/c
print(thickness_to_chord)
chord_x = [x_smallest,x_biggest]
chord_y = [Y[x_smallest_pos],Y[x_biggest_pos]]
plt.axvline(x = x_1,color = 'red')
plt.axvline(x = x_2,color = 'red')

plt.plot(X, Y)
plt.plot(chord_x,chord_y)
plt.show()