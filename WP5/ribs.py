from WP4.USE_THIS_loads_main import get_integrated
from WP5.wing import wing1


class Rib:
    def __init__(self, number):
        self.moment = get_integrated()[1]
        self.sigma = []
        self.positions = [0]
        #  self.number = len(self.positions)
        self.number = number

    def get_spacing(self):
        distance = wing1.length / self.number
        return distance

    def get_positions(self):
        spacing = self.get_spacing()
        for i in range(len(wing1.span)):
            if wing1.span[i] - self.positions[-1] >= spacing:
                self.positions.append(self.positions[-1] + spacing)
                if self.positions[-1] >= wing1.length:
                    self.positions[-1] = wing1.length
                    break
        return self.positions

    def get_distance(self, point):      # works for an integer might be an issue with passing the list in
        self.positions = self.get_positions()
        if point <= self.positions[-2]:
            distance = self.get_spacing()
            return distance
        else:
            return 22.445 - self.positions[-2]



'''
    def get_sigma(self):
        a = 0
        for i in wing1.span:
            self.sigma.append((self.moment[a] * (wingbox1.get_height(i) / 2)) / wingbox1.get_moment_of_inertia(i))
            a += 1
        return self.sigma

    def get_spacing(self):
        spacing = []
        sigma = [abs(i) for i in self.get_sigma()]
        for i in sigma:
            if i != 0:
                spacing.append((np.abs((wingbox1.stringer.clamp_factor * np.pi ** 2 * wingbox1.stringer.young_modulus
                                        * wingbox1.stringer.get_moment_of_inertia()) / (i * wingbox1.stringer.area))) ** 0.5)
            else:
                spacing.append(100)
        return spacing

    def get_placement(self):
        inter = self.get_spacing()
        while self.positions[-1] < 22.445:
            a = 0
            b = 0
            for i in wing1.span:
                if i > self.positions[-1] and b == 0:
                    self.positions.append(self.positions[-1] + inter[a])
                    b = 1
                a += 1
            if self.positions[-1] > 22.445:
                self.positions[-1] = 22.445
        return self.positions

    def get_distance(self, point):
        self.positions = self.get_placement()
        for i in range(1, len(self.positions)):
            if self.positions[i] >= point >= self.positions[i - 1]:
                return self.positions[i] - self.positions[i - 1]

'''


rib1 = Rib(25)



