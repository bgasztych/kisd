import Utils


class GFSimple:

    def __init__(self, value, p):
        if isinstance(value, int) and isinstance(p, int):
            if value >= 0 and p > 0 and value < p and Utils.Utils.is_prime(p):
                self.__value = value
                self.__p = p
            else:
                raise ValueError("value must be >=0 and p must be > 0 and value must < p and p must be prime")

        else:
            raise TypeError("value and p must be integers")

    @property
    def p(self):
        return self.__p

    @property
    def value(self):
        return self.__value

    def __add__(self, other):
        if isinstance(other, GFSimple):
            if self.__p == other.p:
                value = (self.__value + other.value) % self.__p
                p = self.__p
                return GFSimple(value, p)
            else:
                raise ValueError("p must be equal")
        else:
            raise TypeError("argument must be GFSimple type")

    def __mul__(self, other):
        if isinstance(other, GFSimple):
            if self.__p == other.p:
                value = (self.__value * other.value) % self.__p
                p = self.__p
                return GFSimple(value, p)
            else:
                raise ValueError("p must be equal")
        else:
            raise TypeError("argument must be GFSimple type")

    def __invert__(self):
        g, x, _ = Utils.Utils.egcd(self.__value, self.__p)
        if g == 1:
            value = x % self.__p
            p = self.__p
            return GFSimple(value, p)
        else:
            return None

    def __str__(self):
        return "Value= %d p= %d" % (self.__value, self.__p)

    @staticmethod
    def get_generators(p):
        res = []
        exp = [None] * (p - 1)

        for i in range(2, p - 1):
            start = 1
            flag = True

            j = 0
            for j in range(0, (p - 1) // 2):
                start = (start * i) % p
                exp[j] = start
                if (start % p) == 1:
                    flag = False
                    break
            if flag:
                for k in range(j + 1, (p - 1)):
                    start = (start * i) % p
                    exp[k] = start
                # print("exp: %s" % exp)
                for m in range(1, p - 1):
                    b, _, __ = Utils.Utils.egcd(m, p - 1)
                    if b == 1:
                        res.append(exp[m - 1])
                return sorted(res)

    @staticmethod
    def get_generators_amount(p):
        return Utils.Utils.euler_totient(p - 1)


class GFExtended:
    generators = [[], [], [1, 1, 1], [1, 1, 0, 1], [1, 1, 0, 0, 1], [1, 0, 1, 0, 0, 1]]

    def __init__(self, m, element):
        if m == len(element):
            self.m = m
            self.element = element
        else:
            raise ValueError("m must be equal element size")

    def __str__(self):
        out = ""
        for i in range(len(self.element) - 1, -1, -1):
            if self.element[i] != 0:
                value_str = "x^" + str(i)
                if i == 0:
                    value_str = "1"
                elif i == 1:
                    value_str = "x"
                if len(out) == 0:
                    out = value_str
                else:
                    out += " + " + value_str
        return out

    def __add__(self, other):
        if isinstance(other, GFExtended):
            if self.m == other.m:
                result = []
                for i in range(0, len(self.element)):
                    gfs1 = GFSimple(self.element[i], 2)
                    gfs2 = GFSimple(other.element[i], 2)
                    gfs_result = gfs1 + gfs2
                    result.append(gfs_result.value)
                return GFExtended(self.m, result)
            else:
                raise ValueError("m must be equal")
        else:
            raise TypeError("argument must be GFExtended type")

    def __mul__(self, other):
        if isinstance(other, GFExtended):
            if self.m == other.m:
                result = [0] * (len(self.element) * 2 + 1)
                for i in range(len(self.element) - 1, -1, -1):
                    for j in range(len(other.element) - 1, -1, -1):
                        gfs1 = GFSimple(self.element[i], 2)
                        gfs2 = GFSimple(other.element[j], 2)
                        gfs_result = gfs1 * gfs2
                        if gfs_result.value > 0:
                            power = i + j
                            gfs2 = GFSimple(result[power], 2)
                            gfs_result2 = gfs2 + gfs_result
                            result[power] = gfs_result2.value
                return result
            else:
                raise ValueError("m must be equal")
        else:
            raise TypeError("argument must be GFExtended type")

    def __invert__(self):
        result = []
        for i in range(len(self.element) - 1, -1, -1):
            result.append(self.element[i])
        return GFExtended(self.m, result)

    @staticmethod
    def get_biggest_power(poly):
        for i in range(len(poly) - 1, -1, -1):
            if poly[i] == 1:
                return i
        return None

    @staticmethod
    def generate_sequence(generator):
        gfe_size = 2 ** (len(generator) - 1)
        power = GFExtended.get_biggest_power(generator)
        s_array = []

        # Wypełnienie tablicy s indeksami poteg z generatora
        for i in range(power):
            if generator[i] != 0:
                s_array.append(i)

        # Wypełnienie sekwencji ciągiem początkowym
        sequence = [None] * ((gfe_size - 1) * power)
        sequence[0] = 1
        for i in range(1, power):
            sequence[i] = 0

        # Dodawanie elementów sekwencji zgodnie ze wzorem rekurencyjnym
        for i in range(power, len(sequence)):
            gfs = GFSimple(0, 2)
            for j in range(len(s_array)):
                gfs += GFSimple(sequence[i - power + s_array[j]], 2)
            sequence[i] = gfs.value
        return sequence

    @staticmethod
    def get_all_elements(m):
        out = []
        generator = GFExtended.generators[m]
        sequence = GFExtended.generate_sequence(generator)

        # Tworzenie z sekwencji kolejnych elementów ciała
        for i in range(0, len(sequence), m):
            out.append(GFExtended(m, sequence[i:i + m]))
        return out
