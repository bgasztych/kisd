import Utils


class GFSimple:

    def __init__(self, value, p):
        if isinstance(value, int) and isinstance(p, int):
            if value >= 0 and p > 0 and value < p and Utils.Utils.is_prime(p):
                self.__value = value
                self.__p = p
            else:
                raise ValueError("value must be >=0 and p must be > 0 and value must <p and p must be prime")

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
            raise TypeError("argument must be GF type")

    def __mul__(self, other):
        if isinstance(other, GFSimple):
            if self.__p == other.p:
                value = (self.__value * other.value) % self.__p
                p = self.__p
                return GFSimple(value, p)
            else:
                raise ValueError("p must be equal")
        else:
            raise TypeError("argument must be GF type")

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
