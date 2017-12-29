import random
from functools import reduce
from operator import mul, mod

import GF


class Utils:

    @staticmethod
    def rabin_miller(num):
        # Returns True if num is a prime number.

        s = num - 1
        t = 0
        while s % 2 == 0:
            # keep halving s while it is even (and use t
            # to count how many times we halve s)
            s = s // 2
            t += 1

        for trials in range(5):  # try to falsify num's primality 5 times
            a = random.randrange(2, num - 1)
            v = pow(a, s, num)
            if v != 1:  # this test does not apply if v is 1.
                i = 0
                while v != (num - 1):
                    if i == t - 1:
                        return False
                    else:
                        i = i + 1
                        v = (v ** 2) % num
        return True

    # http://primos.mat.br/primeiros_10000_primos.txt
    @staticmethod
    def is_prime(num):
        # Return True if num is a prime number. This function does a quicker
        # prime number check before calling rabinMiller().

        if num < 2:
            return False  # 0, 1, and negative numbers are not prime

        # About 1/3 of the time we can quickly determine if num is not prime
        # by dividing by the first few dozen prime numbers. This is quicker
        # than rabinMiller(), but unlike rabinMiller() is not guaranteed to
        # prove that a number is prime.
        low_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                      101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
                      197,
                      199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
                      313,
                      317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433,
                      439,
                      443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569,
                      571,
                      577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,
                      691,
                      701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827,
                      829,
                      839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971,
                      977,
                      983, 991, 997]

        if num in low_primes:
            return True

        # See if any of the low prime numbers can divide num
        for prime in low_primes:
            if num % prime == 0:
                return False

        # If all else fails, call rabinMiller() to determine if num is a prime.
        return Utils.rabin_miller(num)

    @staticmethod
    def generate_large_prime(keysize=1024):
        # Return a random prime number of keysize bits in size.
        while True:
            num = random.randrange(2 ** (keysize - 1), 2 ** (keysize))
            if Utils.is_prime(num):
                return num

    @staticmethod
    def egcd(b, n):
        x0, x1, y0, y1 = 1, 0, 0, 1
        while n != 0:
            q, b, n = b // n, n, b % n
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
        return b, x0, y0

    # x ≡ 3 ( mod 4)
    # x ≡ 4 ( mod 5)
    # x ≡ 1 ( mod 7)
    #
    # m = [4, 5, 7]
    # a = [3, 4, 1]
    @staticmethod
    def crt(a, m):
        M = reduce(mul, m)  # the product of m elements
        m_i = [M / item for item in m]
        b = map(mod, m_i, m)
        g, k, l = map(Utils.egcd, b, m)
        g, k, l = zip(g, k, l)  # transpose g, k and l arrays
        t = map(mod, k, m)
        e = map(mul, m_i, t)

        x_sum = sum(map(mul, a, e))
        x = x_sum % M

        return x

    @staticmethod
    def euler_totient(p):
        amount = 0

        for k in range(1, p + 1):
            b, _, __ = Utils.egcd(p, k)
            if b == 1:
                amount += 1

        return amount

    @staticmethod
    def get_biggest_power(poly):
        for i in range(len(poly) - 1, -1, -1):
            if poly[i] == 1:
                return i
        return None

    @staticmethod
    def divide_with_remainder(dividend, divider):
        biggest_dividend_power = Utils.get_biggest_power(dividend)
        biggest_divider_power = Utils.get_biggest_power(divider)

        # Jeśli dzielnik mniejszy od dzielnej lub równy 0 nie trzeba dzielić
        if biggest_dividend_power < biggest_divider_power or (
                biggest_dividend_power >= 1 and biggest_divider_power == 0):
            return [[0], dividend]

        max_power = max(biggest_dividend_power, biggest_divider_power)
        result = [[0] * max_power] * 2
        new_dividend = list(dividend)
        while True:
            i = len(new_dividend) - 1
            while new_dividend[i] == 0 and i > 0:
                i -= 1
            if i < 0:
                return result

            # Dzielenie przez najwyższą potęgę dzielnej
            current_power = i - biggest_divider_power

            # Jeśli potęga mniejsza od zera to koniec dzielenia
            if current_power < 0:
                if i == 0:
                    return [dividend, [0] * max_power]
                else:
                    return result

            # Dzielenie największego wyrazu dzielnej przez największy wyraz dzielnika
            div_result = [0] * (current_power + 1)
            gfs1 = GF.GFSimple(new_dividend[i], 2)
            gfs2 = GF.GFSimple(divider[biggest_divider_power], 2)
            gfs_result = gfs1 * (~gfs2)
            div_result[current_power] = gfs_result.value

            # Dodanie wyniku dzielenia do result
            result[0] = Utils.add_polynomials(result[0], div_result)

            # Mnożenie dzielnika przez wynik dzielenia
            subtrahend = Utils.mul_polynomials(divider, div_result)

            # Odejmowanie wymnożonej wartości od dzielnej
            new_dividend = Utils.sub_polynomials(new_dividend, subtrahend)

            # Zwrócenie wyniku jeśli nie można dalej dzielić
            if Utils.get_biggest_power(new_dividend) < biggest_divider_power:
                result[1] = new_dividend
                return result

    @staticmethod
    def add_polynomials(poly1, poly2):
        if len(poly1) >= len(poly2):
            for i in range(0, len(poly2)):
                gfs1 = GF.GFSimple(poly1[i], 2)
                gfs2 = GF.GFSimple(poly2[i], 2)
                gfs_result = gfs1 + gfs2
                poly1[i] = gfs_result.value
            return poly1
        else:
            for i in range(0, len(poly1)):
                gfs1 = GF.GFSimple(poly1[i], 2)
                gfs2 = GF.GFSimple(poly2[i], 2)
                gfs_result = gfs1 + gfs2
                poly2[i] = gfs_result.value
            return poly2

    @staticmethod
    def sub_polynomials(poly1, poly2):
        if len(poly1) >= len(poly2):
            for i in range(0, len(poly2)):
                gfs1 = GF.GFSimple(poly1[i], 2)
                gfs2 = GF.GFSimple(poly2[i], 2)
                gfs_result = gfs1 + (~gfs2)
                poly1[i] = gfs_result.value
            return poly1
        else:
            for i in range(0, len(poly1)):
                gfs1 = GF.GFSimple(poly1[i], 2)
                gfs2 = GF.GFSimple(poly2[i], 2)
                gfs_result = gfs1 + (~gfs2)
                poly2[i] = gfs_result.value
            return poly2

    @staticmethod
    def mul_polynomials(poly1, poly2):
        max_power = Utils.get_biggest_power(poly1) + Utils.get_biggest_power(poly2)

        mul_result = [0] * (max_power + 1)
        for i in range(len(poly1) - 1, -1, -1):
            for j in range(len(poly2) - 1, -1, -1):
                gfs1 = GF.GFSimple(poly1[i], 2)
                gfs2 = GF.GFSimple(poly2[j], 2)
                gfs_mul_result = gfs1 * gfs2
                if gfs_mul_result.value > 0:
                    power = i + j
                    gfs3 = GF.GFSimple(mul_result[power], 2)
                    gfs_add_result = gfs3 + gfs_mul_result
                    mul_result[power] = gfs_add_result.value
        return mul_result
