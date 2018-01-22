import math
import random

gf_exp = [0] * 512  # Create list of 512 elements. In Python 2.6+, consider using bytearray
gf_log = [0] * 256
field_size = int(2 ** 8 - 1)


def rwh_primes1(n):
    # http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    ''' Returns  a list of primes < n '''
    sieve = [True] * (n / 2)
    for i in range(3, int(n ** 0.5) + 1, 2):
        if sieve[i // 2]:
            sieve[i * i / 2::i] = [False] * ((n - i * i - 1) / (2 * i) + 1)
    return [2] + [2 * i + 1 for i in range(1, n / 2) if sieve[i]]


def find_prime_polys(generator=2, c_exp=8, fast_primes=False, single=False):
    '''Compute the list of prime polynomials for the given generator and galois field characteristic exponent.'''
    # fast_primes will output less results but will be significantly faster.
    # single will output the first prime polynomial found, so if all you want is to just find one prime polynomial to generate the LUT for Reed-Solomon to work, then just use that.

    # A prime polynomial (necessarily irreducible) is necessary to reduce the multiplications in the Galois Field, so as to avoid overflows.
    # Why do we need a "prime polynomial"? Can't we just reduce modulo 255 (for GF(2^8) for example)? Because we need the values to be unique.
    # For example: if the generator (alpha) = 2 and c_exp = 8 (GF(2^8) == GF(256)), then the generated Galois Field (0, 1, α, α^1, α^2, ..., α^(p-1)) will be galois field it becomes 0, 1, 2, 4, 8, 16, etc. However, upon reaching 128, the next value will be doubled (ie, next power of 2), which will give 256. Then we must reduce, because we have overflowed above the maximum value of 255. But, if we modulo 255, this will generate 256 == 1. Then 2, 4, 8, 16, etc. giving us a repeating pattern of numbers. This is very bad, as it's then not anymore a bijection (ie, a non-zero value doesn't have a unique index). That's why we can't just modulo 255, but we need another number above 255, which is called the prime polynomial.
    # Why so much hassle? Because we are using precomputed look-up tables for multiplication: instead of multiplying a*b, we precompute alpha^a, alpha^b and alpha^(a+b), so that we can just use our lookup table at alpha^(a+b) and get our result. But just like in our original field we had 0,1,2,...,p-1 distinct unique values, in our "LUT" field using alpha we must have unique distinct values (we don't care that they are different from the original field as long as they are unique and distinct). That's why we need to avoid duplicated values, and to avoid duplicated values we need to use a prime irreducible polynomial.

    # Here is implemented a bruteforce approach to find all these prime polynomials, by generating every possible prime polynomials (ie, every integers between field_charac+1 and field_charac*2), and then we build the whole Galois Field, and we reject the candidate prime polynomial if it duplicates even one value or if it generates a value above field_charac (ie, cause an overflow).
    # Note that this algorithm is slow if the field is too big (above 12), because it's an exhaustive search algorithm. There are probabilistic approaches, and almost surely prime approaches, but there is no determistic polynomial time algorithm to find irreducible monic polynomials. More info can be found at: http://people.mpi-inf.mpg.de/~csaha/lectures/lec9.pdf
    # Another faster algorithm may be found at Adleman, Leonard M., and Hendrik W. Lenstra. "Finding irreducible polynomials over finite fields." Proceedings of the eighteenth annual ACM symposium on Theory of computing. ACM, 1986.

    # Prepare the finite field characteristic (2^p - 1), this also represent the maximum possible value in this field
    root_charac = 2  # we're in GF(2)
    field_charac = int(root_charac ** c_exp - 1)
    field_charac_next = int(root_charac ** (c_exp + 1) - 1)

    prim_candidates = []
    if fast_primes:
        prim_candidates = rwh_primes1(
            field_charac_next)  # generate maybe prime polynomials and check later if they really are irreducible
        prim_candidates = [x for x in prim_candidates if x > field_charac]  # filter out too small primes
    else:
        prim_candidates = range(field_charac + 2, field_charac_next,
                                root_charac)  # try each possible prime polynomial, but skip even numbers (because divisible by 2 so necessarily not irreducible)

    # Start of the main loop
    correct_primes = []
    for prim in prim_candidates:  # try potential candidates primitive irreducible polys
        seen = bytearray(
            field_charac + 1)  # memory variable to indicate if a value was already generated in the field (value at index x is set to 1) or not (set to 0 by default)
        conflict = False  # flag to know if there was at least one conflict

        # Second loop, build the whole Galois Field
        x = 1
        for i in range(field_charac):
            # Compute the next value in the field (ie, the next power of alpha/generator)
            x = gf_mult_noLUT(x, generator, prim, field_charac + 1)

            # Rejection criterion: if the value overflowed (above field_charac) or is a duplicate of a previously generated power of alpha, then we reject this polynomial (not prime)
            if x > field_charac or seen[x] == 1:
                conflict = True
                break
            # Else we flag this value as seen (to maybe detect future duplicates), and we continue onto the next power of alpha
            else:
                seen[x] = 1

        # End of the second loop: if there's no conflict (no overflow nor duplicated value), this is a prime polynomial!
        if not conflict:
            correct_primes.append(prim)
            if single: return prim

    # Return the list of all prime polynomials
    return correct_primes  # you can use the following to print the hexadecimal representation of each prime polynomial: print [hex(i) for i in correct_primes]


def print_prime_polys(correct_primes):
    print([bin(i) for i in correct_primes])


def gf_mult_noLUT(x, y, prim=0, field_charac_full=256, carryless=True):
    '''Galois Field integer multiplication using Russian Peasant Multiplication algorithm (faster than the standard multiplication + modular reduction).
    If prim is 0 and carryless=False, then the function produces the result for a standard integers multiplication (no carry-less arithmetics nor modular reduction).'''
    r = 0
    while y:  # while y is above 0
        if y & 1: r = r ^ x if carryless else r + x  # y is odd, then add the corresponding x to r (the sum of all x's corresponding to odd y's will give the final product). Note that since we're in GF(2), the addition is in fact an XOR (very important because in GF(2) the multiplication and additions are carry-less, thus it changes the result!).
        y = y >> 1  # equivalent to y // 2
        x = x << 1  # equivalent to x*2
        if prim > 0 and x & field_charac_full: x = x ^ prim  # GF modulo: if x >= 256 then apply modular reduction using the primitive polynomial (we just subtract, but since the primitive number can be above 256 then we directly XOR).

    return r


def init_tables(prim=0x11d, generator=2, c_exp=8):
    '''Precompute the logarithm and anti-log tables for faster computation later, using the provided primitive polynomial.'''
    # prim is the primitive (binary) polynomial. Since it's a polynomial in the binary sense,
    # it's only in fact a single galois field value between 0 and 255, and not a list of gf values.
    global gf_exp, gf_log, field_size
    field_size = int(2 ** c_exp - 1)
    gf_exp = [0] * (field_size * 2)  # anti-log (exponential) table
    gf_log = [0] * (field_size + 1)  # log table
    # For each possible value in the galois field 2^8, we will pre-compute the logarithm and anti-logarithm (exponential) of this value
    x = 1
    for i in range(0, field_size):
        gf_exp[i] = x  # compute anti-log for this value and store it in a table
        gf_log[x] = i  # compute log at the same time
        x = gf_mult_noLUT(x, generator, prim, field_size + 1)

        # If you use only generator==2 or a power of 2, you can use the following which is faster than gf_mult_noLUT():
        # x <<= 1 # multiply by 2 (change 1 by another number y to multiply by a power of 2^y)
        # if x & 0x100: # similar to x >= 256, but a lot faster (because 0x100 == 256)
        # x ^= prim # substract the primary polynomial to the current value (instead of 255, so that we get a unique set made of coprime numbers), this is the core of the tables generation

    # Optimization: double the size of the anti-log table so that we don't need to mod 255 to
    # stay inside the bounds (because we will mainly use this table for the multiplication of two GF numbers, no more).
    for i in range(field_size, field_size * 2):
        gf_exp[i] = gf_exp[i - field_size]
    print(gf_exp)
    return [gf_log, gf_exp]


def gf_add(x, y):
    return x ^ y


def gf_sub(x, y):
    return x ^ y  # in binary galois field, subtraction is just the same as addition (since we mod 2)


def gf_poly_add(p, q):
    r = [0] * max(len(p), len(q))
    for i in range(0, len(p)):
        r[i + len(r) - len(p)] = p[i]
    for i in range(0, len(q)):
        r[i + len(r) - len(q)] ^= q[i]
    return r


def gf_div(x, y):
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return gf_exp[(gf_log[x] + 255 - gf_log[y]) % 255]


def gf_mul(x, y):
    if x == 0 or y == 0:
        return 0
    return gf_exp[
        (gf_log[x] + gf_log[y]) % field_size]  # should be gf_exp[(gf_log[x]+gf_log[y])%255] if gf_exp wasn't oversized


def gf_poly_mul(p, q):
    '''Multiply two polynomials, inside Galois Field'''
    # Pre-allocate the result array
    r = [0] * (len(p) + len(q) - 1)
    # Compute the polynomial multiplication (just like the outer product of two vectors,
    # we multiply each coefficients of p with all coefficients of q)
    for j in range(0, len(q)):
        for i in range(0, len(p)):
            r[i + j] ^= gf_mul(p[i], q[j])  # equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j]))
            # -- you can see it's your usual polynomial multiplication
    return r


def gf_poly_div(dividend, divisor):
    '''Fast polynomial division by using Extended Synthetic Division and optimized for GF(2^p) computations
    (doesn't work with standard polynomials outside of this galois field, see the Wikipedia article for generic algorithm).'''
    # CAUTION: this function expects polynomials to follow the opposite convention at decoding:
    # the terms must go from the biggest to lowest degree (while most other functions here expect
    # a list from lowest to biggest degree). eg: 1 + 2x + 5x^2 = [5, 2, 1], NOT [1, 2, 5]

    msg_out = list(dividend)  # Copy the dividend
    # normalizer = divisor[0] # precomputing for performance
    for i in range(0, len(dividend) - (len(divisor) - 1)):
        # msg_out[i] /= normalizer # for general polynomial division (when polynomials are non-monic), the usual way of using
        # synthetic division is to divide the divisor g(x) with its leading coefficient, but not needed here.
        coef = msg_out[i]  # precaching
        if coef != 0:  # log(0) is undefined, so we need to avoid that case explicitly (and it's also a good optimization).
            for j in range(1, len(
                    divisor)):  # in synthetic division, we always skip the first coefficient of the divisior,
                # because it's only used to normalize the dividend coefficient
                if divisor[j] != 0:  # log(0) is undefined
                    msg_out[i + j] ^= gf_mul(divisor[j], coef)  # equivalent to the more mathematically correct
                    # (but xoring directly is faster): msg_out[i + j] += -divisor[j] * coef

    # The resulting msg_out contains both the quotient and the remainder, the remainder being the size of the divisor
    # (the remainder has necessarily the same degree as the divisor -- not length but degree == length-1 -- since it's
    # what we couldn't divide from the dividend), so we compute the index where this separation is, and return the quotient and remainder.
    separator = -(len(divisor) - 1)
    return msg_out[:separator], msg_out[separator:]  # return quotient, remainder.


def gf_pow(x, power):
    return gf_exp[(gf_log[x] * power) % field_size]


def gf_inverse(x):
    return gf_exp[255 - gf_log[x]]  # gf_inverse(x) == gf_div(1, x)


def gf_poly_eval(poly, x):
    '''Evaluates a polynomial in GF(2^p) given the value for x. This is based on Horner's scheme for maximum efficiency.'''
    y = poly[0]
    for i in range(1, len(poly)):
        y = gf_mul(y, x) ^ poly[i]
    return y


def gf_poly_scale(p, x):
    r = [0] * len(p)
    for i in range(0, len(p)):
        r[i] = gf_mul(p[i], x)
    return r


def gf_poly_weight(poly):
    # Obliczenie wagi Hamminga (liczba niezerowych wspolrzednych wielomianu)
    weight = 0
    for i in range(len(poly)):
        if poly[i] != 0:
            weight += 1
    return weight


def gf_shift_poly_right(poly, shifts):
    r = [0] * len(poly)
    shifts %= len(poly)
    for i in range(len(r)):
        position = i + shifts
        if position >= len(r):
            position -= len(poly)
        r[position] = poly[i]
    return r


def gf_shift_poly_left(poly, shifts):
    r = [0] * len(poly)
    shifts %= len(poly)
    for i in range(len(r)):
        position = i - shifts
        if position < 0:
            position += len(poly)
        r[position] = poly[i]
    return r


def gf_solve(A, B):
    n = len(B)
    for p in range(n):
        maximum = p
        for i in range(p + 1, n):
            # print("SOLV: %s | %s" % (A[i][p][0], A[maximum][p][0]))

            if not isinstance(A[i][p], int):
                arg1 = (A[i][p])[0]
            else:
                arg1 = A[i][p]

            if not isinstance(A[maximum][p], int):
                arg2 = (A[maximum][p])[0]
            else:
                arg2 = A[maximum][p]
            # print("ARGS: %s | %s" % (arg1, arg2))

            if abs(arg1) > abs(arg2):
                maximum = i
        temp = A[p]
        A[p] = A[maximum]
        A[maximum] = temp
        t = B[p]
        B[p] = B[maximum]
        B[maximum] = t

        if not isinstance(A[p][p], int):
            arg3 = (A[p][p])[0]
        else:
            arg3 = (A[p][p])

        if abs(arg3) == 0:
            return None

        for i in range(p + 1, n):

            if not isinstance(A[i][p], int):
                arg1 = (A[i][p])[0]
            else:
                arg1 = A[i][p]

            if not isinstance(A[p][p], int):
                arg2 = (A[p][p])[0]
            else:
                arg2 = A[p][p]

            alpha = gf_div(arg1, arg2)

            # print("ALPHA: %s" % B[p][0])
            if not isinstance(B[i], int):
                arg1 = B[i][0]
            else:
                arg1 = B[i]

            if not isinstance(B[p], int):
                arg2 = B[p][0]
            else:
                arg2 = B[p]

            B[i] = gf_add(arg1, gf_mul(alpha, arg2))
            for j in range(p, n):

                if not isinstance(A[i][j], int):
                    arg1 = (A[i][j])[0]
                else:
                    arg1 = A[i][j]

                if not isinstance(A[p][j], int):
                    arg2 = (A[p][j])[0]
                else:
                    arg2 = A[p][j]

                A[i][j] = gf_add(arg1, gf_mul(alpha, arg2))

    x = [0] * n
    for i in range(n - 1, -1, -1):
        suma = 0
        for j in range(i + 1, n):
            suma = gf_add(suma, gf_mul(A[i][j], x[j]))
            x[i] = gf_div(gf_add(B[i], suma), A[i][i])
    return x


def gf_find_roots(error_locator_polynomial):
    v = (len(error_locator_polynomial) - 1)
    result = [0] * v
    last_found_root_index = 0
    for i in range(field_size, 0, -1):
        if gf_poly_eval(error_locator_polynomial, gf_exp[i]) == 0:
            result[last_found_root_index] = ((field_size-i+1) % field_size)
            last_found_root_index += 1
    return result


class RS:
    # GF(2^m)
    # m = 8 => Liczba bitów na symbol
    # n = 2^m - 1 = 255 => Długość wiadomości
    # k = 201 => Długość informacji
    # r = n - k = 255 - 201 = 54 => Liczba symboli kontrolnych
    # floor((n - k) / 2) = 27 => Maksymalna liczba symboli, która może być skorygowana
    def __init__(self, n=255, k=201):
        self.n = n
        self.k = k
        self.r = self.n - self.k
        self.t = int(math.floor(self.r / 2))
        self.generator = self.generator_poly()

    def generator_poly(self):
        # Wyznaczenie wielomianu generujacego
        # g(x) = (x+a)(x+a^2)...(x+a^r)
        g = [1]
        for i in range(0, self.r):
            g = gf_poly_mul(g, [1, gf_pow(2, i)])
        # print("generator_rs: %s\nlen: %d" % (g,  len(g)))
        return g

    def encode(self, msg):
        if len(msg) != self.k:
            raise ValueError("Message length must be equal to k")
        # Dzielimy wielomian informacyjny przez generator
        _, remainder = gf_poly_div(msg + [0] * self.r, self.generator)
        # Doklejamy czesc kontrolna
        msg_out = remainder + msg
        return msg_out

    def decode_simple(self, msg):
        msg_out = list(msg)

        # Wyzerowanie licznika przesuniec
        shifts = 0
        while True:
            # Obliczenie syndromu i jego wagi
            _, synd = gf_poly_div(msg_out, self.generator)
            weight = gf_poly_weight(synd)
            # print("Weight: %d Shifts: %d" % (weight, shifts))

            # Sprawdzenie warunku korekcji
            if weight <= self.t:
                # Waga syndromu jest mniejsza lub rowna zdolnosci korekcyjnej, wiec
                # bledy znajduja sie w czesci kontrolnej
                msg_out = gf_poly_add(msg_out, synd)
                # Przesuwamy wektor kodowy w lewo tyle razy ile zostal on przesuniety w prawo
                # aby odtworzyc jego pierwotna postac
                msg_out = gf_shift_poly_left(msg_out, shifts)
                # print("Shift LEFT: %d" % shifts)
                return msg_out
            else:
                # Waga syndromu wieksza od zdolnosci korekcyjnej, wiec
                # bledy znajduja sie w czesci informacyjnej

                # Jesli nastapilo k przesuniec i nie udalo sie skorygowac wektora kodowego
                # to wystapily bledy niekorygowalne
                if shifts == self.k:
                    msg_out = gf_shift_poly_left(msg_out, shifts)
                    return msg_out
                    # raise ValueError("Bledy niekorygowalne")
                else:
                    msg_out = gf_shift_poly_right(msg_out, 1)
                    # print("Shift RIGHT: 1")
                    shifts += 1

    def calc_syndromes(self, msg):
        # Wyliczamy syndrom
        synd = [0] * self.r
        for i in range(0, self.r):
            synd[i] = gf_poly_eval(msg, gf_pow(2, i))
        return synd

    def compose_s_matrix(self, syndromes, v):
        result = [[0 for x in range(v)] for y in range(v)]
        for i in range(v):
            for j in range(v):
                result[i][j] = syndromes[i + v - j - 1]
        return result

    def compose_x_matrix(self, error_locator_roots, v):
        print("COMPOSE: %s" % error_locator_roots)
        result = [[0 for x in range(v)] for y in range(v)]
        for i in range(v):
            alpha = gf_exp[i]
            for j in range(v):
                val = 1
                for a in range(1, error_locator_roots[j] + 1):
                    val = gf_mul(val, alpha)
                result[i][j] = val
        return result

    def compose_error_locator_eq_lhs(self, syndromes, v):
        result = [0] * v
        for i in range(v):
            result[i] = syndromes[v + i]
        return result

    def decode(self, msg):
        if len(msg) != self.n:
            raise ValueError("Message length must be = %d" % self.n)

        # syndromes = [0] * (2 * self.t)
        # for i in range(2 * self.t):
        #     _, synd = gf_poly_div(list(reversed(msg)), list(reversed([gf_exp[i], 1])))
        #     print("SYND: %s" % synd)
        #     syndromes[i] = synd[0]
        #     # print("GFEXP: %s" % [gf_exp[i], 1])
        syndromes = self.calc_syndromes(msg)
        print("SYNDROM: %s" % syndromes)
        syndromes = [221, 11, 48, 148, 182, 103, 154, 205, 132, 2, 247, 158, 61, 124, 232, 200, 1, 148, 204, 88, 40, 95, 67, 207, 117, 159, 80, 131, 166, 129, 217, 87, 87, 64, 43, 55, 112, 151, 138, 165, 228, 20, 9, 161, 95, 37, 206, 163, 79, 38, 25, 46, 0, 165]

        v = self.t
        error_locator_polynomial = None
        while error_locator_polynomial is None and v > 1:
            error_locator_polynomial = gf_solve(self.compose_s_matrix(syndromes, v), self.compose_error_locator_eq_lhs(syndromes, v))
            v -= 1

        error_locator_polynomial = [217, 21, 41, 100, 254, 245, 139]
        v = 7

        if error_locator_polynomial is None:
            return msg
        else:
            print("len erp: %d v: %d" % (len(error_locator_polynomial), v))
            temp = list(error_locator_polynomial)
            error_locator_polynomial = [0] * (v + 1)
            error_locator_polynomial[0] = 1
            print("len erp: %d v: %d len temp: %d" % (len(error_locator_polynomial), v, len(temp)))
            for i in range(len(temp)):
                print("I: %d" % i)
                error_locator_polynomial[i + 1] = temp[i]

            error_locator_roots = gf_find_roots(error_locator_polynomial)

            error_values = gf_solve(self.compose_x_matrix(error_locator_roots, v), syndromes[:v])

            if error_values is None:
                return msg

            msg_out = list(msg)
            for i in range(len(error_locator_roots)):
                location = error_locator_roots[i]
                msg_out[location] = gf_add(msg_out[location], error_values[i])

            return msg_out


def get_errors_percent_corrected(encoded, decoded, damaged_symbols_number):
    counter = 0
    for i in range(len(encoded)):
        if encoded[i] != decoded[i]:
            counter += 1
    percent = (counter / damaged_symbols_number) * 100
    return percent


def generate_errors(msg, t, magnitude, errors_number=0):
    if magnitude > t:
        raise ValueError("Magnitude must be <= t")
    msg_damaged = list(msg)
    if errors_number <= 0:
        errors_number = int(math.floor(t / magnitude))
    elif errors_number > int(math.floor(t / magnitude)):
        raise ValueError("Errors number illegal value")
    print("Errors number: %d" % errors_number)

    max_spacing = int((255 - magnitude) / errors_number)
    min_spacing = int(max_spacing * 0.5)
    print("Max spacing: %d" % max_spacing)
    start_position = int(random.randint(0, int(min_spacing * 0.25)))
    for i in range(errors_number):
        for j in range(magnitude):
            msg_damaged[i + j + start_position] = random.randint(0, 255)
        start_position += magnitude + random.randint(random.randint(1, min_spacing), max_spacing)
    # print(msg)
    # print(msg_damaged)

    positions = []
    for i in range(len(msg)):
        if msg[i] != msg_damaged[i]:
            positions.append(i)
    print("Errors positions: %s Count: %d" % (positions, len(positions)))
    return msg_damaged, len(positions)


# print("\033[0;31;48m" + str(e) + "\x1b[0m")
def main():
    # print("=====KOD RS=====\n")

    info = "One morning, when Gregor Samsa woke from troubled dreams, he found himself transformed in his bed into a horrible vermin. He lay on his armour-like back, and if he lifted his head a little he could see"
    info_in_unicode = [ord(x) for x in info]
    print("Info: %s LEN: %d" % (info, len(info)))
    print("Info ascii:      %s" % info_in_unicode)

    init_tables(0x11d, 2, 8)
    rs = RS()
    encoded = rs.encode(info_in_unicode)
    print("Encoded:         %s" % encoded)

    encoded_damaged = list(encoded)
    encoded_damaged[0] = 88
    encoded_damaged[2] = 88
    encoded_damaged[3] = 88
    encoded_damaged[10] = 88
    encoded_damaged[47] = 88
    encoded_damaged[205] = 88
    encoded_damaged[225] = 88
    damaged_symbols_number = 7

    encoded_damaged, damaged_symbols_number = generate_errors(encoded, rs.t, 10, 1)
    print("Encoded damaged: %s" % encoded_damaged)
    print("Info: %s" % ''.join([chr(x) for x in encoded_damaged[rs.r:]]))
    decoded_simple = rs.decode_simple(encoded_damaged)
    print("Decoded simple:  %s" % decoded_simple)
    print("Info: %s" % ''.join([chr(x) for x in decoded_simple[rs.r:]]))
    print("Encoded == Decoded SIMPLE: %s" % (encoded == decoded_simple))
    print("Errors percent corrected [SIMPLE]: %s%%\n" % get_errors_percent_corrected(encoded_damaged, decoded_simple, damaged_symbols_number))

    decoded = rs.decode(encoded_damaged)
    print("Decoded:         %s" % decoded)
    print("Info: %s" % ''.join([chr(x) for x in decoded[rs.r:]]))
    print("Encoded == Decoded: %s" % (encoded == decoded))
    print("Errors percent corrected [EXTENDED]: %s%%" % get_errors_percent_corrected(encoded_damaged, decoded, damaged_symbols_number))



if __name__ == "__main__":
    main()
