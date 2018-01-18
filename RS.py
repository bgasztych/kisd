import math

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
        self.t = math.floor(self.r / 2)
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
        msg_out = msg + remainder
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
                    raise ValueError("Bledy niekorygowalne")
                else:
                    msg_out = gf_shift_poly_right(msg_out, 1)
                    # print("Shift RIGHT: 1")
                    shifts += 1

    def calc_syndromes(self, msg):
        # Wyliczamy syndromy
        #
        synd = [0] * self.r
        for i in range(0, self.r):
            synd[i] = gf_poly_eval(msg, gf_pow(2, i))
        return [0] + synd

    def forney_syndromes(self, synd, pos, nmess):
        # Compute Forney syndromes, which computes a modified syndromes to compute only errors (erasures are trimmed out). Do not confuse this with Forney algorithm, which allows to correct the message based on the location of errors.
        erase_pos_reversed = [nmess - 1 - p for p in
                              pos]  # prepare the coefficient degree positions (instead of the erasures positions)

        # Optimized method, all operations are inlined
        fsynd = list(synd[1:])  # make a copy and trim the first coefficient which is always 0 by definition
        for i in range(0, len(pos)):
            x = gf_pow(2, erase_pos_reversed[i])
            for j in range(0, len(fsynd) - 1):
                fsynd[j] = gf_mul(fsynd[j], x) ^ fsynd[j + 1]

        # Equivalent, theoretical way of computing the modified Forney syndromes: fsynd = (erase_loc * synd) % x^(n-k)
        # See Shao, H. M., Truong, T. K., Deutsch, L. J., & Reed, I. S. (1986, April). A single chip VLSI Reed-Solomon decoder. In Acoustics, Speech, and Signal Processing, IEEE International Conference on ICASSP'86. (Vol. 11, pp. 2151-2154). IEEE.ISO 690
        # erase_loc = rs_find_errata_locator(erase_pos_reversed, generator=generator) # computing the erasures locator polynomial
        # fsynd = gf_poly_mul(erase_loc[::-1], synd[1:]) # then multiply with the syndrome to get the untrimmed forney syndrome
        # fsynd = fsynd[len(pos):] # then trim the first erase_pos coefficients which are useless. Seems to be not necessary, but this reduces the computation time later in BM (thus it's an optimization).

        return fsynd

    def find_error_locator(self, synd, erase_loc=None, erase_count=0):
        '''Find error/errata locator and evaluator polynomials with Berlekamp-Massey algorithm'''
        # The idea is that BM will iteratively estimate the error locator polynomial.
        # To do this, it will compute a Discrepancy term called Delta, which will tell us if the error locator polynomial needs an update or not
        # (hence why it's called discrepancy: it tells us when we are getting off board from the correct value).

        # Init the polynomials
        if erase_loc:  # if the erasure locator polynomial is supplied, we init with its value, so that we include erasures in the final locator polynomial
            err_loc = list(erase_loc)
            old_loc = list(erase_loc)
        else:
            err_loc = [
                1]  # This is the main variable we want to fill, also called Sigma in other notations or more formally the errors/errata locator polynomial.
            old_loc = [
                1]  # BM is an iterative algorithm, and we need the errata locator polynomial of the previous iteration in order to update other necessary variables.
        # L = 0 # update flag variable, not needed here because we use an alternative equivalent way of checking if update is needed (but using the flag could potentially be faster depending on if using length(list) is taking linear time in your language, here in Python it's constant so it's as fast.

        # Fix the syndrome shifting: when computing the syndrome, some implementations may prepend a 0 coefficient for the lowest degree term (the constant). This is a case of syndrome shifting, thus the syndrome will be bigger than the number of ecc symbols (I don't know what purpose serves this shifting). If that's the case, then we need to account for the syndrome shifting when we use the syndrome such as inside BM, by skipping those prepended coefficients.
        # Another way to detect the shifting is to detect the 0 coefficients: by definition, a syndrome does not contain any 0 coefficient (except if there are no errors/erasures, in this case they are all 0). This however doesn't work with the modified Forney syndrome, which set to 0 the coefficients corresponding to erasures, leaving only the coefficients corresponding to errors.
        synd_shift = 0
        if len(synd) > self.r:
            synd_shift = len(synd) - self.r

        for i in range(0,
                       self.r - erase_count):  # generally: nsym-erase_count == len(synd), except when you input a partial erase_loc and using the full syndrome instead of the Forney syndrome, in which case nsym-erase_count is more correct (len(synd) will fail badly with IndexError).
            if erase_loc:  # if an erasures locator polynomial was provided to init the errors locator polynomial, then we must skip the FIRST erase_count iterations (not the last iterations, this is very important!)
                K = erase_count + i + synd_shift
            else:  # if erasures locator is not provided, then either there's no erasures to account or we use the Forney syndromes, so we don't need to use erase_count nor erase_loc (the erasures have been trimmed out of the Forney syndromes).
                K = i + synd_shift

            # Compute the discrepancy Delta
            # Here is the close-to-the-books operation to compute the discrepancy Delta: it's a simple polynomial multiplication of error locator with the syndromes, and then we get the Kth element.
            # delta = gf_poly_mul(err_loc[::-1], synd)[K] # theoretically it should be gf_poly_add(synd[::-1], [1])[::-1] instead of just synd, but it seems it's not absolutely necessary to correctly decode.
            # But this can be optimized: since we only need the Kth element, we don't need to compute the polynomial multiplication for any other element but the Kth. Thus to optimize, we compute the polymul only at the item we need, skipping the rest (avoiding a nested loop, thus we are linear time instead of quadratic).
            # This optimization is actually described in several figures of the book "Algebraic codes for data transmission", Blahut, Richard E., 2003, Cambridge university press.
            delta = synd[K]
            for j in range(1, len(err_loc)):
                delta ^= gf_mul(err_loc[-(j + 1)], synd[
                    K - j])  # delta is also called discrepancy. Here we do a partial polynomial multiplication (ie, we compute the polynomial multiplication only for the term of degree K). Should be equivalent to brownanrs.polynomial.mul_at().
            # print "delta", K, delta, list(gf_poly_mul(err_loc[::-1], synd)) # debugline

            # Shift polynomials to compute the next degree
            old_loc = old_loc + [0]

            # Iteratively estimate the errata locator and evaluator polynomials
            if delta != 0:  # Update only if there's a discrepancy
                if len(old_loc) > len(
                        err_loc):  # Rule B (rule A is implicitly defined because rule A just says that we skip any modification for this iteration)
                    # if 2*L <= K+erase_count: # equivalent to len(old_loc) > len(err_loc), as long as L is correctly computed
                    # Computing errata locator polynomial Sigma
                    new_loc = gf_poly_scale(old_loc, delta)
                    old_loc = gf_poly_scale(err_loc, gf_inverse(
                        delta))  # effectively we are doing err_loc * 1/delta = err_loc // delta
                    err_loc = new_loc
                    # Update the update flag
                    # L = K - L # the update flag L is tricky: in Blahut's schema, it's mandatory to use `L = K - L - erase_count` (and indeed in a previous draft of this function, if you forgot to do `- erase_count` it would lead to correcting only 2*(errors+erasures) <= (n-k) instead of 2*errors+erasures <= (n-k)), but in this latest draft, this will lead to a wrong decoding in some cases where it should correctly decode! Thus you should try with and without `- erase_count` to update L on your own implementation and see which one works OK without producing wrong decoding failures.

                # Update with the discrepancy
                err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta))

        # Check if the result is correct, that there's not too many errors to correct
        while len(err_loc) and err_loc[0] == 0: del err_loc[
            0]  # drop leading 0s, else errs will not be of the correct size
        errs = len(err_loc) - 1
        # TODO Odkomentowac
        # if (errs - erase_count) * 2 + erase_count > self.correction_symbols:
        #     raise ValueError("Too many errors to correct")  # too many errors to correct

        return err_loc

    def find_errors(self, err_loc, nmess):  # nmess is len(msg_in)
        '''Find the roots (ie, where evaluation = zero) of error polynomial by brute-force trial, this is a sort of Chien's search
        (but less efficient, Chien's search is a way to evaluate the polynomial such that each evaluation only takes constant time).'''
        errs = len(err_loc) - 1
        err_pos = []
        for i in range(
                nmess):  # normally we should try all 2^8 possible values, but here we optimize to just check the interesting symbols
            if gf_poly_eval(err_loc,
                            gf_pow(2, i)) == 0:  # It's a 0? Bingo, it's a root of the error locator polynomial,
                # in other terms this is the location of an error
                err_pos.append(nmess - 1 - i)
        # Sanity check: the number of errors/errata positions found should be exactly the same as the length of the errata locator polynomial
        # if len(err_pos) != errs: TODO Odkomenrtować
        # couldn't find error locations
        # raise ValueError("Too many (or few) errors found by Chien Search for the errata locator polynomial!") TODO Odkomenrtować
        return err_pos

    def find_errata_locator(self, e_pos):
        '''Compute the erasures/errors/errata locator polynomial from the erasures/errors/errata positions
           (the positions must be relative to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to "h_ll_ worldxxxxxxxxx"
           with xxxxxxxxx being the ecc of length n-k=9, here the string positions are [1, 4], but the coefficients are reversed
           since the ecc characters are placed as the first coefficients of the polynomial, thus the coefficients of the
           erased characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified as an argument.'''

        e_loc = [
            1]  # just to init because we will multiply, so it must be 1 so that the multiplication starts correctly without nulling any term
        # erasures_loc = product(1 - x*alpha**i) for i in erasures_pos and where alpha is the alpha chosen to evaluate polynomials.
        for i in e_pos:
            e_loc = gf_poly_mul(e_loc, gf_poly_add([1], [gf_pow(2, i), 0]))
        return e_loc

    def find_error_evaluator(self, synd, err_loc, nsym):
        '''Compute the error (or erasures if you supply sigma=erasures locator polynomial, or errata) evaluator polynomial Omega
           from the syndrome and the error/erasures/errata locator Sigma.'''

        # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
        _, remainder = gf_poly_div(gf_poly_mul(synd, err_loc),
                                   ([1] + [0] * (nsym + 1)))  # first multiply syndromes * errata_locator, then do a
        # polynomial division to truncate the polynomial to the
        # required length

        # Faster way that is equivalent
        # remainder = gf_poly_mul(synd, err_loc) # first multiply the syndromes with the errata locator polynomial
        # remainder = remainder[len(remainder)-(nsym+1):] # then slice the list to truncate it (which represents the polynomial), which
        # is equivalent to dividing by a polynomial of the length we want

        return remainder

    def correct_errata(self, msg_in, synd, err_pos):  # err_pos is a list of the positions of the errors/erasures/errata
        '''Forney algorithm, computes the values (error magnitude) to correct the input message.'''
        # calculate errata locator polynomial to correct both errors and erasures (by combining the errors positions given by the error locator polynomial found by BM with the erasures positions given by caller)
        coef_pos = [len(msg_in) - 1 - p for p in
                    err_pos]  # need to convert the positions to coefficients degrees for the errata locator algo to work (eg: instead of [0, 1, 2] it will become [len(msg)-1, len(msg)-2, len(msg) -3])
        err_loc = self.find_errata_locator(coef_pos)
        # calculate errata evaluator polynomial (often called Omega or Gamma in academic papers)
        err_eval = self.find_error_evaluator(synd[::-1], err_loc, len(err_loc) - 1)[::-1]

        # Second part of Chien search to get the error location polynomial X from the error positions in err_pos (the roots of the error locator polynomial, ie, where it evaluates to 0)
        X = []  # will store the position of the errors
        for i in range(0, len(coef_pos)):
            l = 255 - coef_pos[i]
            X.append(gf_pow(2, -l))

        # Forney algorithm: compute the magnitudes
        E = [0] * (len(
            msg_in))  # will store the values that need to be corrected (substracted) to the message containing errors. This is sometimes called the error magnitude polynomial.
        Xlength = len(X)
        for i, Xi in enumerate(X):

            Xi_inv = gf_inverse(Xi)

            # Compute the formal derivative of the error locator polynomial (see Blahut, Algebraic codes for data transmission, pp 196-197).
            # the formal derivative of the errata locator is used as the denominator of the Forney Algorithm, which simply says that the ith error value is given by error_evaluator(gf_inverse(Xi)) / error_locator_derivative(gf_inverse(Xi)). See Blahut, Algebraic codes for data transmission, pp 196-197.
            err_loc_prime_tmp = []
            for j in range(0, Xlength):
                if j != i:
                    err_loc_prime_tmp.append(gf_sub(1, gf_mul(Xi_inv, X[j])))
            # compute the product, which is the denominator of the Forney algorithm (errata locator derivative)
            err_loc_prime = 1
            for coef in err_loc_prime_tmp:
                err_loc_prime = gf_mul(err_loc_prime, coef)
            # equivalent to: err_loc_prime = functools.reduce(gf_mul, err_loc_prime_tmp, 1)

            # Compute y (evaluation of the errata evaluator polynomial)
            # This is a more faithful translation of the theoretical equation contrary to the old forney method. Here it is an exact reproduction:
            # Yl = omega(Xl.inverse()) / prod(1 - Xj*Xl.inverse()) for j in len(X)
            y = gf_poly_eval(err_eval[::-1], Xi_inv)  # numerator of the Forney algorithm (errata evaluator evaluated)
            y = gf_mul(gf_pow(Xi, 1), y)

            # Compute the magnitude
            magnitude = gf_div(y,
                               err_loc_prime)  # magnitude value of the error, calculated by the Forney algorithm (an equation in fact): dividing the errata evaluator with the errata locator derivative gives us the errata magnitude (ie, value to repair) the ith symbol
            E[err_pos[i]] = magnitude  # store the magnitude for this error into the magnitude polynomial

        # Apply the correction of values to get our message corrected! (note that the ecc bytes also gets corrected!)
        # (this isn't the Forney algorithm, we just apply the result of decoding here)
        msg_in = gf_poly_add(msg_in,
                             E)  # equivalent to Ci = Ri - Ei where Ci is the correct message, Ri the received (senseword) message, and Ei the errata magnitudes (minus is replaced by XOR since it's equivalent in GF(2^p)). So in fact here we substract from the received message the errors magnitude, which logically corrects the value to what it should be.
        return msg_in

    def correct_msg(self, msg_in, erase_pos=None):
        '''Reed-Solomon main decoding function'''
        if len(msg_in) > 255:  # can't decode, message is too big
            raise ValueError("Message is too long (%i when max is 255)" % len(msg_in))

        msg_out = list(msg_in)  # copy of message
        # erasures: set them to null bytes for easier decoding (but this is not necessary, they will be corrected anyway, but debugging will be easier with null bytes because the error locator polynomial values will only depend on the errors locations, not their values)
        if erase_pos is None:
            erase_pos = []
        else:
            for e_pos in erase_pos:
                msg_out[e_pos] = 0
        # check if there are too many erasures to correct (beyond the Singleton bound)
        if len(erase_pos) > self.r: raise ValueError("Too many erasures to correct")
        # prepare the syndrome polynomial using only errors (ie: errors = characters that were either replaced by null byte
        # or changed to another character, but we don't know their positions)
        synd = self.calc_syndromes(msg_out)
        # check if there's any error/erasure in the input codeword. If not (all syndromes coefficients are 0), then just return the message as-is.
        if max(synd) == 0:
            return msg_out[:-self.r], msg_out[-self.r:]  # no errors

        # compute the Forney syndromes, which hide the erasures from the original syndrome (so that BM will just have to deal with errors, not erasures)
        fsynd = self.forney_syndromes(synd, erase_pos, len(msg_out))
        # compute the error locator polynomial using Berlekamp-Massey
        err_loc = self.find_error_locator(fsynd, erase_count=len(erase_pos))
        # locate the message errors using Chien search (or brute-force search)
        err_pos = self.find_errors(err_loc[::-1], len(msg_out))
        if err_pos is None:
            raise ValueError("Could not locate error")  # error location failed

        # Find errors values and apply them to correct the message
        # compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
        msg_out = self.correct_errata(msg_out, synd, (
                erase_pos + err_pos))  # note that we here use the original syndrome, not the forney syndrome
        # (because we will correct both errors and erasures, so we need the full syndrome)
        # check if the final message is fully repaired
        synd = self.calc_syndromes(msg_out)
        # if max(synd) > 0: TODO Odkomenrtować
        #     raise ValueError("Could not correct message")  # message could not be repaired TODO Odkomenrtować
        # return the successfully decoded message
        return msg_out[:-self.r], msg_out[-self.r:]  # also return the corrected ecc block so that the user can check()


def main():
    print("=====KOD RS=====\n")

    info = "One morning, when Gregor Samsa woke from troubled dreams, he found himself transformed in his bed into a horrible vermin. He lay on his armour-like back, and if he lifted his head a little he could see"
    info_in_unicode = [ord(x) for x in info]
    print("Info: %s LEN: %d" % (info, len(info)))
    print("Info ascii:      %s" % info_in_unicode)

    init_tables(0x11d, 2, 8)
    rs = RS()
    encoded = rs.encode(info_in_unicode)
    print("Encoded:         %s" % encoded)

    encoded_damaged = list(encoded)
    encoded_damaged[4] = 88
    print("Encoded damaged: %s" % encoded_damaged)
    # print(''.join([chr(x) for x in encoded]))
    decoded_simple = rs.decode_simple(encoded_damaged)
    print("Decoded:         %s" % decoded_simple)
    print("Info: %s" % ''.join([chr(x) for x in decoded_simple[:rs.k]]))
    print("Encoded == Decoded: %s" % (encoded == decoded_simple))


    # decoded, ecc = rs.correct_msg(encoded)
    # print(decoded)
    # print(''.join([chr(x) for x in decoded]))


if __name__ == "__main__":
    main()
