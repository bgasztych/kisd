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
    print([hex(i) for i in correct_primes])


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


def rs_generator_poly(nsym):
    '''Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)'''
    g = [1]
    for i in range(0, nsym):
        g = gf_poly_mul(g, [1, gf_pow(2, i)])
    return g


class RS:
    # GF(2^m)
    # m = 7 => Number of bits per symbol
    # n = 2^m - 1 = 127 => Codeword length
    # k = 5
    def __init__(self, n=127, k=5):
        self.n = n
        self.k = k
        self.generator = rs_generator_poly(self.n - self.k)

    def encode(self, msg):
        if len(msg) != self.k:
            raise ValueError("Message length must be equal to k")

        _, remainder = gf_poly_div(msg + [0] * (len(self.generator) - 1), self.generator)
        msg_out = msg + remainder
        # Return the codeword
        return msg_out


def main():
    print("=====KOD RS=====\n")

    s = "hello"
    b = [ord(x) for x in s]
    print(s)
    print(b)
    # print(b'%x' % s.encode("utf8"))
    # s_2 = [chr(x) for x in b]
    # print("".join(s_2))

    # print_prime_polys(find_prime_polys(2, 7))
    init_tables(0x83, 2, 7)
    rs = RS()
    print(rs.encode(b))


if __name__ == "__main__":
    main()
