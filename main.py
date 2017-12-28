import GF
import Utils


def main():
    print("=====KODOWANIE I SZYFROWANIE DANYCH=====\n")

    # gf1 = GF.GFSimple(1, 7)
    # gf2 = GF.GFSimple(2, 7)
    #
    # print("=====Galois Field Simple=====")
    # print("GF1: %s\nGF2: %s" % (gf1, gf2))
    # print("Sum: %s" % (gf1 + gf2))
    # print("Mul: %s" % (gf1 * gf2))
    # print("Invert gf1: %s" % (~gf1))
    # print("Invert gf2: %s" % (~gf2))
    # print("Generators amount: %d" % GF.GFSimple.get_generators_amount(7))
    # print("Generators: %s" % GF.GFSimple.get_generators(7))
    #
    # a = [3, 4, 1]
    # m = [4, 5, 7]
    # print("\n=====Chinese Remainder Theorem=====")
    # for i in range(0, a.__len__()):
    #     print("x ≡ %d ( mod %d)" % (a[i], m[i]))
    #
    # print("Chinese remainder theorem: %d" % Utils.Utils.crt(a, m))

    # gfe = GF.GFExtended(5, [1, 0, 1, 0, 1])
    # print(gfe)
    # elements = GF.GFExtended.get_all_elements(3)
    # for e in elements:
    #     print(e)
    # gfe1 = GF.GFExtended(3, [1, 0, 1])
    # gfe2 = GF.GFExtended(3, [1, 1, 1])
    # print(gfe1 * gfe2)
    gfe1 = GF.GFExtended(3, [1, 0, 0])
    print(gfe1)
    print(~gfe1)


if __name__ == "__main__":
    main()
