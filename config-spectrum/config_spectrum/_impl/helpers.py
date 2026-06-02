




def ratio(fp):
    """
    Parse a Python `float`-type floating point expression
    to produce an initializer expression for a std::ratio.
    There are occasions where an
    unrecognizable std::ratio is generated, due to base-ten roundoff error.
    For example, 2000.003 might become 20000030000000004/100000000000000
    instead of 2000003/1000.
    :param fp: scalar
    :return: string
    """
    from math import modf as math_modf
    denom = 1
    numer = abs(fp)
    sgn = 1 if fp >= 0 else -1
    fp, ip = math_modf(numer)
    while fp != 0:
        denom *= 10
        numer *= 10
        fp, ip = math_modf(numer)
    return f"std::ratio<{int(sgn*numer)},{int(denom)}>"



