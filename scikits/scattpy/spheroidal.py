from numpy import *
from scipy.special import obl_cv, pro_cv, lpmv
from scipy.misc.common import factorial

import f_spheroid


def get_norm_cv(m, n, c, cv, type):
    eps = 1E-14
    value = 1
    sum = 0
    k = 0
    d = f_spheroid.sdmn(m, n, c, cv, type)

    isEven = (n - m) % 2
    if(isEven == 1):
        k = 1

    while (value > eps):
        if (type == 1):
            value = factorial(k + 2 * m) * d[k/2] * d[k/2] / ((2 * k + 2 * m + 1) * factorial(k))
        else:
            value = factorial(k + 2 * m) * d[k/2] * d[k/2] / ((2 * k + 2 * m + 1) * factorial(k))
        sum += value
        k += 2

    sum = sqrt(2 * sum)
    return sum


def get_pro_norm_cv(m, n, c, cv):
    return get_norm_cv(m, n, c, cv, 1)


def get_obl_norm_cv(m, n, c, cv):
    return get_norm_cv(m, n, c, cv, -1)


def get_pro_norm(m, n, c):
    return get_pro_norm_cv(m, n, c, pro_cv(m, n, c))


def get_obl_norm(m, n, c):
    return get_obl_norm_cv(m, n, c, obl_cv(m, n, c))
