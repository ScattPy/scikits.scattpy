import spheroidal

from numpy import pi, real, power, cos, absolute

#according to (67)
def getCext(particle, alpha, k, b, nmax):
    sum = 0
    for i in range(nmax):
        l = i + 1
        c = particle.c
        cv = spheroidal.get_cv(1, l, c, particle.type)
        sum += power(1j, -l) * b[i] *\
               spheroidal.ang1_cv(1, l, c, cv, particle.type, cos([alpha]))[0]
    return 4 * pi / (k * k) * real(sum)

#according to (68)
def getCsca(k, b, nmax):
    sum = 0
    for i in range(nmax):
        sum += 2 * absolute(b[i] * b[i])
    return pi / (k * k) * sum
  