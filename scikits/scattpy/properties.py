import spheroidal

from numpy import pi, real, power, cos, absolute

#according to (67)
def getCext(particle, alpha, b, nmax):
    sum = 0
    k = particle.k
    for l in range(1, nmax+1):
        c = particle.c1
        cv = spheroidal.get_cv(1, l, c, particle.type)
        sum += power(1j, -l) * b[l-1] *\
               spheroidal.ang1_cv(1, l, c, cv, particle.type, cos([alpha]))[0]
    return 4 * pi / (k * k) * real(sum)

#according to (68)
def getCsca(particle, b, nmax):
    sum = 0
    k = particle.k
    for i in range(nmax):
        sum += 2 * absolute(b[i] * b[i])
    return pi / (k * k) * sum
  