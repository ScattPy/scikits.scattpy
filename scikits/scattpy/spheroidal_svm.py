from numpy import bmat
import spheroidal

#according to (85)
def get_fullB(particle, c, nmax):
    return bmat([[spheroidal.get_A10(particle, c, nmax)], [spheroidal.get_A20(particle, c, nmax)]])


def get_fullA(particle, c1, c2, nmax):
    A11 = spheroidal.get_A11(particle, c1, nmax)
    A12 = spheroidal.get_A12(particle, c1, c2, nmax)
    A21 = spheroidal.get_A21(particle, c1, nmax)
    A22 = spheroidal.get_A22(particle, c1, c2, nmax)

    return bmat([[A11, A12], [A21, A22]])