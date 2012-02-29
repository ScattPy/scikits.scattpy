from numpy import *
from scikits.scattpy import *
P=ProlateSpheroid(ab=1.5,xv=1.0,m=1.33+0.07j)
LAB=Lab(P, pi/4)
RES=ebcm(LAB)

Csca_tm,Qsca_tm = LAB.get_Csca(RES.c_sca_tm)
print Csca_tm,Qsca_tm
Theta = linspace(0,pi,1000)
A = LAB.get_amplitude_matrix(RES.c_sca_tm,RES.c_sca_te,Theta,0)
S11g,S21_S11 = LAB.get_int_plr(A)
from matplotlib import pylab
pylab.semilogy(Theta*180/pi, S11g)
pylab.ylabel("S11/g")
pylab.xlabel("Theta")
pylab.title("Scattering field intencity")
pylab.show()
pylab.close()

pylab.plot(Theta*180/pi, S21_S11);
pylab.ylabel("S21/S11");
pylab.xlabel("Theta");
pylab.title("Scattering field degree of linear polarisation");
pylab.show()