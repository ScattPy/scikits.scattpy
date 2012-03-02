import time

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from spheroidal import *
from spheroidal_particles import *
from properties import *

def plot_graphics():
    y=[]
    execution_time=[]
    c1 = 1.3
    c2 = 2
    alpha = pi / 4
    k = 1
    particle = ProlateSpheroid(psi=2,c=c1,derivative=0,eps=1)
    for i in range(6,52,2):
        nmax = i
        print nmax
        start = time.time()
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        b_sca = getSolution(svm,particle, TMInputWave(alpha), c1, c2, nmax)[0]
        C_ext = -getCext(particle, alpha, k, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
        delta=(C_ext-C_sca)/(C_ext+C_sca)
        execution_time.append(time.time() - start)
        y.append(delta)
    x = np.arange(6,52,2)
    y= np.fabs(y)
    plt.plot(x,np.log10(y))
    plt.savefig('error')
    plt.clf()
    plt.plot(x,np.log10(execution_time))
    plt.savefig('time')

if __name__ == '__main__':
    plot_graphics()


