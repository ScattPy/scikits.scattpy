import time

import matplotlib
from spheroidal_ebcm import SpheroidalEBCM

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from pylab import load
from pylab import save

from spheroidal import *
from spheroidal_particles import *
from properties import *

def generate_data():
    y=[]
    execution_time=[]
    c1 = 1
    c2 = 6
    alpha = pi / 4
    k = 1
    particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1.5)
    n_min = 6
    n_max = 52

    for i in range(n_min,n_max,2):
        nmax = i
        print nmax
        start = time.time()
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        b_sca = svm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, k, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
        delta=(C_ext-C_sca)/(C_ext+C_sca)
        execution_time.append(time.time() - start)
        y.append(delta)

    y= np.fabs(y)
    save("svm_delta",y)
    save("svm_time",execution_time)
    y=[]
    execution_time=[]

    for i in range(n_min,n_max,2):
        nmax = i
        print nmax
        start = time.time()
        ebcm = SpheroidalEBCM(particle,c2,c1,nmax)
        b_sca = ebcm.getMatrixSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, k, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
        delta=(C_ext-C_sca)/(C_ext+C_sca)
        execution_time.append(time.time() - start)
        y.append(delta)

    y= np.fabs(y)
    save("ebcm_m_delta",y)
    save("ebcm_m_time",execution_time)
    y=[]
    execution_time=[]

    for i in range(n_min,n_max,2):
        nmax = i
        print nmax
        start = time.time()
        ebcm = SpheroidalEBCM(particle,c2,c1,nmax)
        b_sca = ebcm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, k, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
        delta=(C_ext-C_sca)/(C_ext+C_sca)
        execution_time.append(time.time() - start)
        y.append(delta)

    y= np.fabs(y)
    save("ebcm_delta",y)
    save("ebcm_time",execution_time)
    #x = np.arange(n_min,n_max,2)

def plot_graphics():
    svm_delta=load('svm_delta.npy')
    svm_time=load('svm_time.npy')
    ebcm_m_delta=load('ebcm_m_delta.npy')
    ebcm_m_time=load('ebcm_m_time.npy')
    ebcm_delta=load('ebcm_delta.npy')
    ebcm_time=load('ebcm_time.npy')
    x=np.arange(6,52,2)

    fig = plt.figure(figsize=(10,10))
    plt.title('c1=1,c2=6,psi=2,eps=1.5,alpha=pi/4')
    plt.grid(True)
    plt.plot(x,svm_delta,'k-',label='SVM')
    plt.plot(x,ebcm_m_delta,'k--',label='EBCMm')
    plt.plot(x,ebcm_delta,'k-.',label='EBCMi')
    plt.yscale('log')
    plt.xlabel('N')
    plt.ylabel('Error')
    plt.legend(loc=4)
    plt.savefig('comp_delta.png')
    plt.clf()

    fig = plt.figure(figsize=(10,10))
    plt.title('c1=1,c2=6,psi=2,eps=1.5,alpha=pi/4')
    plt.grid(True)
    plt.plot(x,svm_time,'k-',label='SVM')
    plt.plot(x,ebcm_m_time,'k--',label='EBCMm')
    plt.plot(x,ebcm_time,'k-.',label='EBCMi')
    plt.yscale('log')
    plt.xlabel('N')
    plt.ylabel('Time')
    plt.legend(loc=4)
    plt.savefig('comp_time.png')


if __name__ == '__main__':
    plot_graphics()


