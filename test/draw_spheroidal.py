import time

import matplotlib
from spheroidal_ebcm import SpheroidalEBCM
from spheroidal_pmm import SpheroidalPMM

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from pylab import load
from pylab import save

from spheroidal import *
from spheroidal_particles import *
from properties import *

import sys

m = 1.5
alpha = pi / 4

n_min = 6
n_max = 52

def generate_data(a_b,xl,type):
    y=[]
    execution_time=[]

    particle = Spheroid(m,a_b,type)
    particle.set_xl(xl)



    id = str(a_b) + " " + str(xl) + " " + str(type)

    for i in range(n_min,n_max,2):
        nmax = i
        print nmax
        start = time.time()
        svm = SpheroidalSVM(particle,nmax)
        b_sca = svm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, b_sca, nmax)[0]
        C_sca = getCsca(particle, b_sca, nmax)[0]
        delta=(C_ext-C_sca)/(C_ext+C_sca)
        execution_time.append(time.time() - start)
        y.append(delta)

    y= np.fabs(y)
    save("svm_delta"+id,y)
    save("svm_time"+id,execution_time)
    y=[]
    execution_time=[]

    for i in range(n_min,n_max,2):
        nmax = i
        print nmax
        start = time.time()
        ebcm = SpheroidalEBCM(particle,nmax)
        b_sca = ebcm.getMatrixSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, b_sca, nmax)[0]
        C_sca = getCsca(particle, b_sca, nmax)[0]
        delta=(C_ext-C_sca)/(C_ext+C_sca)
        execution_time.append(time.time() - start)
        y.append(delta)

    y= np.fabs(y)
    save("ebcm_m_delta"+id,y)
    save("ebcm_m_time"+id,execution_time)
    y=[]
    execution_time=[]

    for i in range(n_min,n_max,2):
        nmax = i
        print nmax
        start = time.time()
        ebcm = SpheroidalEBCM(particle,nmax)
        b_sca = ebcm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, b_sca, nmax)[0]
        C_sca = getCsca(particle, b_sca, nmax)[0]
        delta=(C_ext-C_sca)/(C_ext+C_sca)
        execution_time.append(time.time() - start)
        y.append(delta)

    y= np.fabs(y)
    save("ebcm_delta"+id,y)
    save("ebcm_time"+id,execution_time)
    y=[]
    execution_time=[]

#    for i in range(n_min,n_max,2):
#        nmax = i
#        print nmax
#        start = time.time()
#        pmm = SpheroidalPMM(particle,c2,c1,nmax)
#        b_sca = pmm.getMatrixSolution(TMInputWave(alpha))[0]
#        C_ext = getCext(particle, alpha, k, b_sca, nmax)[0]
#        C_sca = getCsca(k, b_sca, nmax)[0]
#        delta=(C_ext-C_sca)/(C_ext+C_sca)
#        execution_time.append(time.time() - start)
#        y.append(delta)
#
#    y= np.fabs(y)
#    save("pmm_m_delta",y)
#    save("pmm_m_time",execution_time)
#    y=[]
#    execution_time=[]
#
#    for i in range(n_min,n_max,2):
#        nmax = i
#        print nmax
#        start = time.time()
#        pmm = SpheroidalPMM(particle,c2,c1,nmax)
#        b_sca = pmm.getSolution(TMInputWave(alpha))[0]
#        C_ext = getCext(particle, alpha, k, b_sca, nmax)[0]
#        C_sca = getCsca(k, b_sca, nmax)[0]
#        delta=(C_ext-C_sca)/(C_ext+C_sca)
#        execution_time.append(time.time() - start)
#        y.append(delta)
#
#    y= np.fabs(y)
#    save("pmm_delta",y)
#    save("pmm_time",execution_time)

def plot_graphics(a_b,xl,type):
    particle = Spheroid(m,a_b,type)
    particle.set_xl(xl)
    id = str(a_b) + " " + str(xl) + " " + str(type)

    svm_delta=load('svm_delta'+id+'.npy')
    svm_time=load('svm_time'+id+'.npy')
    ebcm_m_delta=load('ebcm_m_delta'+id+'.npy')
    ebcm_m_time=load('ebcm_m_time'+id+'.npy')
    ebcm_delta=load('ebcm_delta'+id+'.npy')
    ebcm_time=load('ebcm_time'+id+'.npy')
    x=np.arange(n_min,n_max,2)

    title = 'm = '+str(particle.m)+' a/b = ' + str(a_b) + ' type = '+str(type)+ ' xl = ' + str(xl) + ' xv = '+str(particle.xv) + ' alpha = ' + str(alpha)

    fig = plt.figure(figsize=(10,10))
    plt.title(title)
    plt.grid(True)
    plt.plot(x,svm_delta,'k-',label='SVM')
    plt.plot(x,ebcm_m_delta,'b--',label='EBCMm')
    plt.plot(x,ebcm_delta,'r-.',label='EBCMi')
    plt.yscale('log')
    plt.xlabel('N')
    plt.ylabel('Error')
    plt.legend(loc=4)
    plt.savefig('comp_delta'+id+'.png')
    plt.clf()

    fig = plt.figure(figsize=(10,10))
    plt.title(title)
    plt.grid(True)
    plt.plot(x,svm_time,'k-',label='SVM')
    plt.plot(x,ebcm_m_time,'b--',label='EBCMm')
    plt.plot(x,ebcm_time,'r-.',label='EBCMi')
    plt.yscale('log')
    plt.xlabel('N')
    plt.ylabel('Time')
    plt.legend(loc=4)
    plt.savefig('comp_time'+id+'.png')


if __name__ == '__main__':
    a_b = float(sys.argv[1])
    xl = float(sys.argv[2])
    type = int(sys.argv[3])
    generate_data(a_b,xl,type)
    plot_graphics(a_b,xl,type)


