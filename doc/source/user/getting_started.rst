.. _getting_started:

***************
Getting started
***************

.. _getting_started_installation:

Installation
============

To install ScattPy one needs:

* Python 2.5 - 2.7
* F2PY with a configured FORTRAN77 compiler
* Python packages:

  - NumPy
  - SciPy
  - NumdiffTools

Ubuntu
------

.. code-block:: bash

   $ sudo apt-get install python-numpy python-scipy python-dev python-setuptools ipython
   $ sudo easy_install numdifftools
   $ sudo easy_install scikits.scattpy

Windows
-------

Instructions for installing ScattPy under Windows operating system will be added later

.. _getting_started_sample_session:

Sample Session
==============

The most simple way of trying ScattPy is to run an interactive python session with IPython.

.. code-block:: bash

   $ ipython

Import required packages.

.. ipython::

   In [1]: from numpy import * ;

   In [2]: from scikits.scattpy import * ;

Define a prolate spheroidal particle with the size parameter :math:`x_{\rm V}=2`, aspect ratio :math:`a/b=4` and complex refractive index :math:`m`.

.. ipython::

   In [3]: P = ProlateSpheroid(ab=4., xv=2., m=1.33+0.2j)

Define a laboratory object that desribes the incident wave and the scatterer. We will consider a plane wave propagating at the angle :math:`\alpha=\pi/4` to the particle symmetry axis.

.. ipython::

   In [1]: LAB = Lab(P, alpha=pi/4)

Now we can use any of the ScattPy's numerical methods to obtain scattering properties of the defined model. In this example we apply the extended boundary conditions method (EBCM) that is best suited for homogeneous spheroids.

.. ipython::

   In [1]: RES = ebcm(LAB)

The returned `RES` object contains the expansion coefficients of the scattered field for the TM and TE modes. These coefficients can be used to obtain optical characteristics such as scattering cross-sections and efficiency factors

.. ipython::

   In [1]: Csca_tm,Qsca_tm = LAB.get_Csca(RES.c_sca_tm) ;

   In [2]: print Csca_tm,Qsca_tm

One can also calculate amplitude ans scattering matrix elements for defined ranges of scattering angles, e.g. :math:`\Theta\in[0,\pi],\;\varphi=0`:

.. ipython::

   In [2]: Theta = linspace(0,pi,1000) ;

   In [3]: A = LAB.get_amplitude_matrix(RES.c_sca_tm,RES.c_sca_te,Theta,0) ;
 
   In [4]: S11g,S21_S11 = LAB.get_int_plr(A) ;

Using ScattPy together with Python data visualisation packages one can obtain plots of the scattering matrix elements.

.. ipython::

   In [1]: from matplotlib import pylab

   In [5]: pylab.semilogy(Theta*180/pi, S11g);

   In [5]: pylab.ylabel("S11/g");

   In [5]: pylab.xlabel("Theta");

   In [5]: pylab.title("Scattering field intencity");

   @savefig pylab/getting_started_S11.png
   In [5]: pylab.show()

.. ipython::

   In [5]: pylab.close()

   In [5]: pylab.plot(Theta*180/pi, S21_S11);

   In [5]: pylab.ylabel("S21/S11");

   In [5]: pylab.xlabel("Theta");

   In [5]: pylab.title("Scattering field degree of linear polarisation");

   @savefig pylab/getting_started_S21.png
   In [5]: pylab.show()
