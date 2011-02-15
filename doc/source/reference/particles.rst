.. _particles:

######
Particles
######

.. currentmodule:: scikits.scattpy

A particle is defined by its surface shape and the complex refractive index of its medium.
In the case of a multilayered particle, each layer's refractive index and surface shape
should be defined.

Homogeneous particles
=====================

A homogeneous particle can be defined as an instance of the :class:`HomogeneousParticle`.
There are also several predefined classes for the most often used partiles.

.. autosummary::
   :toctree: generated/

   HomogeneousParticle
   Sphere
   ProlateSpheroid
   OblateSpheroid
   ChebParticle

Layered particles
=====================

A layered particle can be defined as an instance of the :class:`LayeredParticle`.
There are also several predefined classes for the most often used multilayered partiles.

.. autosummary::
   :toctree: generated/

   LayeredParticle
   Layered_EqShape_Particle
   Layered_EqShapeEqVol_Particle
   LayeredConfocalSpheroid

Effective Medium Theory	
=====================

Inhomogeneous scatterers can also be treated with the effective medium theory (EMT).
An EMT particle is an instance of a class derrivered from the :class:`EMT_Particle`.
There are several classes implementing EMT with the most widely used EMT rules.

.. autosummary::
   :toctree: generated/

   EMT_Particle
   EMT_MGarn_Particle
   EMT_IMGarn_Particle
   EMT_Brugg_Particle
