"""ScattPy is a SciKits package providing numerical methods of solving
light scattering by non-spherical particles problem."""
from methods import svm,ebcm,pmm
from laboratory import \
		Shape,\
		ShapeSphere,\
		ShapeSpheroid,\
		ShapeChebyshev,\
		Particle,\
		HomogeneousParticle,\
		Sphere,\
		ProlateSpheroid,\
		OblateSpheroid,\
		ChebParticle,\
		LayeredParticle,\
		Layered_EqShape_Particle,\
		Layered_EqShape_EqVolume_Particle,\
		EffMedium_Particle,\
		EffMedium_MaxwellGarnett_Particle,\
		EffMedium_InvMaxwellGarnett_Particle,\
		EffMedium_Bruggeman_Particle,\
		Lab
from version import __version__

