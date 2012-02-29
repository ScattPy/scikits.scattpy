"""ScattPy is a SciKits package providing numerical methods of solving
light scattering by non-spherical particles problem."""
from methods import svm, ebcm, pmm
from laboratory import\
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
Layered_EqShapeEqVol_Particle,\
LayeredConfocalSpheroid,\
EMT_Particle,\
EMT_MGarn_Particle,\
EMT_IMGarn_Particle,\
EMT_Brugg_Particle,\
Lab
from version import __version__

