"""A package for handling matrices and scalars that have different values for
different values of a parameter.
The package is usefull for numerical integration."""
from numpy import array,matrix,isscalar,shape
import numpy
import kmatrix_f as F
class kmatrix(object):
	def __init__(self,arr):
		arr=array(arr)
		if len(arr.shape)==2:
			self._a = F.vec2kmatrix(arr)
		elif len(arr.shape)==3:
			self._a = array(arr)
		elif issubclass(arr.__class__,kmatrix):
			self=arr
		else:
			raise TypeError,"Unexpected dimension of array: "+\
					arr.shape

	def __mul__(self,other):
		if isscalar(other):
			return kmatrix(self._a*other)
		elif issubclass(other.__class__,kmatrix):
			return kmatrix(F.mat_mat(self._a,other._a))
		elif issubclass(other.__class__,numpy.ndarray) and len(shape(other))<3:
			other = matrix(other)
			return kmatrix(F.mat_fixmat(self._a,other))
		else:
			return NotImplemented
			#raise "Unexpected type: "+str(other.__class__)

	def __rmul__(self,other):
		if isscalar(other):
			return kmatrix(self._a*other)
		#elif len(other.shape)==1 and other.shape[0]==self._a.shape[0]:
		#	return kmatrix(F.coef_mat(other,self._a))
		elif issubclass(other.__class__,numpy.ndarray) and len(shape(other))<3:
			other = matrix(other)
			return kmatrix(F.fixmat_mat(other,self._a))
		else:
			return NotImplemented
	
	def __add__(self,other):
		if issubclass(other.__class__,kmatrix):
			return kmatrix(self._a+other._a)
		elif issubclass(other.__class__,numpy.ndarray) and len(shape(other))<3:
			other = matrix(other)
			return kmatrix(F.sum_fixmat_mat(other,self._a))
		else:
			return NotImplemented

	def __sub__(self,other):
		if issubclass(other.__class__,kmatrix):
			return kmatrix(self._a-other._a)
		elif issubclass(other.__class__,numpy.ndarray) and len(shape(other))<3:
			other = matrix(other)
			return kmatrix(F.sum_fixmat_mat(-other,self._a))
		else:
			return NotImplemented

	def __radd__(self,other):
		if issubclass(other.__class__,kmatrix):
			return kmatrix(self._a+other._a)
		elif issubclass(other.__class__,numpy.ndarray) and len(shape(other))<3:
			other = matrix(other)
			return kmatrix(F.sum_fixmat_mat(other,self._a))
		else:
			return NotImplemented

	def __rsub__(self,other):
		if issubclass(other.__class__,kmatrix):
			return kmatrix(other._a-self._a)
		elif issubclass(other.__class__,numpy.ndarray) and len(shape(other))<3:
			other = matrix(other)
			return kmatrix(F.sum_fixmat_mat(other,-self._a))
		else:
			return NotImplemented

	def __xor__(self,other):
		if not issubclass(other.__class__,kmatrix):
			raise TypeError,"Unexpected type: "+str(other.__class__)
		elif self._a.shape != other._a.shape:
			raise TypeError,"Shapes don't match"
		else:
			return kmatrix(self._a*other._a)

	def _transpose(self):
		return kmatrix(F.mat_transpose(self._a))

	T = property(_transpose)

	def integrate(self,weights):
		return matrix(F.mat_integrate(self._a,weights))


class kscalar(object):
	def __init__(self,arr):
		arr=array(arr)
		if len(arr.shape)==0:
			self = arr
		elif len(arr.shape)==1:
			self._a = arr
		elif issubclass(arr.__class__,karray):
			self=arr
		else:
			raise TypeError,"Unexpected dimension of array"

	def __mul__(self,other):
		if isscalar(other):
			return kscalar(self._a*other)
		elif issubclass(other.__class__,kscalar):
			return kscalar(self._a*other._a)
		elif issubclass(other.__class__,kmatrix):
			return kmatrix(F.coef_mat(self._a,other._a))
		elif issubclass(other.__class__,matrix):
			return 
		else:
			return NotImplemented

	def __div__(self,other):
		if isscalar(other):
			return kscalar(self._a/other)
		elif issubclass(other.__class__,kscalar):
			return kscalar(self._a/other._a)
		else:
			return NotImplemented

	def __rmul__(self,other):
		if isscalar(other):
			return kscalar(self._a*other)
		elif issubclass(other.__class__,kmatrix):
			return kmatrix(F.coef_mat(self._a,other._a))
		else:
			return NotImplemented
	
	def __rdiv__(self,other):
		if isscalar(other):
			return kscalar(other/self._a)
		elif issubclass(other.__class__,kmatrix):
			return kmatrix(F.coef_mat(1./self._a,other._a))
		else:
			return NotImplemented
	
	def __add__(self,other):
		if issubclass(other.__class__,kscalar):
			return kscalar(self._a+other._a)
		elif isscalar(other):
			return kscalar(self._a+other)
		else:
			return NotImplemented

	def __sub__(self,other):
		if issubclass(other.__class__,kscalar):
			return kscalar(self._a-other._a)
		elif isscalar(other):
			return kscalar(self._a-other)
		else:
			return NotImplemented

	def __radd__(self,other):
		if isscalar(other):
			return kscalar(self._a+other)
		else:
			return NotImplemented

	def __rsub__(self,other):
		if isscalar(other):
			return kscalar(other-self._a)
		else:
			return NotImplemented

	def __pow__(self,power):
		return kscalar(self._a**power)

	def __iter__(self):
		return iter(self._a)

	def integrate(self,weights):
		return sum(self._a*weights)

def kfunc(numpyfunc):
	def func(x):
		if issubclass(x.__class__,kscalar):
			return kscalar(numpyfunc(x._a))
		elif issubclass(x.__class__,kmatrix):
			return kmatrix(numpyfunc(x._a))
		else:
			return numpyfunc(x)
	return func

sin = kfunc(numpy.sin)
cos = kfunc(numpy.cos)
tan = kfunc(numpy.tan)
sqrt= kfunc(numpy.sqrt)


		


