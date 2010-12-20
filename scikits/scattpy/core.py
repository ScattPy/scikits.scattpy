# Core functions
from numpy import *
from scipy.special import sph_jn,sph_jnyn,lpmn
from scipy.misc.common import factorial
from math import atan,acos,asin
#from IPython.Debugger import Tracer; debug_here = Tracer()

Pna=0


def get_Pmn(m,n,x):
	return array([lpmn(m,n,xl) for xl in x])

def norm(m,l):
	return sqrt((2*l+1)/2. *factorial(l-m)/factorial(l+m))

def get_Pmn_normed(m,n,x):
	P = get_Pmn(m,n,x)
	for M in xrange(m+1):
		for L in xrange(n+1):
			P[:,:,M,L] = P[:,:,M,L]*norm(M,L)
	return P

#------------- Returning vector fields -------------------------------------

def coord_cartesian2spherical(cPoint):
	x,y,z = list(cPoint)
	r=sqrt(x**2+y**2+z**2)
	if r==0:
	   t=p=0.0
	else:
	   t=acos(z/r)
	   if x==0:
	      p=asin(sign(y))
	   else:
	      p=atan(y/abs(x))
	      if x<0:
	         p = pi-p
	return array([r,t,p])

def coord_spherical2cartesian(sPoint):
	r,t,p = list(sPoint)
	x=r*sin(t)*cos(p)
	y=r*sin(t)*sin(p)
	z=r*cos(t)
	return array([x,y,z])

def vector_cartesian2spherical(cVec,sPoint):
	r,t,p = list(sPoint)
	matr = matrix(\
	  [[sin(t)*cos(p), sin(t)*sin(p),  cos(t)],\
	   [cos(t)*cos(p), cos(t)*sin(p), -sin(t)],\
	   [      -sin(p),        cos(p),       0]])
	sVec = array((matr*matrix(cVec).T).T)[0]
	return sVec

def vector_spherical2cartesian(sVec,sPoint):
	r,t,p = list(sPoint)
	matr = matrix(\
	  [[sin(t)*cos(p), cos(t)*cos(p), -sin(p)],\
	   [sin(t)*sin(p), cos(t)*sin(p),  cos(p)],\
	   [       cos(t),        sin(t),       0]])
	cVec = array((matr*matrix(sVec).T).T)[0]
	return cVec

def get_vector(cPoint,c_sca,k):
	sPoint = coord_cartesian2spherical(cPoint)
	r,t,p = list(sPoint)
	sVec = array([0.,0.,0.])

	m=1

	n=len(c_sca)/2+m-1
	P,Pd = get_Pmn(m,n,[cos(t)])[0]
	Pml  = P[m,m:n+1]
	Pdml = Pd[m,m:n+1]
	Bess,Hank = get_JnHn(n,[k*r])
	Hankl = Hank[0,0,m:]
	Hankdl= Hank[0,1,m:]
	a = c_sca[:n-m+1]
	b = c_sca[n-m+1:]
	l = arange(m,n+1)

	if t==0 or t==pi:
	   vec_r=vec_p=vec_t=0
	else:
	   vec_r = (1./r * m*sin(m*p) * sum( a*Hankl*Pml ) ).real
	   vec_t =-(1./(r*sin(t)) * m*sin(m*p) \
	            * sum( (cos(t)*a+b)*Hankl*Pml ) ).real
	   vec_p = (-1./r *sin(t)*cos(m*p) \
	            * sum( a*Pml*(Hankl+k*r*Hankdl) ) \
	            + sin(t)*cos(m*p)* sum( a*Hankl*Pml )\
	            - cos(m*p)* sum( (cos(t)*a+b)*Hankl*Pdml ) ).real

	sVec += array([vec_r,vec_t,vec_p])
	cVec = vector_spherical2cartesian(sVec,sPoint)
	return cVec
