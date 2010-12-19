# Core functions
from numpy import *
from scipy.special import sph_jn,sph_jnyn,lpmn
from scipy.special.orthogonal import ps_roots
from scipy.misc.common import factorial
from math import atan,acos,asin
#from IPython.Debugger import Tracer; debug_here = Tracer()
import f_utils

knots   = 0
weights = 0
ngauss  = 0
sint = 0
cost = 0
tgt  = 0
ctgt = 0
g = 0
r = 0
rd= 0
rdd=0
rdr=0
Pna=0

def set_consts(Rs):
	global r,rd,rdd,rdr
	r,rd,rdd = Rs
	rdr = rd/r

def get_Jn(n,x):
	return array([sph_jn(n,xl) for xl in x])

def get_JnHn(n,x):
	JnYn = array([sph_jnyn(n,xl) for xl in x])
	return JnYn[:,:2,:],JnYn[:,:2,:]+1j*JnYn[:,2:,:]

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

def set_ngauss(ng):
	global ngauss,knots,weights,cost,sint,tgt,ctgt
	ngauss = ng
	x,w = ps_roots(ng)
	knots = pi*x
	weights = w*pi
	sint = sin(knots)
	cost = cos(knots)
	tgt  = tan(knots)
	ctgt = 1/tgt
	return knots

def mat_integrate(func):
	kr = xrange(len(weights))
	sum=0
	for k in kr:
		sum += func(k)*weights[k]
	return mat(sum)

def matA0(n,m,fRad,fAng,coef):
	Rad = fRad[:,0,m:]
	Angm = fAng[:,0,m,m:]
	#func = lambda k:  outer( Rad[k]*Angm[k], coef[k]*Angm[k])
	#return mat_integrate(func)
	return f_utils.mat_a0(Rad,Angm,coef,weights)

def matA(n,m,fRad,fAng):
	return matA0(n,m,fRad,fAng,sint)

def matB(n,m,fRad,fAng,ki):
	Rad = fRad[:,0,m:]
	Radd= fRad[:,1,m:]
	Angm = fAng[:,0,m,m:]
	Angmd= fAng[:,1,m,m:]
	#func = lambda k: \
	#   outer(ki*r[k]*Radd[k]*Angm[k]+rdr[k]*sint[k]*Rad[k]*Angmd[k],\
	#         sint[k]*Angm[k] )
	#return mat_integrate(func)
	return f_utils.mat_b(ki,Rad,Radd,Angm,Angmd,r,rdr,sint,weights)

def matC(n,m,fRad,fAng,ki,e12, B=None):
	if B is None: B = matB(n,m,fRad,fAng,ki)
	ff = (1.-ctgt*rdr)*sint
	A = matA0(n,m,fRad,fAng,coef=ff)
	return e12*B + (e12-1.)*A

def matD0(n,m,fRad,fAng,ki,coef):
	Rad = fRad[:,0,m:]
	Radd= fRad[:,1,m:]
	Angm = fAng[:,0,m,m:]
	Angmd= fAng[:,1,m,m:]
	#func = lambda k: \
	#   outer(ki*r[k]*cost[k]*Radd[k]*Angm[k]+sint[k]**2*Rad[k]*Angmd[k],\
	#         Angm[k]*coef[k] )
	#return mat_integrate(func)
	return f_utils.mat_d0\
			(ki,Rad,Radd,Angm,Angmd,r,sint,cost,coef,weights)

def matE0(n,m,fRad,fAng,ki,coef):
	Rad = fRad[:,0,m:]
	Radd= fRad[:,1,m:]
	Angm = fAng[:,0,m,m:]
	Angmd= fAng[:,1,m,m:]
	#func = lambda k: \
	#   outer((ki*r[k]*Radd[k]+Rad[k])*Angm[k], Angm[k]*coef[k] )
	#return mat_integrate(func)
	return f_utils.mat_e0(ki,Rad,Radd,Angm,r,coef,weights)

def matG0(n,m,fRad,fAng,ki,coef):
	Rad = fRad[:,0,m:]
	Radd= fRad[:,1,m:]
	Angm = fAng[:,0,m,m:]
	Angmd= fAng[:,1,m,m:]
	#func = lambda k: \
	#  outer(ki*rd[k]*Radd[k]*Angm[k] - sint[k]*Rad[k]*Angmd[k],\
	#         Angm[k]*coef[k] )
	#return mat_integrate(func)
	return f_utils.mat_g0(ki,Rad,Radd,Angm,Angmd,r,rd,sint,coef,weights)

def matD(n,m,fRad,fAng,ki,e12,B=None):
	if B is None: B = matB(n,m,fRad,fAng,ki)
	fd = rdr
	D0 = matD0(n,m,fRad,fAng,ki,coef=fd)
	return B + (e12-1.)*D0

def matE(n,m,fRad,fAng,ki,e12):
	fe = rd
	E0 = matE0(n,m,fRad,fAng,ki,coef=fe)
	return (e12-1.)*E0

def matF(n,m,fRad,fAng,ki,e12):
	fd = (rd*cost-r*sint)/r**2
	D0 = matD0(n,m,fRad,fAng,ki,coef=fd)
	return -(e12-1.)*D0

def matG(n,m,fRad,fAng,ki,e12,B=None):
	if B is None: B = matB(n,m,fRad,fAng,ki)
	fe = (rd*cost-r*sint)/r
	E0 = matE0(n,m,fRad,fAng,ki,coef=fe)
	return B - (e12-1.)*E0

def matA12(n,m,fRad,fAng,ki,e21,A=None):
	if A is None: A = matA(n,m,fRad,fAng)
	fa = r*(rd*cost-r*sint)/(r**2+rd**2)
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	return A - (e21-1)*A0

def matA14(n,m,fRad,fAng,ki,e21):
	fa = (r**2*rd)/(r**2+rd**2)
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	return -(e21-1)*A0

def matA32(n,m,fRad,fAng,ki,e21):
	fa = (rd*sint+r*cost)*(rd*cost-r*sint)/(r*(r**2+rd**2))
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	return (e21-1)*A0

def matA34(n,m,fRad,fAng,ki,e21,A=None):
	if A is None: A = matA(n,m,fRad,fAng)
	fa = rd*(rd*sint+r*cost)/(r**2+rd**2)
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	return A + (e21-1)*A0

def matA22(n,m,fRad,fAng,ki,e21,B=None):
	if B is None: B=matB(n,m,fRad,fAng,ki)
	fg = rd*(rd*cost-r*sint)/(r**2+rd**2)
	fa = (r**2 - r*rdd + 2*rd**2) \
	    *( r*(rd*cost-r*sint) + rd*(rd*sint+r*cost) ) \
	    /(r**2+rd**2)**2
	#fa = f_utils.f1(r,rd,rdd,sint,cost)
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	G0 = matG0(n,m,fRad,fAng,ki,coef=fg)
	return B - (e21-1)*(G0-A0)

def matA24(n,m,fRad,fAng,ki,e21):
	fg = r*rd**2/(r**2+rd**2)
	fa = rd*(r**4 - 2*r**3*rdd + 2*r**2*rd**2 - rd**4) \
	    /(r**2+rd**2)**2
	#fa = f_utils.f2(r,rd,rdd) * rd
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	G0 = matG0(n,m,fRad,fAng,ki,coef=fg)
	return -(e21-1)*(G0-A0)

def matA42(n,m,fRad,fAng,ki,e21):
	fg = (rd*cost-r*sint)**2/(r*(r**2+rd**2))
	fa = 2*(r**2 - r*rdd + 2*rd**2) \
	    *(rd*cost-r*sint)*(rd*sint+r*cost) \
	    /(r*(r**2+rd**2)**2)
	#fa = f_utils.f3(r,rd,rdd,sint,cost) / r
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	G0 = matG0(n,m,fRad,fAng,ki,coef=fg)
	return (e21-1)*(G0-A0)

def matA44(n,m,fRad,fAng,ki,e21,B=None):
	if B is None: B = matB(n,m,fRad,fAng,ki)
	fg = rd*(rd*cost-r*sint)/(r**2+rd**2)
	fa =((r**3*rdd + r**2*rd**2 - r*rd**2*rdd +3*rd**4)*sint\
	    +(rd/r)*cost*(r**4 - 2*r**3*rdd + 2*r**2*rd**2 - rd**4)) \
	    /(r**2+rd**2)**2
	#fa = f_utils.f4(r,rd,rdd,sint,cost)
	A0 = matA0(n,m,fRad,fAng,coef=fa)
	G0 = matG0(n,m,fRad,fAng,ki,coef=fg)
	return B + (e21-1)*(G0-A0)

def get_delta(Q_ext,Q_sca):
	delta = abs(Q_ext-Q_sca)/abs(Q_ext+Q_sca)
	return delta


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
