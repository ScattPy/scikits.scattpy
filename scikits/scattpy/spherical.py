from numpy import *
from scipy.special.orthogonal import ps_roots
from scipy import special
from kmatrix import *

# (r, theta, phi)
class spherical_utilities(object):
	def __init__(self,ng,n,lab):
		self.n=n
		self.ng=ng
		self.set_ngauss()
		self.set_funcs_ang()
		self.set_all_layers(lab)

	def set_ngauss(self):
		k,w = ps_roots(self.ng)
		self.knots,self.weights = k*pi,w*pi
		self.thetas = kscalar(self.knots)
		self.sint=sin(self.thetas)
		self.cost=cos(self.thetas)
		self.tgt = tan(self.thetas)
		self.ctgt=1/self.tgt
	
	def set_funcs_ang(self):
		P = array([special.lpmn(self.n,self.n,ct) for ct in self.cost])
		self.data_Ang,self.data_Angd = P[:,0,:,:],P[:,1,:,:]
	
	def set_all_layers(self,lab):
		self.data_layers=[0]
		for bnd in lab.boundaries():
			R,Rd,Rdd = bnd.shape.R(self.knots)
			r=kscalar(R)
			rd=kscalar(Rd)
			rdd=kscalar(Rdd)
			rdr=rd/r
			r2rd2=r**2+rd**2
			Rad =[0,{},{}]
			Radd=[0,{},{}]
			kis=[0,bnd.k1,bnd.k2]
			for i in [1,2]:
				krs = kis[i]*r
				JY = array([special.sph_jnyn(self.n,kr) for kr in krs])[:,:,:]
				Rad[i] ={'j':JY[:,0,:], 'h':JY[:,0,:]+1j*JY[:,2,:]}
				Radd[i]={'j':JY[:,1,:], 'h':JY[:,1,:]+1j*JY[:,3,:]}

			self.data_layers.append({	'ki':kis,\
							'r':r,\
							'rd':rd,\
							'rdd':rdd,\
							'rdr':rdr,\
							'r2rd2':r2rd2,\
							'Rad':Rad,\
							'Radd':Radd })
	def set_layer_no(self,lay):
		self._lay=lay

	def _get_const(self,name):
		return self.data_layers[self._lay][name]

	def _get_r(self):return self._get_const("r")
	r   = property(fget=_get_r)

	def _get_rd(self):return self._get_const("rd")
	rd  = property(_get_rd)

	def _get_rdd(self):return self._get_const("rdd")
	rdd  = property(_get_rdd)

	def _get_rdr(self):return self._get_const("rdr")
	rdr  = property(_get_rdr)

	def _get_r2rd2(self):return self._get_const("r2rd2")
	r2rd2  = property(_get_r2rd2)

	def _get_ki(self):return self._get_const("ki")
	ki  = property(_get_ki)

	def Rad(self,m,ij,i):
		return kmatrix(self.data_layers[self._lay]['Rad'][i][ij][:,m:])

	def Radd(self,m,ij,i):
		return kmatrix(self.data_layers[self._lay]['Radd'][i][ij][:,m:])

	def Ang(self,m):
		return kmatrix(self.data_Ang[:,m,m:])

	def Angd(self,m):
		return kmatrix(self.data_Angd[:,m,m:])


def func_a(C,m,ij,i):
	return C.Rad(m,ij,i)^C.Ang(m)

def func_b(C,m,ij,i):
	ki=C.ki[i]
	return (ki*C.r*C.Radd(m,ij,i)^C.Ang(m)) \
		+(C.rdr*C.sint*C.Rad(m,ij,i)^C.Angd(m))

def func_d(C,m,ij,i):
	ki=C.ki[i]
	return (ki*C.r*C.cost*C.Radd(m,ij,i)^C.Ang(m)) \
		+(C.sint**2*C.Rad(m,ij,i)^C.Angd(m))

def func_e(C,m,ij,i):
	ki=C.ki[i]
	return ((ki*C.r*C.Radd(m,ij,i))+C.Rad(m,ij,i))^C.Ang(m)

def func_g(C,m,ij,i):
	ki=C.ki[i]
	return (ki*C.rd*C.Radd(m,ij,i)^C.Ang(m)) \
		-(C.sint*C.Rad(m,ij,i)^C.Angd(m))

def matA0(C,m,ij,i,coef):
	return ( (func_a(C,m,ij,i)).T * (coef*C.Ang(m)) ).integrate(C.weights)

def matA(C,m,ij,i):
	return matA0(C,m,ij,i,C.sint)

def matB(C,m,ij,i):
	return ( (func_b(C,m,ij,i)).T * (C.sint*C.Ang(m)) ).integrate(C.weights)

def matC(C,m,ij,i,e12, B=None):
	if B is None: B = matB(C,m,ij,i)
	ff = (1.-C.ctgt*C.rdr)*C.sint
	A = matA0(C,m,ij,i,coef=ff)
	return e12*B + (e12-1.)*A

def matD0(C,m,ij,i,coef):
	return ( (func_d(C,m,ij,i)).T * (coef*C.Ang(m)) ).integrate(C.weights)
def matE0(C,m,ij,i,coef):
	return ( (func_e(C,m,ij,i)).T * (coef*C.Ang(m)) ).integrate(C.weights)
def matG0(C,m,ij,i,coef):
	return ( (func_g(C,m,ij,i)).T * (coef*C.Ang(m)) ).integrate(C.weights)


def matD(C,m,ij,i,e12,B=None):
	if B is None: B = matB(C,m,ij,i)
	fd = C.rdr
	D0 = matD0(C,m,ij,i,coef=fd)
	return B + (e12-1.)*D0

def matE(C,m,ij,i,e12):
	fe = C.rd
	E0 = matE0(C,m,ij,i,coef=fe)
	return (e12-1.)*E0

def matF(C,m,ij,i,e12):
	fd = (C.rd*C.cost-C.r*C.sint)/C.r**2
	D0 = matD0(C,m,ij,i,coef=fd)
	return -(e12-1.)*D0

def matG(C,m,ij,i,e12,B=None):
	if B is None: B = matB(C,m,ij,i)
	fe = (C.rd*C.cost-C.r*C.sint)/C.r
	E0 = matE0(C,m,ij,i,coef=fe)
	return B - (e12-1.)*E0

def matA12(C,m,ij,i,e21,A=None):
	if A is None: A = matA(C,m,ij,i)
	fa = C.r*(C.rd*C.cost-C.r*C.sint)/C.r2rd2
	A0 = matA0(C,m,ij,i,coef=fa)
	return A - (e21-1)*A0

def matA14(C,m,ij,i,e21):
	fa = (C.r**2*C.rd)/C.r2rd2
	A0 = matA0(C,m,ij,i,coef=fa)
	return -(e21-1)*A0

def matA32(C,m,ij,i,e21):
	fa = (C.rd*C.sint+C.r*C.cost)*(C.rd*C.cost-C.r*C.sint)/(C.r*C.r2rd2)
	A0 = matA0(C,m,ij,i,coef=fa)
	return (e21-1)*A0

def matA34(C,m,ij,i,e21,A=None):
	if A is None: A = matA(C,m,ij,i)
	fa = C.rd*(C.rd*C.sint+C.r*C.cost)/C.r2rd2
	A0 = matA0(C,m,ij,i,coef=fa)
	return A + (e21-1)*A0

def matA22(C,m,ij,i,e21,B=None):
	if B is None: B=matB(C,m,ij,i)
	fg = C.rd*(C.rd*C.cost-C.r*C.sint)/C.r2rd2
	fa = (C.r**2 - C.r*C.rdd + 2*C.rd**2) \
	    *( C.r*(C.rd*C.cost-C.r*C.sint) + C.rd*(C.rd*C.sint+C.r*C.cost) ) \
	    /(C.r2rd2)**2
	A0 = matA0(C,m,ij,i,coef=fa)
	G0 = matG0(C,m,ij,i,coef=fg)
	return B - (e21-1)*(G0-A0)

def matA24(C,m,ij,i,e21):
	fg = C.r*C.rd**2/C.r2rd2
	fa = C.rd*(C.r**4 - 2*C.r**3*C.rdd + 2*C.r**2*C.rd**2 - C.rd**4) \
	    /(C.r2rd2)**2
	A0 = matA0(C,m,ij,i,coef=fa)
	G0 = matG0(C,m,ij,i,coef=fg)
	return -(e21-1)*(G0-A0)

def matA42(C,m,ij,i,e21):
	fg = (C.rd*C.cost-C.r*C.sint)**2/(C.r*C.r2rd2)
	fa = 2*(C.r**2 - C.r*C.rdd + 2*C.rd**2) \
	    *(C.rd*C.cost-C.r*C.sint)*(C.rd*C.sint+C.r*C.cost) \
	    /(C.r*C.r2rd2**2)
	A0 = matA0(C,m,ij,i,coef=fa)
	G0 = matG0(C,m,ij,i,coef=fg)
	return (e21-1)*(G0-A0)

def matA44(C,m,ij,i,e21,B=None):
	if B is None: B = matB(C,m,ij,i)
	fg = C.rd*(C.rd*C.cost-C.r*C.sint)/C.r2rd2
	fa =((C.r**3*C.rdd + C.r**2*C.rd**2 - C.r*C.rd**2*C.rdd +3*C.rd**4)*C.sint\
	    +(C.rdr)*C.cost*(C.r**4 - 2*C.r**3*C.rdd + 2*C.r**2*C.rd**2 - C.rd**4)) \
	    /(C.r2rd2)**2
	#fa = f_utils.f4(r,rd,rdd,sint,cost)
	A0 = matA0(C,m,ij,i,coef=fa)
	G0 = matG0(C,m,ij,i,coef=fg)
	return B + (e21-1)*(G0-A0)
