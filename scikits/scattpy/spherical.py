from numpy import *
from scipy.special.orthogonal import ps_roots
from scipy import special
import f_utils

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
		self.thetas = self.knots
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
			r,rd,rdd = bnd.shape.R(self.knots)
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
		return self.data_layers[self._lay]['Rad'][i][ij][:,m:]

	def Radd(self,m,ij,i):
		return self.data_layers[self._lay]['Radd'][i][ij][:,m:]

	def Ang(self,m):
		return self.data_Ang[:,m,m:]

	def Angd(self,m):
		return self.data_Angd[:,m,m:]


def matA0(C,m,jh,i,coef):
	Rad = C.Rad(m,jh,i)
	Angm = C.Ang(m)
	#func = lambda k:  outer( Rad[k]*Angm[k], coef[k]*Angm[k])
	#return mat_integrate(func)
	return f_utils.mat_a0(Rad,Angm,coef,C.weights)

def get_matX(func):
	def matX(C,m,jh,i,bnd):
	  Rad = C.Rad(m,jh,i)
	  Radd= C.Radd(m,jh,i)
	  Ang = C.Ang(m)
	  Angd= C.Angd(m)
	  return f_utils.mat_x(func,C.ki[i],bnd.e12,bnd.e21,\
                              Rad,Radd,Ang,Angd,C.r,C.rd,C.rdd,C.sint,C.cost,C.ctgt,C.weights)
	return matX

def get_matXX(func):
	def matX(C,m,jh,i,bnd):
	  Rad = C.Rad(m,jh,i)
	  Radd= C.Radd(m,jh,i)
	  Ang = C.Ang(m)
	  Angd= C.Angd(m)
	  return func(C.ki[i],bnd.e12,bnd.e21,\
                              Rad,Radd,Ang,Angd,C.r,C.rd,C.rdd,C.sint,C.cost,C.ctgt,C.weights)
	return matX

matA = get_matXX(f_utils.mat__a)
matB = get_matXX(f_utils.mat__b)
matC = get_matXX(f_utils.mat__c_tm)
matDtm = get_matXX(f_utils.mat__d_tm)
matEtm = get_matXX(f_utils.mat__e_tm)
matFtm = get_matXX(f_utils.mat__f_tm)
matGtm = get_matXX(f_utils.mat__g_tm)
matDte = get_matXX(f_utils.mat__d_te)
matEte = get_matXX(f_utils.mat__e_te)
matFte = get_matXX(f_utils.mat__f_te)
matGte = get_matXX(f_utils.mat__g_te)
matAa = get_matXX(f_utils.mat__aa)
matAb = get_matXX(f_utils.mat__ab)
matAc = get_matXX(f_utils.mat__ac)
matAd = get_matXX(f_utils.mat__ad)

def get_matXQ(func):
	def matXQ(C,m,jh1,jh2,bnd):
	  Rad1 = C.Rad(m,jh1,1)
	  Radd1= C.Radd(m,jh1,1)
	  Ang = C.Ang(m)
	  Angd= C.Angd(m)
	  Rad2 = C.Rad(m,jh2,2)
	  Radd2= C.Radd(m,jh2,2)
	  return func(bnd.k1,bnd.k2,bnd.e12,bnd.e21,\
                              Rad1,Radd1,Rad2,Radd2,Ang,Angd,\
			      C.r,C.rd,C.rdd,C.sint,C.cost,C.ctgt,C.weights)
	return matXQ

matQtm = get_matXQ(f_utils.mat_q_tm)
matQte = get_matXQ(f_utils.mat_q_te)
matQ11tm = get_matXQ(f_utils.mat_q11_tm)
matQ12tm = get_matXQ(f_utils.mat_q12_tm)
matQ21tm = get_matXQ(f_utils.mat_q21_tm)
matQ22tm = get_matXQ(f_utils.mat_q22_tm)
matQ11te = get_matXQ(f_utils.mat_q11_te)
matQ12te = get_matXQ(f_utils.mat_q12_te)
matQ21te = get_matXQ(f_utils.mat_q21_te)
matQ22te = get_matXQ(f_utils.mat_q22_te)
