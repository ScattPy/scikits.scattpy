import core
from numpy import *
from scipy.special import sph_jn,sph_jnyn,lpmn
from scipy.misc.common import factorial
#from IPython.Debugger import Tracer; debug_here = Tracer()

class ShapeSpheroid:
	def __init__(self,ab,xv):
		self.ab = ab
		self.xv = xv
	def __str__(self):
		return "spheroid with ab="+str(self.ab)+", xv="+str(self.xv)
	def copy(self):
		return self.__class__(self.ab,self.xv)
	def R(self,theta):
		sint = sin(theta)
		cost = cos(theta)
		xb=self.xv/self.ab**(1./3.)
		ba=1./self.ab
		r=xb/sqrt(1.+(ba**2-1.)*cost**2)
		rd=xb*(ba**2-1.)*sint*cost/(sqrt(1.+(ba**2-1.)*cost**2))**3

		f = xb * (ba**2 -1) * sint * cost
		g = (sqrt(1. + (ba**2 - 1) * cost**2))**3
		f1 = xb * (ba**2 -1) * (cost**2 - sint**2)
		g1 = - 3./2. * sqrt(1. + (ba**2 - 1) * cost**2) \
		             * (ba**2 - 1.) * 2 * cost * sint
		rdd = f1/g - f*g1/g**2
		return r,rd,rdd
	def _get_a(self):
		return self.ab*self.b
	def _get_b(self):
		return self.xv/self.ab**(1./3.)
	a = property(fget=_get_a)
	b = property(fget=_get_b)

class ShapeChebParticle:
	def __init__(self,N,eps,xv):
		self.N = N
		self.eps = eps
		self.xv = xv
	def __str__(self):
		return "Cheb"+str(self.N)+" with eps="+str(self.eps)+", xv="+str(self.xv)
	def copy(self):
		return self.__class__(self.N,self.eps,self.xv)
	def R(self,theta):
		xv = self.xv
		eps = self.eps
		N = self.N

		r = xv*(1+eps*cos(N*theta))
		rd = -xv*eps*N*sin(N*theta)
		rdd = -xv*eps*N*N*cos(N*theta)
		return r,rd,rdd



class Layer:
	def __init__(self,shape,m):
		self.m = m
		self.shape = shape.copy()
	def __str__(self):
		return str(self.shape)+", m="+str(self.m)
	def copy(self):
		return Layer(self.shape.copy(),self.m)
	m=None
	shape=None

class Particle(object):
	def __str__(self):
		Nlays = len(self.layers)
		if Nlays==1:
			s = "homogeneous "+str(self.layers[0])+"\n"
		else:
			s = ""
			for l in xrange(Nlays):
			  s +="layer "+str(l)+": "+str(self.layers[l])+"\n"
		return s
	def copy(self):
		kwargs = self.copy_args
		del kwargs['self']
		return self.__class__(**kwargs)
	layers = []
	def get_xv(self):
		return self.layers[0].shape.xv
	xv = property(fget=get_xv)
	def get_g(self):
		return pi*self.xv**2
	g  = property(fget=get_g)

	def get_Nlayers(self):
		return len(self.layers)
	Nlayers = property(fget=get_Nlayers)


class SpheroidHom(Particle):
	def __init__(self,ab,xv,m):
		self.copy_args = locals()
		self.layers = [Layer(ShapeSpheroid(ab,xv),m)]

class ChebParticleHom(Particle):
	def __init__(self,N,eps,xv,m):
		self.copy_args = locals()
		self.layers = [Layer(ShapeChebParticle(N,eps,xv),m)]

class Layered_EqShape_Particle(Particle):
	def __init__(self,shape0,ms,volumes,Nlayers):
		self.copy_args = locals()
		# Update Nlayers to match series
		series_len = len(ms)
		Nseries = int(ceil(Nlayers/(series_len+0.)))
		Nlayers1 = Nseries*series_len
		if not Nlayers==Nlayers1:
			print "Updating number of layers"\
				+" from "+str(Nlayers)\
				+" to "  +str(Nlayers1)+"\n"\
				+" to match "+str(Nseries)+" series"\
				+" of "+str(series_len)+" layers."
			Nlayers = Nlayers1

		self.layers =[]
		l = 0
		# derrivation of xv_k for EqVolume
		# R_k - equivolume sphere radius
		# sphere volume is 4/3*pi*R**3
		# (R_1)**3 - (R_2)**3 =...= (R_{n-1})**3-(R_n)**3 = (R_n)**3
		# (R_{n-k})**3 = (k+1)(R_n)**3
		# (R_1)**3 = n(R_n)**3 => (R_n)**3 = (R_1)**3/n
		# (R_k) = ((n-k)/n)**(1/3) R_1
		# xv_k = 2pi*R_k/lambda
		# (xv_k) = ((n-k)n)**(1/3) xv_1
		coefs = zeros(series_len)
		volumes = volumes/(sum(volumes)+0.) #normalize volume fract-s
		self.ms = ms
		self.volumes = volumes*100
		for i in xrange(series_len):
			coefs[i]=sum(volumes[i:])

		serie_no = 0
		while len(self.layers)<Nlayers:
			serie_no = serie_no+1
			for i in xrange(series_len):
				m = ms[i]
				l = l+1
				shape = shape0.copy()
				if l>1:
				   shape.xv *= ((coefs[i]+Nseries-serie_no)\
				   		 /Nseries)**(1/3.)
				self.layers.append(Layer(shape,m))
				if l==Nlayers: break
	def __str__(self):
		ms_str = "[ "
		for m in self.ms:
			ms_str = ms_str + str(m) +"; "
		ms_str = ms_str[:-2] + " ]"

		vols_str = "[ "
		for v in self.volumes:
			vols_str = vols_str + str(v) +"% ; "
		vols_str = vols_str[:-2] + " ]"

		return str(self.Nlayers)+"-layered "\
				+str(self.layers[0].shape)+"\n"\
				+" with layer series "+ms_str+"\n"\
				+" with volume fractions "+vols_str

class Layered_EqShape_EqVolume_Particle(Layered_EqShape_Particle):
	def __init__(self,shape0,ms,Nlayers):
		self.copy_args = locals()
		volumes = zeros(len(ms))+1
		super(Layered_EqShape_EqVolume_Particle,self)\
				.__init__(shape0,ms,volumes,Nlayers)

class EffMedium_Particle(Particle):
	mixing_rule_str = "Average mixing rule"
	def __init__(self,shape0,ms,volumes):
		self.copy_args = locals()
		self.ms = ms
		self.volumes = volumes/(sum(volumes)+0.)*100
		m_mixed = self.mixing_rule(ms,volumes)
		self.layers = [Layer(shape0,m_mixed)]
	def __str__(self):
		ms_str = "[ "
		for m in self.ms:
			ms_str = ms_str + str(m) +"; "
		ms_str = ms_str[:-2] + " ]"

		vols_str = "[ "
		for v in self.volumes:
			vols_str = vols_str + str(v) +"% ; "
		vols_str = vols_str[:-2] + " ]"

		return "Effective medium ("+self.mixing_rule_str+")\n "\
				+str(self.layers[0].shape)+"\n"\
				+" with layer series "+ms_str+"\n"\
				+" with volume fractions "+vols_str
	def mixing_rule(self,ms,vols):
		return sum(array(ms)*array(vols))/sum(array(vols))


class EffMedium_MaxwellGarnett_Particle(EffMedium_Particle):
	mixing_rule_str = "Maxwell-Garnett"
	def mixing_rule(self,ms,vols):
		e1 = ms[0]**2
		e2 = ms[1]**2
		p = vols[0]/(sum(vols)+0.)
		ee = e2*(1 +(3*p*(e1-e2)/(e1+2*e2))/(1-p*(e1-e2)/(e1+2.*e2)))
		return sqrt(ee)

class EffMedium_InvMaxwellGarnett_Particle(EffMedium_Particle):
	mixing_rule_str = "Inverse Maxwell-Garnett"
	def mixing_rule(self,ms,vols):
		e1 = ms[1]**2
		e2 = ms[0]**2
		p = vols[1]/(sum(vols)+0.)
		ee = e2*(1 +(3*p*(e1-e2)/(e1+2*e2))/(1-p*(e1-e2)/(e1+2.*e2)))
		return sqrt(ee)

class EffMedium_Bruggeman_Particle(EffMedium_Particle):
	mixing_rule_str = "Bruggeman"
	def mixing_rule(self,ms,vols):
		e1 = ms[0]**2
		e2 = ms[1]**2
		p1 = vols[0]/(sum(vols)+0.)
		p2 = vols[1]/(sum(vols)+0.)

		a = 2*(p1+p2)
		b = (p2-2*p1)*e1 + (p1-2*p2)*e2
		c = -(p1+p2)*e1*e2
		ee = (-b+sqrt(b**2-4*a*c))/(2*a)
		return sqrt(ee)

class Boundary:
	def __init__(self,lab,l):
		lay1 = lab.particle.layers[l]
		shape = lay1.shape
		m2 = lay1.m
		if l==0:
			m1 = lab.m1
		else:
			m1 = lab.particle.layers[l-1].m
		is_last = (l==len(lab.particle.layers)-1)
		is_first= (l==0)
		self.layer_no = l+1
		self.parentLab = lab
		self.shape = shape
		self.k0 = lab.k0
		self.m1 = m1
		self.m2 = m2
		self.k1 = self.k0*m1
		self.k2 = self.k0*m2
		self.e1 = m1**2
		self.e2 = m2**2
		self.e12 = self.e1/self.e2
		self.e21 = self.e2/self.e1
		self.is_last = is_last
		self.is_first= is_first
	def __str__(self):
		return "m1="+str(self.m1)+"; m2="+str(self.m2)\
		      +"; "+str(self.shape)
	layer_no = None
	parentLab = None
	m1 = None
	m2 = None
	k1 = None
	k2 = None
	k0 = None
	e1 = None
	e2 = None
	e12= None
	e21= None
	shape = None
	is_last = None
	is_first=None

class Lab:
	def __init__(self,particle,alpha,m1=1):
		self.particle = particle.copy()
		self.alpha = alpha
		self.m1 = m1
		self.k1 = self.k0*m1

	def __str__(self):
		return "m1="+str(self.m1)\
				+", alpha="+str(self.alpha*180/pi)+"\n"\
				+str(self.particle)
	
	particle=None
	alpha = None
	k0 = 1
	m1 = None
	k1 = None
	Pna= None
	def boundary(self,l=0):
		return Boundary(self,l)

	def boundaries(self):
		for l in xrange(len(self.particle.layers)):
			yield self.boundary(l)

	def boundaries_reversed(self):
		for l in reversed(xrange(len(self.particle.layers))):
			yield self.boundary(l)

	def get_inc(self,n,m,axisymm=False):
		#if (self.Pna is None) or (len(self.Pna)!=n):
		self.Pna,self.Pdna = lpmn(n,n,cos(self.alpha))[:]
		self.sina = sin(self.alpha)
		sina = self.sina
		Pna = self.Pna[m,m:]
		l = arange(m,n+1)
		if axisymm:
			c_inc = 1j**l * (2*l+1.)/(l*(l+1.)) * Pna
		else:
			Pdna= self.Pdna[m,m:]
			c_inc = zeros(2*(n-m+1),dtype=complex)
			if sina==0. and m==1:
				c_inc[:n-m+1] = -1j**(l-1)*(2*l+1)
			else:
				c_inc[:n-m+1] =  -1j**(l-1)/sina \
				          * 2*(2*l+1)\
				          *factorial(l-m)/factorial(l+m) \
				          * Pna
		return c_inc

	def get_Cext(self,c_sca):
		Cext = self.get_Cext_m(1,c_sca[0],True)
		for m in xrange(1,len(c_sca)):
			Cext += self.get_Cext_m(m,c_sca[m],False)

		g = self.particle.g
		Qext = Cext/g
		return Cext,Qext

	def get_Cext_m(self,m,c_sca,axisymm=False):
		sina = self.sina
		k1 = self.k1
		if axisymm:
			n=len(c_sca)
			Pna = self.Pna[1,1:n+1]
			l = arange(1,n+1)
			Cext = sum((1j**(-l)*c_sca).real*Pna)
		else:
			n=len(c_sca)/2+m-1
			Pna  = self.Pna[m,m:n+1]
			Pdna = self.Pdna[m,m:n+1]
			a = c_sca[:n-m+1]
			b = c_sca[n-m+1:]
			l = arange(m,n+1)
			if self.alpha==0 and m==1:
				Cext = -sum((1j**(-l)*b*l*(l+1)/2.).real)
			else:
				Cext =-sum( ( \
				      1j**(-l+1)*(k1*a*Pna+1j*b*Pdna)*sina \
				      ).real )
		
		Cext = -4*pi/k1**2 * Cext
		return Cext

	def get_Csca(self,c_sca):
		Csca = self.get_Csca_m(1,c_sca[0],True)
		for m in xrange(1,len(c_sca)):
			Csca += self.get_Csca_m(m,c_sca[m],False)

		g = self.particle.g
		Qsca = Csca/g
		return Csca,Qsca

	def get_Csca_m(self,m,c_sca,axisymm=False):
		k1 = self.k1
		if axisymm:
			n=len(c_sca)
			l = arange(1,n+1)
			Csca = 4*sum(abs(c_sca)**2 * l*(l+1)/(2*l+1.))
		else:
			n=len(c_sca)/2+m-1
			a = c_sca[:n-m+1]
			b = c_sca[n-m+1:]
			l = arange(m,n+1)
			om_ln,kap_ln,kap_nl,tau_ln = get_omegakappatau(n,m)
			Csca = sum((
			       1j**fromfunction(lambda l,n:n-l,shape(om_ln))
			      *( outer( k1**2*a,conj(a))* om_ln	\
				+outer( 1j*k1*b,conj(a))*kap_ln	\
				-outer( 1j*k1*a,conj(b))*kap_nl	\
				+outer(       b,conj(b))*tau_ln)).real )
		Csca = pi/k1**2 * Csca
		return Csca

	def get_amplitude_matrix(self,c_sca_tm,c_sca_te,Theta,phi):
		n_max = max(len(c_sca_tm[0]),len(c_sca_te[0]))
		m_max = max(len(c_sca_tm),len(c_sca_te))
		A2,A4 = self.__get_amplitude_matrix2(c_sca_tm,Theta,phi)
		A1,A3 = self.__get_amplitude_matrix2(c_sca_te,Theta,phi)
		return array([[-A2,A3],[A4,A1]])

	def __get_amplitude_matrix2(self,c_sca,Theta,phi):
		#debug_here()
		n = len(c_sca[0])
		ms = len(c_sca)
		P = core.get_Pmn(ms-1,n,cos(self.alpha+Theta))
		PnT = P[:,0,1,1:]
		l = arange(1,n+1)
		A1 = sum( 1j**(-l)*c_sca[0]*PnT, axis=1) #axes: 0-Theta, 1-l
		A2 = zeros_like(A1)
		sinT = sin(self.alpha+Theta)
		cosT = cos(self.alpha+Theta)
		for m in xrange(1,ms):
			PnT = P[:,0,m,m:]
			PdnT= P[:,1,m,m:]
			l = arange(m,n+1)
			a = c_sca[m][:n-m+1]
			b = c_sca[m][n-m+1:]
			A1 -= sinT*cos(m*phi)\
			     *sum(1j**(-l+1)*(a*PnT+1j*b*PdnT), axis=1)
			#A2 += m*sin(m*phi)/sinT\
			#     *sum(1j**(-l)*b*Pnt, axis=1)
			A2 += m*sin(m*phi)*cosT\
			     *sum(((l-m)*1j**(-l))*b*P[:,0,m-1,m:], axis=1)
			A2 += m*sin(m*phi)\
			     *sum(((l+m)*1j**(-l))*b*P[:,0,m-1,m-1:size(P,3)-1], axis=1)
		# PnT/sinT = null when SinT==0
		# to solve this we can use:
		#	sqrt(1-x**2)P^m_l(x) 
		#		= (l-m)xP^{m-1)_l(x) - (l+m)P^{m-1}_{l-1}(x)
		# Here x will be cos(theta), sqrt(1-x**2)=cos(theta)
		return A1,A2

	def get_int_plr(self,A):
		F11 = ( abs(A[0,0])**2 + abs(A[0,1])**2 \
		       +abs(A[1,0])**2 + abs(A[1,1])**2)/2.
	
		F12 = ( abs(A[0,0])**2 - abs(A[0,1])**2 \
		       +abs(A[1,0])**2 - abs(A[1,1])**2)/2.

		g = self.particle.g
		intensity = F11/g
		linpol_degree = -F12/F11

		intensity = where(intensity==inf,nan,intensity)
		intensity = where(intensity==0,  nan,intensity)
		linpol_degree = where(linpol_degree==inf,nan,linpol_degree)
		linpol_degree = where(linpol_degree==0,  nan,linpol_degree)
		return intensity,linpol_degree

	def get_vector_field_sca(self,c_sca,X,Y,Z):
		vX = zeros_like(X)+.0
		vY = zeros_like(X)+.0
		vZ = zeros_like(X)+.0
		for i in xrange(size(X,0)):
		  for j in xrange(size(Y,1)):
		     for k in xrange(size(Z,2)):
			     cPoint = array([X[i,j,k],Y[i,j,k],Z[i,j,k]])
			     sPoint = core.coord_cartesian2spherical(cPoint)
			     r,t,p = list(sPoint)
			     if r >= self.particle.layers[0].shape.R(t)[0]:
			         vec = core.get_vector(cPoint,c_sca,self.k1)
				 if abs(vec[0])<1e-16:vec[0]=.0
				 if abs(vec[1])<1e-16:vec[1]=.0
				 if abs(vec[2])<1e-16:vec[2]=.0
			         vX[i,j,k],vY[i,j,k],vZ[i,j,k] = list(vec)
		return vX,vY,vZ




def set_cheb_particle(n,eps,xv):
	r = xv*(1+eps*cos(n*knots))
	rd = -n*xv*eps*sin(n*knots)
	set_r(xv,r,rd)
	return r,rd

def get_omegakappatau(n,m):
	dim = n-m+1
	delta11 = identity(dim)
	if dim-1>0:
		delta21 = diag(ones(dim-1),-1)
		delta12 = diag(ones(dim-1), 1)
	else:
		delta21 = 0
		delta12 = 0
	if dim-2>0:
		delta31 = diag(ones(dim-2),-2)
		delta13 = diag(ones(dim-2), 2)
	else:
		delta31 = 0
		delta13 = 0
	l = arange(m,n+1)
	ff = factorial(l+m)/factorial(l-m) *2/(2*l+1.)
	ff = transpose(repeat([ff],len(l),0))
	l  = transpose(repeat([ l],len(l),0))

	omega = ( 2*(l*(l+1)+m**2-1)/((2*l+3)*(2*l-1.)) * delta11 \
		 -(l-m)*(l-m-1)     /((2*l-1)*(2*l-3.)) * delta31    \
	         -(l+m+1)*(l+m+2)   /((2*l+3)*(2*l+5.)) * delta13 )  \
	       * ff
	kappa = ( (l+1)*(l-m)/(2*l-1.) * delta21    \
		 -l*(l+m+1)/(2*l+3.)   * delta12 )  \
	       * ff
	tau   = l*(l+1) * ff * delta11
	return omega,kappa,transpose(kappa),tau
