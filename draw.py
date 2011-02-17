from numpy import *
from scikits import scattpy
from scikits.scattpy import f_coords
import pylab
from scipy.special import sph_jn,sph_jnyn,lpmn

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

def get_JnHn(n,x):
	JnYn = array([sph_jnyn(n,xl) for xl in x])
	return JnYn[:,:2,:],JnYn[:,:2,:]+1j*JnYn[:,2:,:]

def get_Ui_Vr(xs,ys,zs,coefs,ki,jh):
	Rx = zeros([size(xs),size(ys),size(zs)],complex)
	Ry = zeros([size(xs),size(ys),size(zs)],complex)
	Rz = zeros([size(xs),size(ys),size(zs)],complex)
	
	for i,x in enumerate(xs):
	  print x
	  for j,y in enumerate(ys):
	    for k,z in enumerate(zs):
		U=V=0.
		cPoint = [x,y,z]
		sPoint = f_coords.point_c2s(cPoint)
		r,t,p = sPoint
		cost = cos(t)
		n  = shape(coefs)[2]
		for m,coefs_m in enumerate(coefs):
			if m==0:continue # skip axisymmetric part
		        P,Pd = get_Pmn(m,n,[cost])[0]
        		Pml  = P[m,m:n+1]
		        Pdml = Pd[m,m:n+1]
        		Bess,Hank = get_JnHn(n,[ki*r])
			if jh=='j':
				Rad = Bess[0,0,m:]
				Radd= Bess[0,1,m:]
			else:
	        		Rad = Hank[0,0,m:]
        			Radd= Hank[0,1,m:]
			Nsize = n-m+1
		        a = coefs_m[0,-Nsize:]
        		b = coefs_m[1,-Nsize:]
	        	l = arange(m,n+1)

			U += sum(a*Rad*Pml)*cos(m*p)
			V += sum(b*Rad*Pml)*cos(m*p)

		sVec = (V,0.,0.)
		cVec = f_coords.vector_s2c(sVec,sPoint)
		cVec[2] += U
		Rx[i,j,k],Ry[i,j,k],Rz[i,j,k] = cVec
	return Rx,Ry,Rz

def get_I(particle,xs,ys,zs):
	I = zeros([len(xs),len(ys),len(zs)],int)

	for i,x in enumerate(xs):
	  print x
	  for j,y in enumerate(ys):
	    for k,z in enumerate(zs):
			cPoint = [x,y,z]
			sPoint = f_coords.point_c2s(cPoint)
			r,t,p = sPoint
			R,Rd,Rdd = particle.layers[0].shape.R(t)
			if r>=R:
				I[i,j,k]=1
	return I

def get_If(particle,xs,ys,zs):
	def fR(t):
		return particle.layers[0].shape.R(t)[0]
	return array(f_coords.get_i(fR,xs,ys,zs))

def get_all(particle,xs,ys,zs,res=None):
	global Hmax,Hmin,Emax,Emin,Pmax,Pmin
	LAB = scattpy.Lab(particle,0.)
	if res:
		RESULTS = res
	else:
		print "solve"
		RESULTS = scattpy.svm(LAB,allfields=True)

	print "indicator function"
	dx=xs[1]-xs[0]
	dy=ys[1]-ys[0]
	dz=zs[1]-zs[0]
	I = get_I(particle,xs,ys,zs)

	I1 = where(I==0,1,0)

	print "scattered vector field"
	Rsca = get_Ui_Vr(xs,ys,zs,RESULTS.c_sca_tm,LAB.k1,'h')
	Hsca = curl(Rsca,dx,dy,dz)
	Esca = array(curl(Hsca,dx,dy,dz))*(-1./1j)
	Hsca = [real(Hsca[0])*I,real(Hsca[1])*I,real(Hsca[2])*I]
	Esca = [real(Esca[0])*I,real(Esca[1])*I,real(Esca[2])*I]
	Psca = vector_cross(Esca,Hsca)

	print "incident vector field"
	Rinc = get_Ui_Vr(xs,ys,zs,RESULTS._c_all_tm[0,0],LAB.k1,'j')
	Hinc = curl(Rinc,dx,dy,dz)
	Einc = array(curl(Hinc,dx,dy,dz))*(-1./1j)
	Hinc = [real(Hinc[0])*I,real(Hinc[1])*I,real(Hinc[2])*I]
	Einc = [real(Einc[0])*I,real(Einc[1])*I,real(Einc[2])*I]
	Pinc = vector_cross(Einc,Hinc)

	print "internal vector field"
	Rint = get_Ui_Vr(xs,ys,zs,RESULTS._c_all_tm[1,0],LAB.boundary(0).k2,'j')
	Hint = curl(Rint,dx,dy,dz)
	Eint = array(curl(Hint,dx,dy,dz))*(-1./1j)
	Hint = [real(Hint[0])*I1,real(Hint[1])*I1,real(Hint[2])*I1]
	Eint = [real(Eint[0])*I1,real(Eint[1])*I1,real(Eint[2])*I1]
	Pint = vector_cross(Eint,Hint)

	Emax = max(sqrt(Einc[0]**2+Einc[1]**2+Einc[2]**2).reshape(len(xs)*len(ys)*len(zs)))*2
	Hmax = max(sqrt(Hinc[0]**2+Hinc[1]**2+Hinc[2]**2).reshape(len(xs)*len(ys)*len(zs)))*2
	Pmax = max(sqrt(Pinc[0]**2+Pinc[1]**2+Pinc[2]**2).reshape(len(xs)*len(ys)*len(zs)))*2
	Hmin=Emin=Pmin = 0

	return Hsca,Esca,Psca,Hinc,Einc,Pinc,Hint,Eint,Pint

def plot_vf(vf,xs,ys,zs,str,quiver_step,Fmin=None,Fmax=None):
	Fx,Fy,Fz = vf
	n = (len(xs)-1)/2
	qs = quiver_step
	Fabs = sqrt(Fx**2+Fy**2+Fz**2)
	if str == 'xz':
		X,Z = meshgrid(xs,zs)
		pylab.pcolor(X,Z,log10(Fabs)[:,n,:].T,vmin=Fmin,vmax=Fmax)
		pylab.colorbar()
		pylab.show()
#		pylab.contourf(X,Z,sqrt(Fx**2+Fy**2+Fz**2)[:,n,:])
		pylab.quiver(X[::qs,::qs],Z[::qs,::qs],\
				(Fx/Fabs)[::qs,n,::qs],(Fz/Fabs)[::qs,n,::qs],color='w')
	elif str == 'yz':
		Y,Z = meshgrid(ys,zs)
		pylab.pcolor(Y[,Z,log10(Fabs)[n,:,:].T,vmin=Fmin,vmax=Fmax)
		pylab.colorbar()
		pylab.show()
#		pylab.contourf(Y,Z,sqrt(Fx**2+Fy**2+Fz**2)[n,:,:])
		pylab.quiver(Y[::qs,::qs],Z[::qs,::qs],\
				(Fy/Fabs)[n,::qs,::qs],(Fz/Fabs)[n,::qs,::qs])
	elif str == 'xy':
		X,Y = meshgrid(xs,ys)
		pylab.pcolor(X,Y,log10(Fabs)[:,:,n].T,vmin=Fmin,vmax=Fmax)
		pylab.colorbar()
		pylab.show()
#		pylab.contourf(X,Y,sqrt(Fx**2+Fy**2+Fz**2)[:,:,n])
		pylab.quiver(X[::qs,::qs],Y[::qs,::qs],\
				(Fx/Fabs)[::qs,::qs,n],(Fy/Fabs)[::qs,::qs,n])

def plot_vf_lic(F,xs,ys,zs,str,nsize=401,cmap='hot',ftype=None,Fmin=None,Fmax=None):
	"""Plot vector field using colored linear integral convolution"""
	from scikits import vectorplot as vp
	from scipy import interpolate
	global Hmax,Hmin,Emax,Emin,Pmax,Pmin
	if ftype == "H":
		Fmax = Hmax
		Fmin = Hmin
	elif ftype == "E":
		Fmax = Emax
		Fmin = Emin
	elif ftype == "P":
		Fmax = Pmax
		Fmin = Pmin
	else:
		if not Fmax:
			Fmax = None
		if not Fmin:
			Fmin = None
	n = (len(xs)-1)/2
	if str == 'xz':
		x1=xs
		x2=zs
		F1=F[0][:,n,:]
		F2=F[2][:,n,:]
	elif str == 'yz':
		x1=ys
		x2=zs
		F1=F[1][n,:,:]
		F2=F[2][n,:,:]
	elif str == 'xy':
		x1=xs
		x2=ys
		F1=F[0][:,:,n]
		F2=F[1][:,:,n]
	spl_F1 = interpolate.RectBivariateSpline(x1,x2,F1)
	spl_F2 = interpolate.RectBivariateSpline(x1,x2,F2)
	nx1 = linspace(min(x1),max(x1),nsize)
	nx2 = linspace(min(x2),max(x2),nsize)
	nF1 = spl_F1(nx1,nx2).T
	nF2 = spl_F2(nx1,nx2).T

	kernellen=31
	kernel = sin(arange(kernellen)*pi/kernellen)
	kernel = kernel.astype(float32)
	texture = random.rand(nsize,nsize).astype(float32)
	image = vp.line_integral_convolution(nF1.astype(float32), nF2.astype(float32), texture, kernel)

	imax = max(image.reshape(nsize**2))
	pylab.figimage(log10(sqrt(nF1**2+nF2**2))*(image+imax)/imax/2,cmap=cmap,vmin=Fmin,vmax=Fmax)
	pylab.show()

def curl(F,dx,dy,dz):
	Fx,Fy,Fz = F
	Fxx,Fxy,Fxz = gradient(Fx,dx,dy,dz)
	Fyx,Fyy,Fyz = gradient(Fy,dx,dy,dz)
	Fzx,Fzy,Fzz = gradient(Fz,dx,dy,dz)
	return [ Fzy-Fyz , Fxz-Fzx, Fyx-Fxy ]

def vector_cross(A,B):
	Ax,Ay,Az = A
	Bx,By,Bz = B
	return [ Ay*Bz-Az*By , Az*Bx-Ax*Bz , Ax*By-Ay*Bx ]
