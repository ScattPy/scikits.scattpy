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
		for m,coefs_m in enumerate(coefs):
			if m==0:continue # skip axisymmetric part
			n=len(coefs_m)/2+m-1
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
		        a = coefs_m[:n-m+1]
        		b = coefs_m[n-m+1:]
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

def get_all(particle,xs,ys,zs):
	LAB = scattpy.Lab(particle,0.)
	dx=xs[1]-xs[0]
	dy=ys[1]-ys[0]
	dz=zs[1]-zs[0]
	I = get_I(particle,xs,ys,zs)

	print "solve"
	r = scattpy.svm(LAB,arange(10,40,2),ngauss=200)
	from scikits.scattpy.svm import RESULTS

	print "scattered vector field"
	Rs = get_Ui_Vr(xs,ys,zs,RESULTS.c_sca_tm,LAB.k1,'h')
	Hs = curl(Rs,dx,dy,dz)
	Es = array(curl(Hs,dx,dy,dz))*(-1./1j)
	Hs = [real(Hs[0])*I,real(Hs[1])*I,real(Hs[2])*I]
	Es = [real(Es[0])*I,real(Es[1])*I,real(Es[2])*I]
	Ps = vector_cross(Es,Hs)

	print "incident vector field"
	c_inc = [[0],LAB.get_inc(len(RESULTS.c_sca_tm[1])/2,1)]
	Ri = get_Ui_Vr(xs,ys,zs,c_inc,LAB.k1,'j')
	Hi = curl(Ri,dx,dy,dz)
	Ei = array(curl(Hi,dx,dy,dz))*(-1./1j)
	Hi = [real(Hi[0])*I,real(Hi[1])*I,real(Hi[2])*I]
	Ei = [real(Ei[0])*I,real(Ei[1])*I,real(Ei[2])*I]
	Pi = vector_cross(Ei,Hi)

	return Hs,Es,Ps,Hi,Ei,Pi

def plot_vf(vf,xs,ys,zs,str):
	Fx,Fy,Fz = vf
	n = (len(xs)-1)/2
	if str == 'xz':
		X,Z = meshgrid(xs,zs)
		pylab.pcolor(X,Z,sqrt(Fx**2+Fy**2+Fz**2)[:,n,:].T)
		pylab.colorbar()
		pylab.show()
#		pylab.contourf(X,Z,sqrt(Fx**2+Fy**2+Fz**2)[:,n,:])
#		pylab.quiver(X,Z,(Fx)[:,n,:],(Fz)[:,n,:])
	elif str == 'yz':
		Y,Z = meshgrid(ys,zs)
		pylab.pcolor(Y,Z,sqrt(Fx**2+Fy**2+Fz**2)[n,:,:].T)
		pylab.colorbar()
		pylab.show()
#		pylab.contourf(Y,Z,sqrt(Fx**2+Fy**2+Fz**2)[n,:,:])
#		pylab.quiver(Y,Z,(Fy)[n,:,:],(Fz)[n,:,:])
	elif str == 'xy':
		X,Y = meshgrid(xs,ys)
		pylab.pcolor(X,Y,sqrt(Fx**2+Fy**2+Fz**2)[:,:,n].T)
		pylab.colorbar()
		pylab.show()
#		pylab.contourf(X,Y,sqrt(Fx**2+Fy**2+Fz**2)[:,:,n])
#		pylab.quiver(X,Y,(Fx)[:,:,n],(Fy)[:,:,n])

def plot_vf_lic(F,xs,ys,zs,str,nsize=401,cmap='hot'):
	"""Plot vector field using colored linear integral convolution"""
	from scikits import vectorplot as vp
	from scipy import interpolate
	Fx,Fy,Fz = F
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
	pylab.figimage(sqrt(nF1**2+nF2**2)*(image+imax)/imax/2,cmap=cmap)
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
