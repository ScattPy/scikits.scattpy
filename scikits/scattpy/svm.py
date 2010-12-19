from numpy import *
import core
import laboratory
from scipy import sparse
#from scipy.sparse.linalg.dsolve import linsolve
#import pdb
#from IPython.Debugger import Tracer; debug_here = Tracer()
#from IPython.Shell import IPShellEmbed; ishell = IPShellEmbed()
#import pylab
#import plotting

ab = 2.#sqrt(2.)
#xv = 2.
#m2 = 1.5+0.0j
#alpha = pi/4
#delta_c=1e-16
#core.set_ngauss(200)
#nrange = range(4,80,2)

eps = 0.21
xv = 1.
m2 = 1.5+0.0j
alpha = 10*pi/180.
delta_c=1e-16
ngauss_coef = 6
Nrange = range(10,40,2)
Ms = None

a = ab*xv/ab**(1./3.)
b = a*ab

RESULTS = None

class Results:
	def __init__(self,c_sca_te,delta_te,N_te,convlog_te,\
	                  c_sca_tm,delta_tm,N_tm,convlog_tm):
		self.c_sca_te = c_sca_te
		self.delta_te = delta_te
		self.N_te = N_te
		self.convlog_te = convlog_te

		self.c_sca_tm = c_sca_tm
		self.delta_tm = delta_tm
		self.N_tm = N_tm
		self.convlog_tm = convlog_tm

def svm(lab,nrange,accuracyLimit=None,ngauss="auto",conv_stop=True,conv_test=True,iterative=True):
	print "\n","*"*60
	print lab
	print "*"*60

	if conv_test:
		nrange = svm_test_M0(lab,nrange,accuracyLimit,ngauss,conv_stop,iterative)

	convlog_te = {'N':[],'delta':[]}
	convlog_tm = {'N':[],'delta':[]}

	# first we need at least one run of solver to have initial values
	global Ms
	Cext_tm,c_sca_tm, Cext_te,c_sca_te = svm_n(lab,nrange[0],ngauss,Ms,iterative)

	delta_tm = delta_te = 1e16 # vars for best accuracies obtained
	Cext_tm_p = Cext_tm  # value of Cext_tm at previous iteration
	Cext_te_p = Cext_te  # value of Cext_te at previous iteration
	N_te = N_tm = nrange[0]

	errors_tm = False
	errors_te = False

	print "N \t\t err TM\t\t err TE"
	for ni,n in enumerate(nrange[1:]):
		try:
		  Cext_tm_n,c_sca_tm_n, Cext_te_n,c_sca_te_n \
				  = svm_n(lab,n,ngauss,Ms,iterative)
		except linalg.linalg.LinAlgError, e:
			print "Terminating: '"+str(e)+"' error received"
			break
		delta_tm_n = abs(Cext_tm_n-Cext_tm_p)/Cext_tm_n
		delta_te_n = abs(Cext_te_n-Cext_te_p)/Cext_te_n

		convlog_tm['N'].append(n)
		convlog_tm['delta'].append(delta_tm_n)
		convlog_te['N'].append(n)
		convlog_te['delta'].append(delta_te_n)
		print "%(n)s\t\t%(delta_tm_n) 5.1e\t%(delta_te_n) 5.1e" \
				% locals()

		if Cext_tm_n<0: errors_tm = True
		if Cext_te_n<0: errors_te = True

		if conv_stop and errors_tm and errors_te:
			print "Terminating: no more convergence expected"
			break

		if 0<delta_tm_n<delta_tm:
			delta_tm = delta_tm_n
			c_sca_tm = c_sca_tm_n
			N_tm = n

		if 0<delta_te_n<delta_te:
			delta_te = delta_te_n
			c_sca_te = c_sca_te_n
			N_te = n

		if conv_stop and\
		   delta_tm<1e-2 and delta_tm_n/delta_tm > 1e2 and\
		   delta_te<1e-2 and delta_te_n/delta_te > 1e2:
			print "Terminating: no more convergence expected"
			break

		if not accuracyLimit==None and\
		   delta_tm>=0 and delta_te>=0 and\
		   delta_tm+delta_te<accuracyLimit:
			print "Terminating: accuracy limit obtained"
			break

		Cext_tm_p = Cext_tm_n
		Cext_te_p = Cext_te_n
		if conv_stop and\
		   (Cext_tm_n<0 or Cext_te_n<0) and\
		   delta_tm+delta_te<2e-2:
			print "Terminating: no more convergence expected"
			break

		if conv_stop and\
		   ni>10 and (delta_tm+delta_te>1e-2):
			print "Terminating: no more convergence expected"
			break

	
	# Print results
	Cext_tm,Qext_tm = lab.get_Cext(c_sca_tm)
	Cext_te,Qext_te = lab.get_Cext(c_sca_te)
	Csca_tm,Qsca_tm = lab.get_Csca(c_sca_tm)
	Csca_te,Qsca_te = lab.get_Csca(c_sca_te)
	n_tm = len(c_sca_tm[0])
	n_te = len(c_sca_te[0])

	print "TM mode: Qext=%(Qext_tm)20.16e , Qsca=%(Qsca_tm)20.16e, delta=%(delta_tm)5.1e , n=%(n_tm)s" % locals()

	print "TE mode: Qext=%(Qext_te)20.16e , Qsca=%(Qsca_te)20.16e, delta=%(delta_te)5.1e , n=%(n_te)s" % locals()

#	Theta_d = arange(0,180,10)
#	Theta = Theta_d*pi/180
#	#TODO: fix with phi != 0
#	phi = 0
#	A = lab.get_amplitude_matrix(c_sca_tm,c_sca_te,Theta,phi)
#	intensity,linpol_degree = lab.get_int_plr(A)
#	pylab.semilogy(Theta_d,intensity)
#	pylab.show()
#	pylab.plot(Theta_d,linpol_degree)
#	pylab.show()

	#print "Theta\tF11\t\tF12"
	#for k in xrange(len(Theta)):
	#	print Theta_d[k],"\t",intensity[k],"\t",linpol_degree[k]
	#print "\nF11\n"
	#print mat([Theta_d,intensity]).T
	#print "\nF21/F11\n"
	#print mat([Theta_d,linpol_degree]).T
	#debug_here()
	global RESULTS
	RESULTS = Results(c_sca_te,delta_te,N_te,convlog_te,\
	                  c_sca_tm,delta_tm,N_tm,convlog_tm)
	return [c_sca_tm,delta_tm],[c_sca_te,delta_te]

def svm_n(lab,n,ngauss,ms=None,iterative=True):
	global th
	if ngauss=="auto":
		core.set_ngauss(ngauss_coef*n)
	else:
		core.set_ngauss(ngauss)
	Pn = core.get_Pmn(n,n,core.cost)
	#Pn = core.get_Pmn_normed(n,n,core.cost)

	JnHn1l = [] # Bessel and Hankel functions for all layers
	JnHn2l = []
	for bnd in lab.boundaries():
		Rs = bnd.shape.R(core.knots) #here Rs = [r,rd,rdd]
		core.set_consts(Rs)
		JnHn1l.append( core.get_JnHn(n,bnd.k1*core.r) )
		JnHn2l.append( core.get_JnHn(n,bnd.k2*core.r) )

	Cext_tm = Cext_te = 0
	c_sca_tm = []
	c_sca_te = []
	if lab.alpha==0:
		ms = [1]
		c_sca_tm.append(zeros(n))
		c_sca_te.append(zeros(n))
	else:
		if not ms:
			ms = xrange(n+1)

	for m in ms:
		# Instead of m=0 we start with m=-1 (axisymmetric part)
		if m==0:
			m=1
			axisymm = True
		else:
			axisymm = False
	
		c_inc = lab.get_inc(n,m,axisymm)

		Nsize = len(c_inc)
		Nlays = lab.particle.Nlayers

		if not iterative:
			AAtm = sparse.lil_matrix((Nlays*Nsize*2,Nlays*Nsize*2),dtype=complex)
			AAte = sparse.lil_matrix((Nlays*Nsize*2,Nlays*Nsize*2),dtype=complex)
			BB   = zeros((Nlays*Nsize*2,1),dtype=complex)

		for bnd in lab.boundaries_reversed():
			Rs = bnd.shape.R(core.knots) #here Rs = [r,rd,rdd]
			core.set_consts(Rs)
			Jn1,Hn1 = JnHn1l[bnd.layer_no-1]
			Jn2,Hn2 = JnHn2l[bnd.layer_no-1]
			#Jn1,Hn1 = core.get_JnHn(n,bnd.k1*core.r)
			#Jn2,Hn2 = core.get_JnHn(n,bnd.k2*core.r)

			K11,K12,K21tm,K21te,K22tm,K22te =\
			   svm_bnd_system(n,m,bnd,Pn,Jn1,Jn2,Hn1,Hn2,axisymm)
			if not iterative:
			   # Single matrix method
			   Nlay = bnd.layer_no
			   Npos1a = (Nlay-1)*Nsize*2
			   Npos1b = Nlay*Nsize*2
			   Npos2a = Nsize+(Nlay-2)*Nsize*2
			   Npos2b = Nsize+(Nlay-1)*Nsize*2
			   Npos2c = Nsize+(Nlay-0)*Nsize*2
			   if bnd.is_last:
				AAtm[Npos1a:,Npos2b:] = K21tm
				AAte[Npos1a:,Npos2b:] = K21te
			   else:
				AAtm[Npos1a:Npos1b,Npos2b:Npos2c]= c_[K21tm,K22tm]
				AAte[Npos1a:Npos1b,Npos2b:Npos2c]= c_[K21te,K22te]

			   if not bnd.is_first:
				AAtm[Npos1a:Npos1b,Npos2a:Npos2b]= -c_[K11,K12]
				AAte[Npos1a:Npos1b,Npos2a:Npos2b]= -c_[K11,K12]
			   else:
				AAtm[:Nsize*2,:Nsize] = -K12
				AAte[:Nsize*2,:Nsize] = -K12
				BB[:Nsize*2,:] = K11*mat(c_inc).T
				#c_sca_tm_m = sparse.linalg.spsolve(AAtm.tocsr(),BB)[:Nsize]
				#c_sca_te_m = sparse.linalg.spsolve(AAte.tocsr(),BB)[:Nsize]
				c_sca_tm_m = sparse.linalg.bicg(AAtm.tocsr(),BB,tol=1e-10)[:Nsize]
				c_sca_te_m = sparse.linalg.bicg(AAte.tocsr(),BB,tol=1e-10)[:Nsize]
			else:
			   # Iterative method
			   if bnd.is_last:
				P0tm = K21tm
				P0te = K21te
			   else:
				P0tm = c_[K21tm,K22tm]*Ptm
				P0te = c_[K21te,K22te]*Pte
			
			   if not bnd.is_first:
				Ptm = linalg.solve(c_[K11,K12],P0tm)
				Pte = linalg.solve(c_[K11,K12],P0te)
			   else:
				BB = -K11*mat(c_inc).T
				c_sca_tm_m = linalg.solve(c_[K12,P0tm],BB)
				c_sca_te_m = linalg.solve(c_[K12,P0te],BB)

				#cast to array
				c_sca_tm_m = array(c_sca_tm_m)[:len(c_inc),0] 
				c_sca_te_m = array(c_sca_te_m)[:len(c_inc),0]

		CextM_tm_m = lab.get_Cext_m(m,c_sca_tm_m,axisymm)
		CextM_te_m = lab.get_Cext_m(m,c_sca_te_m,axisymm)
		Cext_tm +=CextM_tm_m
		Cext_te +=CextM_te_m
		if CextM_tm_m/Cext_tm < delta_c and \
		   CextM_te_m/Cext_te < delta_c       : break
		c_sca_tm.append(c_sca_tm_m)
		c_sca_te.append(c_sca_te_m)

	return Cext_tm,c_sca_tm, Cext_te,c_sca_te

def svm_bnd_system(n,m,bnd,Pn,Jn1,Jn2,Hn1,Hn2,axisymm=False):
	Aj1 = core.matA(n,m,Jn1,Pn)
	Aj2 = core.matA(n,m,Jn2,Pn)
	Ah1 = core.matA(n,m,Hn1,Pn)
	Bj1 = core.matB(n,m,Jn1,Pn,bnd.k1)
	Bj2 = core.matB(n,m,Jn2,Pn,bnd.k2)
	Bh1 = core.matB(n,m,Hn1,Pn,bnd.k1)

	# Axisymmetric part
	if axisymm:
		Cj2 = core.matC(n,m,Jn2,Pn,bnd.k2,bnd.e12,B=Bj2)

		K11 = bmat([[ Aj1.T ],\
			    [ Bj1.T ] ])

		K12 = bmat([[ Ah1.T ],\
			    [ Bh1.T ] ])

		K21tm = bmat([[ Aj2.T ],\
			      [ Cj2.T ] ])

		K21te = bmat([[ Aj2.T ],\
			      [ Bj2.T ] ])

		# For not-last boundaries we need Ah2,Bh2, etc.
		if not bnd.is_last:
			Ah2 = core.matA(n,m,Hn2,Pn)
			Bh2 = core.matB(n,m,Hn2,Pn,bnd.k2)
			Ch2 = core.matC(n,m,Hn2,Pn,bnd.k2,bnd.e12,B=Bh2)

			K22tm = bmat([[ Ah2.T ],\
				      [ Ch2.T ] ])

			K22te = bmat([[ Ah2.T ],\
				      [ Bh2.T ] ])
		else:
			K22tm=K22te=mat(zeros_like(K11))

	# Non-axisymmetric part
	else:
		Dj2 = core.matD(n,m,Jn2,Pn,bnd.k2,bnd.e12,B=Bj2)
		Ej2 = core.matE(n,m,Jn2,Pn,bnd.k2,bnd.e12)
		Fj2 = core.matF(n,m,Jn2,Pn,bnd.k2,bnd.e12)
		Gj2 = core.matG(n,m,Jn2,Pn,bnd.k2,bnd.e12,B=Bj2)
		
		Aj12 = core.matA12(n,m,Jn2,Pn,bnd.k2,bnd.e21,A=Aj2)
		Aj22 = core.matA22(n,m,Jn2,Pn,bnd.k2,bnd.e21,B=Bj2)
		Aj32 = core.matA32(n,m,Jn2,Pn,bnd.k2,bnd.e21)
		Aj42 = core.matA42(n,m,Jn2,Pn,bnd.k2,bnd.e21)
	
		Aj14 = core.matA14(n,m,Jn2,Pn,bnd.k2,bnd.e21)
		Aj24 = core.matA24(n,m,Jn2,Pn,bnd.k2,bnd.e21)
		Aj34 = core.matA34(n,m,Jn2,Pn,bnd.k2,bnd.e21,A=Aj2)
		Aj44 = core.matA44(n,m,Jn2,Pn,bnd.k2,bnd.e21,B=Bj2)
	
		c0 = mat(zeros_like(Aj2))

		K11 = bmat( [[ Aj1.T,   c0  ],\
			     [ Bj1.T,   c0  ],\
			     [    c0, Aj1.T ],\
			     [    c0, Bj1.T ]] )

		K12 = bmat( [[ Ah1.T,   c0  ],\
			     [ Bh1.T,   c0  ],\
			     [    c0, Ah1.T ],\
			     [    c0, Bh1.T ]] )

		K21tm = bmat( [[ Aj2.T,   c0  ],\
			       [ Dj2.T, Ej2.T ],\
			       [    c0, Aj2.T ],\
			       [ Fj2.T, Gj2.T ]] )

		K21te = bmat( [[ Aj12.T, Aj14.T ],\
			       [ Aj22.T, Aj24.T ],\
			       [ Aj32.T, Aj34.T ],\
			       [ Aj42.T, Aj44.T ]] )

		# For not-last boundaries we need Ah2,Bh2, etc.
		if not bnd.is_last:
			Ah2 = core.matA(n,m,Hn2,Pn)
			Bh2 = core.matB(n,m,Hn2,Pn,bnd.k2)

			Dh2 = core.matD(n,m,Hn2,Pn,bnd.k2,bnd.e12,B=Bh2)
			Eh2 = core.matE(n,m,Hn2,Pn,bnd.k2,bnd.e12)
			Fh2 = core.matF(n,m,Hn2,Pn,bnd.k2,bnd.e12)
			Gh2 = core.matG(n,m,Hn2,Pn,bnd.k2,bnd.e12,B=Bh2)
			
			Ah12 = core.matA12(n,m,Hn2,Pn,bnd.k2,bnd.e21,A=Ah2)
			Ah22 = core.matA22(n,m,Hn2,Pn,bnd.k2,bnd.e21,B=Bh2)
			Ah32 = core.matA32(n,m,Hn2,Pn,bnd.k2,bnd.e21)
			Ah42 = core.matA42(n,m,Hn2,Pn,bnd.k2,bnd.e21)
		
			Ah14 = core.matA14(n,m,Hn2,Pn,bnd.k2,bnd.e21)
			Ah24 = core.matA24(n,m,Hn2,Pn,bnd.k2,bnd.e21)
			Ah34 = core.matA34(n,m,Hn2,Pn,bnd.k2,bnd.e21,A=Ah2)
			Ah44 = core.matA44(n,m,Hn2,Pn,bnd.k2,bnd.e21,B=Bh2)

			K22tm = bmat( [[ Ah2.T,   c0  ],\
				       [ Dh2.T, Eh2.T ],\
				       [    c0, Ah2.T ],\
				       [ Fh2.T, Gh2.T ]] )

			K22te = bmat( [[ Ah12.T, Ah14.T ],\
				       [ Ah22.T, Ah24.T ],\
				       [ Ah32.T, Ah34.T ],\
				       [ Ah42.T, Ah44.T ]] )

		else:
			K22tm = K22te = mat(zeros_like(K11))

	return K11,K12,K21tm,K21te,K22tm,K22te

def get_Bs_Br(bnd,Jn1,Jn2,Hn1,Pn):
	Bs = Br = 0
	k1 = bnd.k1
	k2 = bnd.k2
	e12 = bnd.e12
	for k in xrange(len(core.weights)):
		jn1 = Jn1[k,0,1:]
		jdn1 = Jn1[k,1,1:]
		jn2 = Jn2[k,0,1:]
		jdn2 = Jn2[k,1,1:]
		hn1 = Hn1[k,0,1:]
		hdn1 = Hn1[k,1,1:]
		pn = Pn[k,0,1,1:]
		pdn = Pn[k,1,1,1:]

		r = core.r[k]
		rdr = core.rdr[k]
		sint = core.sint[k]
		ctgt = core.ctgt[k]

		Aj1 = mat(outer( jn1*pn, sint*pn))
		Aj2 = mat(outer( jn2*pn, sint*pn))
		Ah1 = mat(outer( hn1*pn, sint*pn))
		Bj1 = mat(outer(k1*r*jdn1*pn+rdr*sint*jn1*pdn, sint*pn ))
		Bj2 = mat(outer(k2*r*jdn2*pn+rdr*sint*jn2*pdn, sint*pn ))
		Bh1 = mat(outer(k1*r*hdn1*pn+rdr*sint*hn1*pdn, sint*pn ))
		Cj2 = e12*Bj2 + (e12-1.)*(1.-ctgt*rdr)*Aj2

		Bs += 1j*( Bh1*Aj2.T - Ah1*Cj2.T)*core.weights[k]
		Br += 1j*( Bj1*Aj2.T - Aj1*Cj2.T)*core.weights[k]

	return Bs,Br

def svm_test_M0(lab,nrange,accuracyLimit=None,ngauss="auto",conv_stop=True,iterative=True):
	print "\n","*"*60
	print "Convergence test: axisymmetric part"
	print "*"*60

	convlog_te = {'N':[],'delta':[]}
	convlog_tm = {'N':[],'delta':[]}

	# first we need at least one run of solver to have initial values
	Ms = [0]
	Cext_tm,c_sca_tm, Cext_te,c_sca_te = svm_n(lab,nrange[0],ngauss,Ms,iterative)

	delta_tm = delta_te = 1e16 # vars for best accuracies obtained
	Cext_tm_p = Cext_tm  # value of Cext_tm at previous iteration
	Cext_te_p = Cext_te  # value of Cext_te at previous iteration
	nrange_tm = list(nrange)[0:1]
	nrange_te = list(nrange)[0:1]

	print "N \t err TM\t err TE"
	for n in nrange[1:]:
		try:
		  Cext_tm_n,c_sca_tm_n, Cext_te_n,c_sca_te_n \
				  = svm_n(lab,n,ngauss,Ms,iterative)
		except linalg.linalg.LinAlgError, e:
			print "Terminating: '"+str(e)+"' error received"
			break
		delta_tm_n = abs(Cext_tm_n-Cext_tm_p)/Cext_tm_n
		delta_te_n = abs(Cext_te_n-Cext_te_p)/Cext_te_n

		convlog_tm['N'].append(n)
		convlog_tm['delta'].append(delta_tm_n)
		convlog_te['N'].append(n)
		convlog_te['delta'].append(delta_te_n)
		print "%(n)s\t\t%(delta_tm_n) 5.1e\t%(delta_te_n) 5.1e" \
				% locals()

		nrange_tm+=[n]
		nrange_te+=[n]

		if 0<delta_tm_n<delta_tm:
			delta_tm = delta_tm_n
			c_sca_tm = c_sca_tm_n
			del(nrange_tm[:-2])

		if 0<delta_te_n<delta_te:
			delta_te = delta_te_n
			c_sca_te = c_sca_te_n
			del(nrange_te[:-2])

		if conv_stop and\
		   delta_tm<1e-3 and delta_tm_n/delta_tm > 1e2 and\
		   delta_te<1e-3 and delta_te_n/delta_te > 1e2:
			   print "Terminating: no more convergence expected"
			   break

		if not accuracyLimit==None:
			if delta_tm+delta_te<accuracyLimit:
				print "Terminating: accuracy limit obtained"
				break

		Cext_tm_p = Cext_tm_n
		Cext_te_p = Cext_te_n
	
	# Print results
	Cext_tm,Qext_tm = lab.get_Cext(c_sca_tm)
	Cext_te,Qext_te = lab.get_Cext(c_sca_te)
	Csca_tm,Qsca_tm = lab.get_Csca(c_sca_tm)
	Csca_te,Qsca_te = lab.get_Csca(c_sca_te)
	n_tm = len(c_sca_tm[0])
	n_te = len(c_sca_te[0])

	print "TM mode: Qext=%(Qext_tm)20.16e , Qsca=%(Qsca_tm)20.16e, delta=%(delta_tm)5.1e , n=%(n_tm)s" % locals()

	print "TE mode: Qext=%(Qext_te)20.16e , Qsca=%(Qsca_te)20.16e, delta=%(delta_te)5.1e , n=%(n_te)s" % locals()

	print list(unique1d(list(nrange_tm)+list(nrange_te)))
	print nrange_tm,nrange_te
	return list(unique1d(list(nrange_tm)+list(nrange_te)))


def test_svm(nr=xrange(4,50,2),ngauss=80):
	ptcl = laboratory.SpheroidHom(ab,xv,m2)
	lab = laboratory.Lab(ptcl,alpha)
	svm(lab)

def test_cheb(nr=xrange(5,50,5),ngauss=80):
	ptcl = laboratory.ChebParticleHom(5,eps,xv,m2)
	lab = laboratory.Lab(ptcl,alpha)
	svm(lab)

def test_cprof():
	import cProfile
	#cProfile.run('test_svm
