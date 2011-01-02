from numpy import *
#import core
import laboratory
from scipy import sparse
#from scipy.sparse.linalg.dsolve import linsolve
import pdb
#from IPython.Debugger import Tracer; debug_here = Tracer()
#from IPython.Shell import IPShellEmbed; ishell = IPShellEmbed()
#import pylab
#import plotting
import spherical

delta_c=1e-16
ngauss_coef = 6
nrange_default = arange(2,61)
Ms=[]


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


def methods_factory(meth_bnd_system):
   def meth(lab,nrange=None,accuracyLimit=None,ngauss=200,conv_stop=True,conv_test=False,iterative=True):
	print "\n","*"*60
	print lab
	print "*"*60

	if conv_test:
		nrange = svm_test_M0(lab,nrange,accuracyLimit,ngauss,conv_stop,iterative)

	convlog_te = {'N':[],'delta':[]}
	convlog_tm = {'N':[],'delta':[]}

	# get nrange
	global nrange_default
	if nrange==None:
		if lab.particle.layers[0].shape.nrange == None:
			nrange = nrange_default
		else:
			nrange = lab.particle.layers[0].shape.nrange

	# first we need at least one run of solver to have initial values
	global Ms
	Cext_tm,c_sca_tm, Cext_te,c_sca_te = meth_n(meth_bnd_system,lab,nrange[0],ngauss,Ms,iterative)

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
				  = meth_n(meth_bnd_system,lab,n,ngauss,Ms,iterative)
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

	global RESULTS
	RESULTS = Results(c_sca_te,delta_te,N_te,convlog_te,\
	                  c_sca_tm,delta_tm,N_tm,convlog_tm)
	return RESULTS
   return meth

def meth_n(meth_bnd_system,lab,n,ngauss,ms=None,iterative=True):
	global th
	if ngauss=="auto":
		ng=ngauss_coef*n
	else:
		ng=ngauss
	C=spherical.spherical_utilities(ng,n,lab)

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
			K11,K12,K21tm,K21te,K22tm,K22te =\
			   meth_bnd_system(C,m,bnd,axisymm)
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

def svm_bnd_system(C,m,bnd,axisymm=False):
	#pdb.set_trace()
	lay = bnd.layer_no
	C.set_layer_no(lay)
	Aj1 = spherical.matA(C,m,'j',1)
	Aj2 = spherical.matA(C,m,'j',2)
	Ah1 = spherical.matA(C,m,'h',1)
	Bj1 = spherical.matB(C,m,'j',1)
	Bj2 = spherical.matB(C,m,'j',2)
	Bh1 = spherical.matB(C,m,'h',1)

	# Axisymmetric part
	if axisymm:
		Cj2 = spherical.matC(C,m,'j',2,bnd.e12,B=Bj2)

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
			Ah2 = spherical.matA(C,m,'h',2)
			Bh2 = spherical.matB(C,m,'h',2)
			Ch2 = spherical.matC(C,m,'h',2,bnd.e12,B=Bh2)

			K22tm = bmat([[ Ah2.T ],\
				      [ Ch2.T ] ])

			K22te = bmat([[ Ah2.T ],\
				      [ Bh2.T ] ])
		else:
			K22tm=K22te=mat(zeros_like(K11))

	# Non-axisymmetric part
	else:
		Dj2 = spherical.matD(C,m,'j',2,bnd.e12,B=Bj2)
		Ej2 = spherical.matE(C,m,'j',2,bnd.e12)
		Fj2 = spherical.matF(C,m,'j',2,bnd.e12)
		Gj2 = spherical.matG(C,m,'j',2,bnd.e12,B=Bj2)
		
		Aj12 = spherical.matA12(C,m,'j',2,bnd.e21,A=Aj2)
		Aj22 = spherical.matA22(C,m,'j',2,bnd.e21,B=Bj2)
		Aj32 = spherical.matA32(C,m,'j',2,bnd.e21)
		Aj42 = spherical.matA42(C,m,'j',2,bnd.e21)
	
		Aj14 = spherical.matA14(C,m,'j',2,bnd.e21)
		Aj24 = spherical.matA24(C,m,'j',2,bnd.e21)
		Aj34 = spherical.matA34(C,m,'j',2,bnd.e21,A=Aj2)
		Aj44 = spherical.matA44(C,m,'j',2,bnd.e21,B=Bj2)
	
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
			Ah2 = spherical.matA(C,m,'h',2)
			Bh2 = spherical.matB(C,m,'h',2)

			Dh2 = spherical.matD(C,m,'h',2,bnd.e12,B=Bh2)
			Eh2 = spherical.matE(C,m,'h',2,bnd.e12)
			Fh2 = spherical.matF(C,m,'h',2,bnd.e12)
			Gh2 = spherical.matG(C,m,'h',2,bnd.e12,B=Bh2)
			
			Ah12 = spherical.matA12(C,m,'h',2,bnd.e21,A=Ah2)
			Ah22 = spherical.matA22(C,m,'h',2,bnd.e21,B=Bh2)
			Ah32 = spherical.matA32(C,m,'h',2,bnd.e21)
			Ah42 = spherical.matA42(C,m,'h',2,bnd.e21)
		
			Ah14 = spherical.matA14(C,m,'h',2,bnd.e21)
			Ah24 = spherical.matA24(C,m,'h',2,bnd.e21)
			Ah34 = spherical.matA34(C,m,'h',2,bnd.e21,A=Ah2)
			Ah44 = spherical.matA44(C,m,'h',2,bnd.e21,B=Bh2)

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

def ebcm_bnd_system(C,m,bnd,axisymm=False):
	lay = bnd.layer_no
	C.set_layer_no(lay)
	# Axisymmetric part
	if axisymm:
		Ajjtm = spherical.mat_ebcm_axi_tm(C,m,'j','j',bnd.e12,bnd.e21)
		Ahjtm = spherical.mat_ebcm_axi_tm(C,m,'h','j',bnd.e12,bnd.e21)
		Ajjte = spherical.mat_ebcm_axi_te(C,m,'j','j',bnd.e12,bnd.e21)
		Ahjte = spherical.mat_ebcm_axi_te(C,m,'h','j',bnd.e12,bnd.e21)

		nshape = Ajjtm.shape
		nsize = nshape[0]

		K11 = bmat([[ identity(nsize) ],\
			    [ zeros(nshape)   ] ])

		K12 = bmat([[ zeros(nshape)   ],\
			    [ identity(nsize) ] ])

		K21tm = bmat([[ Ahjtm ],\
			      [-Ajjtm ] ])

		K21te = bmat([[ Ahjte ],\
			      [-Ajjte ] ])
		# For not-last boundaries we need Ah2,Bh2, etc.
		if not bnd.is_last:
			Ajhtm = spherical.mat_ebcm_axi_tm(C,m,'j','h',bnd.e12,bnd.e21)
			Ahhtm = spherical.mat_ebcm_axi_tm(C,m,'h','h',bnd.e12,bnd.e21)
			Ajhte = spherical.mat_ebcm_axi_te(C,m,'j','h',bnd.e12,bnd.e21)
			Ahhte = spherical.mat_ebcm_axi_te(C,m,'h','h',bnd.e12,bnd.e21)

			K22tm = bmat([[ Ahhtm ],\
				      [-Ajhtm ] ])

			K22te = bmat([[ Ahhte ],\
				      [-Ajhte ] ])
		else:
			K22tm=K22te=mat(zeros_like(K11))

	# Non-axisymmetric part
	else:
		Ajjtm = spherical.mat_ebcm_naxi_tm(C,m,'j','j',bnd.e12,bnd.e21)
		Ahjtm = spherical.mat_ebcm_naxi_tm(C,m,'h','j',bnd.e12,bnd.e21)
		Ajjte = spherical.mat_ebcm_naxi_te(C,m,'j','j',bnd.e12,bnd.e21)
		Ahjte = spherical.mat_ebcm_naxi_te(C,m,'h','j',bnd.e12,bnd.e21)

		nshape = Ajjtm.shape
		nsize = nshape[0]

		K11 = bmat([[ identity(nsize) ],\
			    [ zeros(nshape)   ] ])

		K12 = bmat([[ zeros(nshape)   ],\
			    [ identity(nsize) ] ])

		K21tm = bmat([[ Ahjtm ],\
			      [-Ajjtm ] ])

		K21te = bmat([[ Ahjte ],\
			      [-Ajjte ] ])

		# For not-last boundaries we need Ah2,Bh2, etc.
		if not bnd.is_last:
			Ajhtm = spherical.mat_ebcm_naxi_tm(C,m,'j','h',bnd.e12,bnd.e21)
			Ahhtm = spherical.mat_ebcm_naxi_tm(C,m,'h','h',bnd.e12,bnd.e21)
			Ajhte = spherical.mat_ebcm_naxi_te(C,m,'j','h',bnd.e12,bnd.e21)
			Ahhte = spherical.mat_ebcm_naxi_te(C,m,'h','h',bnd.e12,bnd.e21)

			K22tm = bmat([[ Ahhtm ],\
				      [-Ajhtm ] ])

			K22te = bmat([[ Ahhte ],\
				      [-Ajhte ] ])
		else:
			K22tm = K22te = mat(zeros_like(K11))

	return K11,K12,K21tm,K21te,K22tm,K22te

def pmm_bnd_system(C,m,bnd,axisymm=False):
	lay = bnd.layer_no
	C.set_layer_no(lay)
	if axisymm:
		Atm = spherical.mat_pmm_axi_tm(C,bnd.e12,bnd.shape.xv)
	else:
		Atm = spherical.mat_pmm_naxi_tm(C,m,bnd.e12,bnd.shape.xv)

	nsize = Atm.shape[0]/2
	nshape = (nsize,nsize)

	K11 = Atm[:,2*nsize:]
	K12 = Atm[:,:nsize]
	K21tm=K21te = Atm[:,nsize:2*nsize]
	#pdb.set_trace()

	if not bnd.is_last:
		K22tm=K22te=mat(zeros_like(K11))
	else:
		K22tm=K22te=mat(zeros_like(K11))

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

def meth_test_M0(meth_bnd_system,lab,nrange,accuracyLimit=None,ngauss="auto",conv_stop=True,iterative=True):
	print "\n","*"*60
	print "Convergence test: axisymmetric part"
	print "*"*60

	convlog_te = {'N':[],'delta':[]}
	convlog_tm = {'N':[],'delta':[]}

	# first we need at least one run of solver to have initial values
	Ms = [0]
	Cext_tm,c_sca_tm, Cext_te,c_sca_te = meth_n(meth_bnd_system,lab,nrange[0],ngauss,Ms,iterative)

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

svm = methods_factory(svm_bnd_system)
ebcm = methods_factory(ebcm_bnd_system)
pmm0 = methods_factory(pmm_bnd_system)

def pmm(lab,nrange=None,accuracyLimit=None,ngauss=200,conv_stop=True,conv_test=False,iterative=True):
	print "Warning: layered structures aren't supported"
	print "Warning: TE mode isn't supported"
	return pmm0(lab,nrange,accuracyLimit,ngauss,conv_stop,conv_test,iterative)

