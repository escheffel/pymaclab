import numpy as np
#from ordqz import dtgsen
#from qz import dgges
from pymaclab.linalg._drivers import qz
from pymaclab.linalg.qz import dgges, zgges
from scipy.io import loadmat

f = loadmat('./ordqzinputs.mat')
f2 = loadmat('./phipsi.mat')
#phiTr = f['phiTr'] 
#psiTr = f['psiTr']
#q = f['q']
#z = f['z']
select = f['selector']

phi = f2['phi']
psi = f2['psi']
npredet = 3
nvars = 8
njumpvars = nvars - npredet

phiTr,psiTr,sdim,alphar,alphai,beta, q, z, work, info = dgges(lambda x,y,z: (x+y)/z <=1.,phi,psi)
E = alphar/beta
E[np.isfinite(E)] = np.sort(E[np.isfinite(E)])
E[np.isinf(E)] = -np.inf
explos_index = np.where(np.abs(E)>1)[0]
if len(explos_index) < njumpvars:
    raise ValueError, "Not enough explosive roots"
elif len(explos_index) > njumpvars:
    raise ValueError, "Too many explosive roots, multiple RE equilibrium"
stable_index = np.abs(E)<1
#select = np.asfortranarray(np.squeeze(np.array(select, dtype=int)))
phiTr2,psiTr2,alphar,alphai,beta,q2,z2,m,pl,pr,dif,info = dtgsen(0,1,1,stable_index,phiTr,psiTr,q,z)
#selct = lambda x,y,z: (x + y)/z > 1
#select = f['selector']
