#########################################################################
# INITIALISATION STUFF START
from __future__ import division
import sys
import os as OPS
if OPS.getenv('MACLAB'):
	HOMEDIR = OPS.getenv('MACLAB')
else:
	HOMEDIR = '/home/escheffel/workspace'
sys.path.append(HOMEDIR)
import numpy as N
import pylab as P
import scipy as S
from numpy import matlib as MAT
from scikits import timeseries as TS
from pymaclab import *
from pymaclab import macrolab as mac
root = OPS.path.join(HOMEDIR,'pymaclab')
datapath = OPS.path.join(root,'data')
modfpath = OPS.path.join(root,'modfiles')
browserpath = '/usr/lib/firefox/firefox'
lapackpath = OPS.path.join(root,'libs/')
lapackname = 'liblapack_macro.so'
mlabpath = OPS.path.join(root,'mlab_files/')
mac.datapath = datapath
mac.modfpath = modfpath
mac.browserpath = browserpath
mac.lapackpath = lapackpath
mac.lapackname = lapackname
mac.txtedpath = 'winefish'
mac.texedpath = 'winefish'
mac.pyedpath = 'winefish'
mac.pdfpath = 'acroread'
mac.mlabpath = mlabpath
mac.locdic = locals()
mac.use_matlab = True
mac.use_anaderiv = True
mac.mk_hessian = True
mac.ncpus = 2
# INITIALISATION STUFF END
###########################################################################
