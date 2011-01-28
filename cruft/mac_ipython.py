#########################################################################
# INITIALISATION STUFF START
from __future__ import division
import sys
sys.path.append('/home/escheffel/phd/Programming/pytest/macroeconomics')
import os as OPS
import numpy as N
import pylab as P
import scipy as S
from numpy import matlib as MAT
from scikits import timeseries as TS
import macrolab
from macrolab import *
from macrolab import macrolab as mac
root = '/home/escheffel/phd/Programming/pytest/macroeconomics'
datapath = '/home/escheffel/phd/Programming/pytest/data'
modfpath = OPS.path.join(root,'macrolab/modfiles')
browserpath = '/usr/lib/firefox/firefox'
lapackpath = OPS.path.join(root,'macrolab/libs/')
lapackname = 'liblapack_macro.so'
mlabpath = OPS.path.join(root,'macrolab/mlab_files/')
mac.datapath = datapath
mac.modfpath = modfpath
mac.browserpath = browserpath
mac.lapackpath = lapackpath
mac.lapackname = lapackname
mac.txtedpath = 'winefish'
mac.texedpath = 'winefish'
mac.mlabpath = mlabpath
mac.locdic = locals()
mac.use_matlab = True
mac.use_anaderiv = False
mac.mk_hessian = False
# INITIALISATION STUFF END
###########################################################################
