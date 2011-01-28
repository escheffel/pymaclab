from __future__ import division
import wx
import wx.grid as gridlib
import trace
import time
import datetime
from .._hpfilter import hpfilt
# compiled hpfilter with 
# f2py -c hpfilter.f -m hpfilter
#import f90mods
# compiled isolab with
# f2py -c isolab.pyf -m isolab solab.f90 isolab.f90 -llapack
from ..isolab import isolab
import subprocess
import numpy.ma as MA
from numpy.ma import nomask
from scikits import timeseries as TSS
import copy as COP
import os as OPS
import re as RE
import sys as SYST
import re as RE
import string as STR
import scipy as S
from scipy import io
from scipy import optimize as O
import scipy.stats
import helpers as HLP
import numpy as N
import pylab as P
from numpy import matlib as MAT
from scipy import linalg as LIN
from scipy import stats as STA
from scikits.timeseries.lib import plotlib as TPL
import macrolab as MACLAB
import errors
from errors import *
import tempfile as TMF
import threading as THR
import glob
try:
    import sympy as SP
except:
    print "You need to install sympy"
try:
    import sympycore as SPC
except:
    print "You need to install sympycore"
    print "svn checkout http://sympycore.googlecode.com/svn/trunk/ sympycore"
import popen2
try:
    import pp as PP
except:
    print "You need to install pp"

#Define a list of the Greek Alphabet for Latex
greek_alph = ['alpha','beta','gamma','delta','epsilon',
	      'varepsilon','zeta','eta','theta','vartheta',
	      'iota','kappa','lambda','mu','nu','xi','pi',
	      'varpi','rho','varrho','sigma','varsigma',
	      'tau','upsilon','phi','varphi','chi','psi',
	      'omega','Gamma','Delta','Theta','Lambda',
	      'Xi','Pi','Sigma','Upsilon','Phi','Psi','Omega']

# Setting paths
root = ''
modfpath = ''
datapath = ''
browserpath = ''
txtedpath = ''
texedpath = ''
pyedpath = ''
pdfpath = ''
lapackpath = ''
lapackname = ''
mlabpath = ''
# Is matlab installed and can it therefore be used?
#TODO: this is not safe at all if matlab is not installed
try:
    import mlabraw
    use_matlab = True
except:
    use_matlab = False
    sess1 = None
# Should the Model Jacobian and Hessian be calculated symbolically?
use_anaderiv = True
# Should the Hessian be computed at all? (For 2nd order accurate methods)
mk_hessian = True
# Number of cores used for calculations
#ncpus = 2
from multiprocessing import cpu_count
ncpus = cpu_count()
ncpus = 1 # debug without parallel stuff
'''Model init level, 0 just pars, 1 pars and also steady state,
2 Jacobian (and perhaps Hessian)'''
initlev = 2
# Update level
updlev = 0

# Empty locdic for the locate helper function
locdic = {}


# Open mlab session
if use_matlab:
	sess1 = mlabraw.open('matlab - nojvm -nosplash')

# Setting standard VAR options
VAR_opts = {}
VAR_opts['IRF_periods'] = 20

###############THE VECTOR AUTOREGRESSION CLASS (WORKS)###############
class VAR:
	"""
	This is the vector autoregression class. It supports estimation,
	doing IRFs and other stuff.
	"""	
	def __init__(self,laglen=1,data='none',set_useconst='const'):
		self.VAR_attr = {}
		self.getdata(data)
		if set_useconst == 'noconst':
			self.setuseconst(0)
		elif set_useconst == 'const':
			self.setuseconst(1)
		self.setlaglen(laglen)

	def getdata(self,data):
		self.data = N.mat(data)
		self.VAR_attr['nuofobs'] = self.data.shape[0]
		self.VAR_attr['veclen'] = self.data.shape[1]

	def setuseconst(self,useconst):
		self.VAR_attr['useconst'] = useconst

	def setlaglen(self,laglen):	
		self.VAR_attr['laglen'] = laglen
		self.VAR_attr['avobs'] = self.VAR_attr['nuofobs'] - self.VAR_attr['laglen']

	#This is the OLS function and does just that, no input
	def ols(self):
		self.ols_results = {}
		data = self.data
		VAR_attr = self.VAR_attr
		veclen = VAR_attr['veclen']
		laglen = VAR_attr['laglen']
		nuofobs = VAR_attr['nuofobs']
		avobs = VAR_attr['avobs']
		useconst = VAR_attr['useconst']		
		y = data[laglen:,:]
		X = MAT.zeros((avobs,veclen*laglen))
		for x1 in range(0,laglen,1):
			X[:,x1*veclen:(x1+1)*veclen] = data[(laglen-1)-x1:(nuofobs-1)-x1,:]
		if self.VAR_attr['useconst'] == 1:
			X = N.hstack((MAT.ones((avobs,1)),X[:,:]))
		self.ols_results['y'] = y
		self.ols_results['X'] = X

		if useconst == 1:
			beta = MAT.zeros((veclen,1+veclen*laglen))
		else:
			beta = MAT.zeros((veclen,veclen*laglen))
		XX = X.T*X
		try:
			iXX = XX.I
		except:
			err = TS_err('Singular Matrix')
			return
		self.iXX = iXX
		for i in range(0,veclen,1):
			Xy = X.T*y[:,i]
			beta[i,:] = (iXX*Xy).T
		self.ols_results['beta'] = beta
		yfit = X*beta.T
		self.ols_results['yfit'] = yfit
		resid = y-yfit
		self.ols_results['resid'] = resid
		#Separate out beta's into B matrices
		self.ols_results['BBs'] = {}
		if self.VAR_attr['useconst'] == 0:
			i1 = 1
			for x1 in range(0,veclen*laglen,laglen):
				self.ols_results['BBs']['BB'+str(i1)] = beta[:,x1:x1+veclen]
				i1 = i1 + 1
		elif self.VAR_attr['useconst'] == 1:
			self.ols_results['BBs']['BB0'] = beta[:,0]
			i1 = 1
			for x in range(0,veclen*laglen,laglen):
				self.ols_results['BBs']['BB'+str(i1)] = beta[:,x+1:x+veclen+1]
				i1 = i1 + 1

		#Make variance-covariance matrix of residuals
		omega = (resid.T*resid)/avobs
		self.ols_results['omega'] = omega
		#Make variance-covariance matrix for est. coefficients under OLS
		omega_beta_ols = (N.multiply(iXX,N.diag(omega)[0]),)
		for i1 in range(1,veclen,1):
			omega_beta_ols = omega_beta_ols+(N.multiply(iXX,N.diag(omega)[i1]),)
		self.ols_results['omega_beta'] = omega_beta_ols
		#Extract just the diagonal variances of omega_beta_ols
		omega_beta_ols_va = (N.diag(omega_beta_ols[0]),)
		for i1 in range(1,veclen,1):
			omega_beta_ols_va = omega_beta_ols_va+(N.diag(omega_beta_ols[i1]),)
		self.ols_results['omega_beta_va'] = omega_beta_ols_va
		#Make variance-covariance matrix of est. coefficients under GLS
		XeeX = X.T*resid[:,0]*resid[:,0].T*X
		omega_beta_gls = (iXX*XeeX*iXX,)
		for i1 in range(1,veclen,1):
			XeeX = X.T*resid[:,i1]*resid[:,i1].T*X
			omega_beta_gls = omega_beta_gls+(iXX*XeeX*iXX,)
		self.ols_results['omega_beta_gls'] = omega_beta_gls
		#Extract just the diagonal variances of omega_beta_gls
		omega_beta_gls_va = N.diag(omega_beta_gls[0])
		for i1 in range(1,veclen,1):
			omega_beta_gls_va = (omega_beta_gls_va,)+(N.diag(omega_beta_gls[i1]),)
		self.ols_results['omega_beta_gls_va'] = omega_beta_gls_va
		sqresid = N.power(resid,2)
		rss = N.sum(sqresid)
		self.ols_results['rss'] = rss
		AIC = S.linalg.det(omega)+(2.0*laglen*veclen**2)/nuofobs
		self.ols_results['AIC'] = AIC
		BIC = S.linalg.det(omega)+(N.log(nuofobs)/nuofobs)*laglen*veclen**2
		self.ols_results['BIC'] = BIC

	def do_irf(self,spos=1,plen=VAR_opts['IRF_periods']):
		self.IRF_attr = {}
		self.IRF_attr['spos'] = spos
		self.IRF_attr['plen'] = plen
		VAR_attr = self.VAR_attr
		veclen = VAR_attr['veclen']
		laglen = VAR_attr['laglen']
		nuofobs = VAR_attr['nuofobs']
		avobs = VAR_attr['avobs']
		useconst = VAR_attr['useconst']
		self.irf_results = {}
		# Strip out the means of vector y
		data = self.data
		dmeans = N.mat(N.average(data,0))
		dmeans = N.mat(dmeans.tolist()*nuofobs)
		self.dmeans = dmeans
		dmdata = data - dmeans
		self.IRF_attr['dmdata'] = dmdata
		# Do OLS on de-meaned series and collect BBs
		if self.VAR_attr['useconst'] == 1:
			self.setuseconst(0)
			self.data = dmdata
			self.ols_comp()
			self.IRF_attr['beta'] = COP.deepcopy(self.ols_comp_results['beta'])
			self.ols()
			omega = COP.deepcopy(self.ols_results['omega'])
			self.IRF_attr['omega'] = COP.deepcopy(self.ols_results['omega'])
			self.setuseconst(1)
			self.data = data
			self.ols_comp()
		elif self.VAR_attr['useconst'] == 0:
			self.data = dmdata
			self.ols_comp()
			self.IRF_attr['beta'] = COP.deepcopy(self.ols_comp_results['beta'])
			self.ols()
			omega = self.ols_results['omega']
			self.IRF_attr['omega'] = COP.deepcopy(self.ols_results['omega'])
			self.data = data
			self.ols_comp()

		# Calculate Cholesky decomposition of omega
		A0 = N.mat(S.linalg.cholesky(omega))
		A0 = A0.T
		A0 = A0.I
		self.IRF_attr['A0'] = A0
		beta = self.IRF_attr['beta']

		#Calculate IRFs using ols_comp_results
		ee = MAT.zeros((veclen*laglen,1))
		ee[spos-1,:] = 1
		shock = N.vstack((A0.I*ee[0:veclen,:],ee[veclen:,:]))
		Gam = shock.T
		Gam_2 = Gam[0,0:veclen]
		for x1 in range(0,plen,1):
			Gam_X = beta*Gam[x1,:].T
			Gam_2 = N.vstack((Gam_2,Gam_X.T[:,0:veclen]))
			Gam = N.vstack((Gam,Gam_X.T))

		self.irf_results['Gammas'] = Gam_2

	#This does OLS in first-order companion form (stacked), no input
	def ols_comp(self):
		self.ols_comp_results = {}
		self.make_y_X_stack()
		VAR_attr = self.VAR_attr
		veclen = VAR_attr['veclen']
		laglen = VAR_attr['laglen']
		nuofobs = VAR_attr['nuofobs']
		avobs = VAR_attr['avobs']
		useconst = VAR_attr['useconst']
		veclen = veclen*laglen
		laglen = 1
		X = self.ols_comp_results['X']
		y = self.ols_comp_results['y']

		if useconst == 1:
			beta = MAT.zeros((veclen,1+veclen*laglen))
		else:
			beta = MAT.zeros((veclen,veclen*laglen))
		XX = X.T*X
		iXX = XX.I
		for i in range(0,veclen,1):
			Xy = X.T*y[:,i]
			beta[i,:] = (iXX*Xy).T

		veclen2 = VAR_attr['veclen']
		for x1 in range(0,beta.shape[0],1):
			for x2 in range(0,beta.shape[1],1):
				if beta[x1,x2] < 1e-7 and x1 > veclen2-1:
					beta[x1,x2] = 0.0

		for x1 in range(0,beta.shape[0],1):
			for x2 in range(0,beta.shape[1],1):
				if  0.9999999 < beta[x1,x2] < 1.0000001 and x1 > veclen2-1:
					beta[x1,x2] = 1.0

		self.ols_comp_results['beta'] = beta
		yfit = X*beta.T
		self.ols_comp_results['yfit'] = yfit
		resid = y-yfit
		self.ols_comp_results['resid'] = resid
		#Make variance-covariance matrix of residuals
		omega = (resid.T*resid)/avobs
		self.ols_comp_results['omega'] = omega
		#Make variance-covariance matrix for est. coefficients under OLS
		omega_beta_ols = (N.diag(omega)[0]*iXX,)
		for i1 in range(1,veclen,1):
			omega_beta_ols = omega_beta_ols + (N.diag(omega)[i1]*iXX,)
		self.ols_comp_results['omega_beta'] = omega_beta_ols
		#Extract just the diagonal variances of omega_beta_ols
		omega_beta_ols_va = (N.diag(omega_beta_ols[0]),)
		for i1 in range(1,veclen,1):
			omega_beta_ols_va = omega_beta_ols_va+(N.diag(omega_beta_ols[i1]),)
		self.ols_comp_results['omega_beta_va'] = omega_beta_ols_va
		#Make variance-covariance matrix of est. coefficients under GLS
		XeeX = X.T*resid[:,0]*resid[:,0].T*X
		omega_beta_gls = (iXX*XeeX*iXX,)
		for i1 in range(1,veclen,1):
			XeeX = X.T*resid[:,i1]*resid[:,i1].T*X
			omega_beta_gls = omega_beta_gls+(iXX*XeeX*iXX,)
		self.ols_comp_results['omega_beta_gls'] = omega_beta_gls
		#Extract just the diagonal variances of omega_beta_gls
		omega_beta_gls_va = (N.diag(omega_beta_gls[0]),)
		for i1 in range(1,veclen,1):
			omega_beta_gls_va = omega_beta_gls_va+(N.diag(omega_beta_gls[i1]),)
		self.ols_comp_results['omega_beta_gls_va'] = omega_beta_gls_va
		sqresid = N.power(resid,2)
		rss = N.sum(sqresid)
		self.ols_comp_results['rss'] = rss
		AIC = S.linalg.det(omega[0:veclen,0:veclen])+(2.0*laglen*veclen**2)/nuofobs
		self.ols_comp_results['AIC'] = AIC
		BIC = S.linalg.det(omega[0:veclen,0:veclen])+(N.log(nuofobs)/nuofobs)*laglen*veclen**2
		self.ols_comp_results['BIC'] = BIC

	#This is a function which sets lag by AIC, input is maximum lags over which to search
	def lag_by_AIC(self,lags):
		laglen = self.VAR_attr['laglen']
		lagrange = N.arange(1,lags+1,1)
		i1 = 0
		for i in lagrange:
			self.laglen = i
			self.ols()
			if i1 == 0:
				AIC = self.AIC
				i1 = i1 + 1
			else:
				AIC = append(AIC,self.AIC)
				i1 = i1 + 1
		index = argsort(AIC)
		i1 = 1
		for element in index:
			if element == 0:
				lag_AIC = i1
				break
			else:
				i1 = i1 + 1
		self.VAR_attr['lag_AIC'] = lag_AIC
		self.VAR_attr['laglen'] = lag_AIC

	#This is a function which sets lag by BIC, input is maximum lags over which to search
	def lag_by_BIC(self,lags):
		laglen = self.VAR_attr['laglen']
		lagrange = N.arange(1,lags+1,1)
		i1 = 0
		for i in lagrange:
			self.laglen = i
			self.ols()
			if i1 == 0:
				BIC = self.BIC
				i1 = i1 + 1
			else:
				BIC = append(BIC,self.BIC)
				i1 = i1 + 1
		index = argsort(BIC)
		i1 = 1
		for element in index:
			if element == 0:
				lag_BIC = i1
				break
			else:
				i1 = i1 + 1
		self.VAR_attr['lag_BIC'] = lag_BIC
		self.VAR_attr['laglen'] = lag_BIC

	#Auxiliary function, creates the y and X matrices for OLS_comp
	def make_y_X_stack(self):
		VAR_attr = self.VAR_attr
		veclen = VAR_attr['veclen']
		laglen = VAR_attr['laglen']
		nuofobs = VAR_attr['nuofobs']
		avobs = VAR_attr['avobs']
		useconst = VAR_attr['useconst']
		data = self.data

		y = data[laglen:,:]
		X = MAT.zeros((avobs,veclen*laglen))
		for x1 in range(0,laglen,1):
			X[:,x1*veclen:(x1+1)*veclen] = data[(laglen-1)-x1:(nuofobs-1)-x1,:]
		if self.VAR_attr['useconst'] == 1:
			X = N.hstack((MAT.ones((avobs,1)),X[:,:]))
		try:
			self.ols_results['y'] = y
			self.ols_results['X'] = X
		except:
			self.ols()
			self.ols_results['y'] = y
			self.ols_results['X'] = X

		if useconst == 0:
			y_stack = MAT.zeros((avobs,veclen*laglen))
			y_stack_1 = MAT.zeros((avobs,veclen*laglen))
			y_stack[:,0:veclen] = y
			y_stack_1 = X

			y_stack = N.hstack((y_stack[:,0:veclen],y_stack_1[:,0:veclen*(laglen-1)]))

			self.ols_comp_results['X'] = y_stack_1
			self.ols_comp_results['y'] = y_stack
		else:
			y_stack = MAT.zeros((avobs,veclen*laglen))
			y_stack_1 = MAT.zeros((avobs,1+veclen*laglen))
			y_stack_1[:,0] = MAT.ones((avobs,1))[:,0]
			y_stack[:,0:veclen] = y
			y_stack_1 = X

			y_stack = N.hstack((y_stack[:,0:veclen],y_stack_1[:,1:veclen*(laglen-1)+1]))

			self.ols_comp_results['X'] = y_stack_1
			self.ols_comp_results['y'] = y_stack
"""***********************************************************"""
################THE TIMESERIES DATABASE CLASS (WORKS)#################
class TSDataBase:
	"""
	This is the Time Series Database class.
	"""
	def __init__(self):
		self.datdic = {}
		self.modif = {}

	def nameimDatStr(self,filename='none'):
		'''
		This imports the DataStream datafile using the
		file's NAME as an argument.
		'''
		str_tmp1 = open(OPS.path.join(datapath,filename),'r')
		flist = str_tmp1.read().splitlines()[:]
		dat_Start = flist[0].split(',')[1].strip()[:]
		dat_End = flist[1].split(',')[1].strip()[:]
		dat_Freq = RE.sub('"','',flist[2].split(',')[1]).strip()[:]
		dat_Name = RE.sub('"','',flist[3].split(',')[1]).strip()[:]
		dat_Code = RE.sub('"','',flist[4].split(',')[1]).strip()[:]		
		dat_Curr = RE.sub('"','',flist[5].split(',')[1]).strip()[:]
		datadate = flist[6:]
		datacsv = ''
		i1=0
		for x in datadate:
			str_tmp1 = x.split(',')[1]
			datacsv = datacsv[:] + str_tmp1[:]+'\n'
			i1 = i1 + 1
		output = open(OPS.getcwd()+'/'+'tmp_001.csv', 'w')
		output.write(datacsv)
		pdata = P.load('tmp_001.csv',delimiter=',')
		OPS.remove(OPS.getcwd()+'/'+'tmp_001.csv')
		if dat_Freq == 'D':
			tsdata = data = TSS.time_series(pdata,start_date=TSS.Date(freq=dat_Freq,
										  year=int(dat_Start[6:10]),
										  month=int(dat_Start[3:5]),
										  day=int(dat_Start[0:2])))
		elif dat_Freq == 'M':
			tsdata = data = TSS.time_series(pdata,start_date=TSS.Date(freq=dat_Freq,
										  year=int(dat_Start[6:10]),
										  month=int(dat_Start[3:5])))
		elif dat_Freq == 'Q':
			tsdata = data = TSS.time_series(pdata,start_date=TSS.Date(freq=dat_Freq,
										  year=int(dat_Start[6:10]),
										  quarter=int(dat_Start[3:4])))

		self.datdic[dat_Code]={}
		self.datdic[dat_Code]['infile'] = tsdata
		self.datdic[dat_Code]['is_DatStr'] = True
		self.datdic[dat_Code]['Dat_Desc'] = dat_Name
		self.datdic[dat_Code]['DatStr_Code'] = dat_Code
		self.datdic[dat_Code]['imTimeStamp'] = HLP.now_is()
		self.datdic[dat_Code]['start_date'] = tsdata.start_date
		self.datdic[dat_Code]['end_date'] = tsdata.end_date
		self.datdic[dat_Code]['freq'] = tsdata.freqstr[0:1]
		self.datdic[dat_Code]['infile'].TSdesc = dat_Name
		self.datdic[dat_Code]['infile'].TSname = dat_Code
		self.datdic[dat_Code]['svalue'] = tsdata.dates[0].value
		self.datdic[dat_Code]['evalue'] = tsdata.dates[-1].value
		self.datdic[dat_Code]['alias'] = False

	def strimDatStr(self,csvfile):
		'''
		This imports the DataStream datafile using the
		entire file as a STRING as an argument.
		'''
		str_tmp1 = COP.deepcopy(csvfile)
		flist = str_tmp1.splitlines()[:]
		dat_Start = flist[0].split(',')[1].strip()[:]
		dat_End = flist[1].split(',')[1].strip()[:]
		dat_Freq = RE.sub('"','',flist[2].split(',')[1]).strip()[:]
		dat_Name = RE.sub('"','',flist[3].split(',')[1]).strip()[:]
		dat_Code = RE.sub('"','',flist[4].split(',')[1]).strip()[:]		
		dat_Curr = RE.sub('"','',flist[5].split(',')[1]).strip()[:]
		datadate = flist[6:]
		datacsv = ''
		i1=0
		for x in datadate:
			if i1 != len(datadate):
				str_tmp1 = x.split(',')[1]
				datacsv = datacsv + str_tmp1+',\n'
			else:
				str_tmp1 = x.split(',')[1]
				datacsv = datacsv + str_tmp1+'\n'
			i1 = i1 + 1
		output = open(OPS.getcwd()+'/'+'tmp_001.csv', 'w')
		output.write(datacsv)
		pdata = P.load('tmp_001.csv',delimiter=',')
		OPS.remove(OPS.getcwd()+'/'+'tmp_001.csv')
		if dat_Freq == 'D':
			tsdata = data = time_series(pdata,start_date=TSS.Date(freq=dat_Freq,
									      year=dat_Start[0:1],
									      month=dat_Start[3:4],
									      day=dat_Start[6:9]))
		elif dat_Freq == 'M':
			tsdata = data = time_series(pdata,start_date=TSS.Date(freq=dat_Freq,
									      year=dat_Start[0:1],
									      month=dat_Start[3:4]))
		elif dat_Freq == 'Q':
			tsdata = data = time_series(pdata,start_date=TSS.Date(freq=dat_Freq,
									      year=dat_Start[0:1],
									      quarter=dat_Start[3:4]))

		self.datdic[dat_Name]={}
		self.datdic[dat_Name]['infile'] = tsdata
		self.datdic[dat_Name]['is_DatStr'] = True
		self.datdic[dat_Name]['DatStr_Code'] = dat_Code
		self.datdic[dat_Name]['imTimeStamp'] = HLP.now_is()

	def tsim(self,datname,tsdesc,tsin):
		self.datdic[datname]={}
		self.datdic[datname]['infile'] = tsin
		self.datdic[datname]['is_DatStr'] = False
		self.datdic[datname]['Dat_Desc'] = tsdesc
		self.datdic[datname]['DatStr_Code'] = None
		self.datdic[datname]['imTimeStamp'] = HLP.now_is()
		self.datdic[datname]['start_date'] = tsin.start_date
		self.datdic[datname]['end_date'] = tsin.end_date
		self.datdic[datname]['freq'] = tsin.freqstr[0:1]
		self.datdic[datname]['infile'].TSdesc = tsdesc
		self.datdic[datname]['infile'].TSname = datname
		self.datdic[datname]['svalue'] = tsin.dates[0].value
		self.datdic[datname]['evalue'] = tsin.dates[-1].value
		self.datdic[datname]['alias'] = False

	def isDatStreamCSV(self,csvfile):
		str_tmp1 = COP.deepcopy(csvfile)
		flist = str_tmp1.splitlines()[:]
		if '"Start"' in flist[0] and \
		   '"End"' in flist[1] and \
		   '"Frequency"' in flines[2] and \
		   '"Names"' in flist[3] and \
		   '"Code"' in flist[4] and \
		   '"CURRENCY"' in flist[5]:
			return True
		else:
			return False

	def mkhpf(self,tsname,tsout,lam=1600):
		tsinf = self.datdic[tsname]['infile']
		tsin = self.datdic[tsname]
		tsoutf = hpfilt(tsinf,N.zeros((N.shape(tsinf)[0],3)),
			N.shape(tsinf)[0],lam,0)
		tsoutf = TSS.time_series\
		       (tsoutf,start_date=tsinf.start_date,freq=tsin['freq'])
		self.tsim(tsout,tsin['Dat_Desc']+',hp-filtered',tsoutf)

	def mklog(self,tsname,tsout):
		tsinf = self.datdic[tsname]['infile']
		tsin = self.datdic[tsname]
		tsoutf = N.log(tsinf)
		self.tsim(tsout,tsin['Dat_Desc']+',logarithmic',tsoutf)

	def mkexp(self,tsname,tsout):
		tsinf = self.datdic[tsname]['infile']
		tsin = self.datdic[tsname]
		tsoutf = N.exp(tsinf)
		self.tsim(tsout,tsin['Dat_Desc']+',to base e',tsoutf)

	def mkalias(self,dbtsname,alname):
		if self.datdic.has_key(dbtsname):
			self.datdic[dbtsname]['alias'] = alname

	def mkmodif(self,freq):
		list_tmp = []
		for x1 in self.datdic.items():
			if x1[1]['alias'] != 'None' and x1[1]['freq'] == freq:
				list_tmp = list_tmp + [[x1[0],x1[1]['alias'],x1[1]['svalue'],x1[1]['evalue']],]
		sval_list = []
		for x1 in list_tmp:
			sval_list = sval_list + [x1[2],]
		eval_list = []
		for x1 in list_tmp:
			eval_list = eval_list + [x1[3],]
		max_sval = max(sval_list)
		min_eval = min(eval_list)
		s_date = TSS.Date(freq,value=max_sval)
		e_date = TSS.Date(freq,value=min_eval)
		list_tmp = []
		for x1 in self.datdic.items():
			if x1[1]['alias'] != 'None' and x1[1]['freq'] == freq:
				self.modif[x1[1]['alias']] = x1[1]['infile'][s_date:e_date]

	def mkvar(self,varlag=1,varord='None',spos='None',useconst='const'):
		tmp_dic={}
		for x1 in varord:
			tmp_dic[x1[0]] = self.modif[x1[0]]
		dat_mat = MAT.zeros((N.shape(N.mat(tmp_dic.items()[0][1].data).T)[0],len(varord)))
		self.datmat = dat_mat
		for x1 in varord:
			self.datmat[:,x1[1]-1] = N.mat(tmp_dic[x1[0]].data).T
		self.VAR = VAR(varlag,self.datmat,useconst)
		self.VAR.ols()
		self.VAR.ols_comp()
		if spos != 'None':
			self.VAR.do_irf(spos)
			self.VAR.IRF_attr['shock_var'] = varord[[x1[1] for x1 in varord].index(spos-1)][0]
"""***********************************************************"""
##################THE DSGE MODEL CLASS (WORKS)#####################
class DSGEmodel:
	'''
	This is the macrolab DSGEmodel class. It is the main class
	of the packages and creates DSGE model instances which have
	lots of interesting solution and other features.
	'''
	# Initializes the absolute basics, errors difficult to occur
	def __init__(self,ffile=None,dbase=None):
		# Set no author
		self.setauthor()
		# Open the switches dic and initiate
		self.switches = {}
		self.switches['errocc'] = '0'
		self.switches['ss_suc'] = ['0','0']
		# Attach some model attributes
		self.modfile = ffile
		self.dbase = dbase
	# Initializes all of the rest, errors can occur here ! (steady state, jacobian, hessian)
	def init2(self):
		# Create None tester regular expression
		_nreg = '^\s*None\s*$'
		nreg = RE.compile(_nreg)

		self.txtpars = TXTparser(self.modfile)
		# Start to populate the model from file
		self.pop(input=self.txtpars)
		# Attach the data from database
		if self.dbase != None:
			self.getdata(dbase=self.dbase)
		if initlev == 0: return #NOTE: this is a global variable. get rid of.
################## STEADY STATE CALCULATIONS !!! ####################
		self.sssolvers = SSsolvers()
		# Solve for steady-state using fsolve
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['ssm'][0]]) == 0:
			intup = (self.ssys_list,self.ssidic,self.paramdic)
			self.sssolvers.fsolve = Fsolve(intup)
			self.sssolvers.fsolve.solve()
			if self.sssolvers.fsolve.ier == 1:
				self.sstate = self.sssolvers.fsolve.fsout
				self.numssdic = self.sssolvers.fsolve.fsout
				self.switches['ss_suc'] = ['1','1']
			else:
				self.switches['ss_suc'] = ['1','0']
		# Solve for steady-state using manss
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['sss'][0]]) == 0:
			if self.switches['ss_suc'] == ['1','1']:
				alldic = {}
				alldic.update(self.sstate)
				alldic.update(self.paramdic)
				intup = (self.manss_sys,alldic)
				self.sssolvers.manss = Manss(intup)
				self.sssolvers.manss.solve()
				self.sstate.update(self.sssolvers.manss.sstate)
			else:
				intup = (self.manss_sys,self.paramdic)
				self.sssolvers.manss = Manss(intup)
				self.sssolvers.manss.solve()
				self.sstate = self.sssolvers.manss.sstate
		if initlev == 1: return

		# No populate more with stuff that needs steady state
		self.pop2(self.txtpars)

		# Open the model solution tree branch
		self.modsolvers = MODsolvers()
		######################## LINEAR METHODS !!! ############################
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['modeq'][0]]) == 0:
			# Open the matlab Uhlig object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1,
				 self.vardic)
			self.modsolvers.matuhlig = MatUhlig(intup)
			# Open the native Uhlig object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1)
			self.modsolvers.pyuhlig = PyUhlig(intup)
			# Open the matlab Klein object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1)
			self.modsolvers.matklein = MatKlein(intup)
			# Open the Fortran Klein object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1)
			self.modsolvers.forklein = ForKlein(intup)
		################## 1ST-ORDER NON-LINEAR METHODS !!! ##################
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['focs'][0]]) == 0:

			# First, create the Jacobian and (possibly-->mk_hessian==True?) Hessian
			if use_anaderiv:
				if ncpus > 1 and mk_hessian:
					self.mkjahepp()
				elif ncpus > 1 and not mk_hessian:
					self.mkjahepp()
				else:
					self.mkjahe()
			else:
				self.mkjahenmat()

			# Open the MatWood object
			intup = (self.jAA,self.jBB,
				 self.nexo,self.ncon,
				 self.nendo,sess1)
			self.modsolvers.matwood = MatWood(intup)
			# Open the Fortran KleinD object
			if 'nlsubsys' in dir(self):
				intup = (self.numj,
					 self.nendo,self.nexo,
					 self.ncon,self.sigma,
					 self.jAA,self.jBB,
					 self.vardic,self.vdic,
					 self.modname,self.audic,
					 self.numjl,
					 self.nother)
			else:
				intup = (self.numj,
					 self.nendo,self.nexo,
					 self.ncon,self.sigma,
					 self.jAA,self.jBB,
					 self.vardic,self.vdic,
					 self.modname,self.audic)
			self.modsolvers.forkleind = ForKleinD(intup)
		################## 2ND-ORDER NON-LINEAR METHODS !!! ##################
			if sum([nreg.search(x)!=None for x in self.txtpars.secs['vcvm'][0]]) == 0 and\
			   'numh' in dir(self):
				# Open the MatKlein2D object
				if 'nlsubsys' in dir(self):
					intup = (self.numj,self.numh,
						 self.nendo,self.nexo,
						 self.ncon,self.sigma,
						 self.jAA,self.jBB,
						 self.vardic,self.vdic,
						 self.modname,self.audic,
						 self.numjl,self.numhl,
						 self.nother,sess1)
				else:
					intup = (self.numj,self.numh,
						 self.nendo,self.nexo,
						 self.ncon,self.sigma,
						 self.jAA,self.jBB,
						 self.vardic,self.vdic,
						 self.modname,self.audic,
						 sess1)
				self.modsolvers.matklein2d = MatKlein2D(intup)
				# Open the PyKlein2D object
				intup = intup[:-1]
				self.modsolvers.pyklein2d = PyKlein2D(intup)

		if 'jAA' in dir(self):
			self.mkeigv()
			
		# Wrap the paramdic at the end of model initialization
		self.params = dicwrap(self,initlev,nreg)

	# Initial population method of model, does NOT need steadt states
	def pop(self,input=None):
		# Create None tester regular expression
		_nreg = '^\s*None\s*$'
		nreg = RE.compile(_nreg)

		# Get all information from the varvec section
		def mkvarinfo():
			vtexp = RE.compile('^.*?\[.*?\]\s*(?P<vari>.+?)\s*:\s*(?P<varn>.+?)\s*\{(?P<vtype>[^\[]+)}\s*(?P<mod>\[.*?\]){0,1}')
			viiexp = RE.compile('^.*?\[.*?\]\s*(?P<vari>.+?)\s*:\s*(?P<viid>.+?)\s*:\s*(?P<varn>.+?)\s*\{(?P<vtype>[^\[]+)}\s*(?P<mod>\[.*?\]){0,1}')
			voexp = RE.compile('^.*?\[.*?\]\s*(?P<vari>@.+?):(?P<varn>[^\[]+)\s*(?P<mod>\[.*?\]){0,1}')
			vardic = {}
			audic = {}
			vardic['endo'] = {}
			audic['endo'] = {}
			vardic['exo'] = {}
			audic['exo'] = {}
			vardic['con'] = {}
			audic['con'] = {}
			vardic['other'] = {}
			vardic['endo']['var'] = []
			vardic['endo']['mod'] = []
			vardic['exo']['var'] = []
			vardic['exo']['mod'] = []
			vardic['con']['var'] = []
			vardic['con']['mod'] = []
			vardic['other']['var'] = []
			vardic['other']['mod'] = []
			audic['endo']['var'] = []
			audic['endo']['mod'] = []
			audic['con']['var'] = []
			audic['con']['mod'] = []
			audic['exo']['var'] = []
			audic['exo']['mod'] = []

			for x in input.secs['varvec'][0]:
				if viiexp.search(x):
					ma = viiexp.search(x)
					vari = ma.group('vari').strip()
					viid = ma.group('viid').strip()
					varn = ma.group('varn').strip()
					vtype = ma.group('vtype').strip()
					mods = ma.group('mod')
					vardic[vtype]['var'].append([vari,varn,viid])
					if mods != None:
						mods = mods.strip()
						if ',' in mods:
							vardic[vtype]['mod'].append(mods[1:-1].strip().split(','))
						else:
							vardic[vtype]['mod'].append([mods[1:-1].strip(),])
					else:
						vardic[vtype]['mod'].append([])

				elif vtexp.search(x):
					ma = vtexp.search(x)
					vari = ma.group('vari').strip()
					varn = ma.group('varn').strip()
					vtype = ma.group('vtype').strip()
					mods = ma.group('mod')
					vardic[vtype]['var'].append([vari,varn])
					if mods != None:
						mods = mods.strip()
						if ',' in mods:
							vardic[vtype]['mod'].append(mods[1:-1].strip().split(','))
						else:
							vardic[vtype]['mod'].append([mods[1:-1].strip(),])
					else:
						vardic[vtype]['mod'].append([])

				elif voexp.search(x):
					ma = voexp.search(x)
					# Slice off the @
					vari = ma.group('vari')[1:].strip()
					varn = ma.group('varn').strip()
					mods = ma.group('mod')
					vardic['other']['var'].append([vari,varn])
					if mods != None:
						mods = mods.strip()
						if ',' in mods:
							vardic['other']['mod'].append(mods[1:-1].strip().split(','))
						else:
							vardic['other']['mod'].append([mods[1:-1].strip(),])
					else:
						vardic['other']['mod'].append([])

			self.nendo = len(vardic['endo']['var'])
			self.nexo = len(vardic['exo']['var'])
			self.ncon = len(vardic['con']['var'])
			self.nother = len(vardic['other']['var'])
			self.nstat = self.nendo+self.nexo
			self.nall = self.nstat+self.ncon
			self.vardic = vardic
			self.audic = audic

		if sum([nreg.search(x)!=None for x in input.secs['varvec'][0]]) == 0:	
			mkvarinfo()

		#Extract the model description into string list
		def mkdesc():
			self.mod_desc = input.secs['mod'][0][0]
		if sum([nreg.search(x)!=None for x in input.secs['mod'][0]]) == 0:
			mkdesc()

		# Extract the model name
		def mkname():
			for x in input.secs['info'][0]:
				if 'Name' in x:
					self.modname = x.split('=')[1].replace(';','').strip()
		if sum([nreg.search(x)!=None for x in input.secs['info'][0]]) == 0:
			mkname()

		# Extract parameters into dictionary
		def mkparam():
			param = {}
			for x in input.secs['para'][0]:
				list_tmp = x.split(';')
				list_tmp = list_tmp[0].split('=')[:]
				str_tmp1 = list_tmp[0].strip()
				str_tmp2 = list_tmp[1].strip()
				param[str_tmp1] = eval(str_tmp2)
				locals()[str_tmp1] = eval(str_tmp2)
			return param
		if sum([nreg.search(x)!=None for x in input.secs['para'][0]]) == 0:
			self.paramdic = mkparam()
		# Collect information on manual (closed-form) steady state
		def colecmanss():
			# Join multiline steady state definitions
			mansys = self.txtpars.secs['sss'][0]
			list_tmp1 = []
			i1=0
			counter=0
			for x in mansys:
				if '...' in x:
					counter = counter + 1
				elif '...' not in x:
					if counter == 0:
						list_tmp1.append(x)
					elif counter > 0:
						str_tmp = ''
						for y in mansys[i1-counter:i1+1]:
							str_tmp = str_tmp + y.replace('...','')
						list_tmp1.append(str_tmp)
						counter = 0	
				i1=i1+1
			return list_tmp1

		if sum([nreg.search(x)!=None for x in input.secs['sss'][0]]) == 0:
			self.manss_sys = colecmanss()
		# Collect info on numerical steady state
		def colecnumss():
			_mreg = '[a-zA-Z]*_bar\s*=\s*[0-9]*\.[0-9]*'
			_mreg2 = '[a-zA-Z]*\s*=\s*[0-9]*\.[0-9]*'
			mreg = RE.compile(_mreg+'|'+_mreg2)
			indx = []
			ssidic={}
			list_tmp = []
			i1=0
			counter=0
			for x in input.secs['ssm'][0]:
				if mreg.search(x):
					ma = mreg.search(x)
					str1 = ma.group().split('=')[0].strip()
					str2 = ma.group().split('=')[1].strip()
					ssidic[str1] = eval(str2)
					indx.append(i1)
				elif not mreg.search(x) and '...' in x:
					counter = counter + 1
				elif not mreg.search(x) and '...' not in x:
					if counter == 0:
						list_tmp.append(x.replace(';','').split(']')[1].split('=')[0].strip())
					elif counter > 0:
						str_tmp = ''
						for y in input.secs['ssm'][0][i1-counter:i1+1]:
							str_tmp = str_tmp + y.replace('...','').strip()
						list_tmp.append(str_tmp.split(']')[1].split('=')[0].replace(';','').strip())
						counter = 0	
				i1=i1+1
			return (ssidic,list_tmp)
		if sum([nreg.search(x)!=None for x in input.secs['ssm'][0]]) == 0:
			self.ssidic,self.ssys_list = colecnumss()
	# 2dn stage population of model, does NEED the steady state
	def pop2(self,input=None):
		# Create None tester regular expression
		_nreg = '^\s*None\s*$'
		nreg = RE.compile(_nreg)

		# Creates nonlinear foc system		
		def mknlinsys():

			# Make the non-linear system by joining lines and stripping
			list_tmp1 = []
			i1=0
			counter=0
			for x in input.secs['focs'][0]:
				if '...' in x:
					counter = counter + 1
				elif '...' not in x:
					if counter == 0:
						list_tmp1.append(x.replace(';','').split(']')[1].strip())
					elif counter > 0:
						str_tmp = ''
						for y in input.secs['focs'][0][i1-counter:i1+1]:
							str_tmp = str_tmp + y.replace('...','').replace(';','').strip()
						list_tmp1.append(str_tmp.split(']')[1].strip())
						counter = 0	
				i1=i1+1
			for i,x in enumerate(list_tmp1):
				list_tmp1[i] = x.split('=')[0].strip()

			# Check and do substitutions
			if sum([nreg.search(x)!=None for x in input.secs['vsfocs'][0]]) == 0:

				# Make the substitution system by joining lines and stripping
				list_tmp2 = []
				i1=0
				counter=0
				for x in input.secs['vsfocs'][0]:
					if '...' in x:
						counter = counter + 1
					elif '...' not in x:
						if counter == 0:
							list_tmp2.append(x.replace(';','').split(']')[1].strip())
						elif counter > 0:
							str_tmp = ''
							for y in input.secs['vsfocs'][0][i1-counter:i1+1]:
								str_tmp = str_tmp + y.replace('...','').replace(';','').strip()
							list_tmp2.append(str_tmp.split(']')[1].strip())
							counter = 0	
					i1=i1+1

				for i,x in enumerate(list_tmp2):
					list_tmp2[i] = [x.split('=')[0].strip(),x.split('=')[1].strip()]

				# Replace substitutions inside substitutions
				mreg = RE.compile('@(E\(t.*?\)\|){0,1}.*?\(t.*?\)')
				for i,x in enumerate(list_tmp2):
					while mreg.search(list_tmp2[i][1]):
						ma = mreg.search(list_tmp2[i][1])
						pos = ma.span()[0]
						poe = ma.span()[1]
						indx = [x[0] for x in list_tmp2].index(ma.group())
						list_tmp2[i][1] = list_tmp2[i][1][:pos]+'('+list_tmp2[indx][1]+')'+list_tmp2[i][1][poe:]


				# substitute out in main nonlinear equation system
				mreg = RE.compile('@(E\(t.*?\)\|){0,1}.*?\(t.*?\)')
				for i1,line in enumerate(list_tmp1):
					while mreg.search(list_tmp1[i1]):
						ma = mreg.search(list_tmp1[i1])
						subv = ma.group()
						pos = ma.span()[0]
						poe = ma.span()[1]
						indx = [x[0] for x in list_tmp2].index(subv)
						list_tmp1[i1] = list_tmp1[i1][:pos]+'('+list_tmp2[indx][1]+')'+list_tmp1[i1][poe:]


			def mkaug1(insys,othersys):			
				# Determine the lengths of augmented vars

				list_tmp2 = COP.deepcopy(insys)
				list_tmp1 = COP.deepcopy(othersys)

				endosp = []
				for x in self.vardic['endo']['var']:
					endosp = endosp + [[x,[-1,0]]]
				exosp = []
				for x in self.vardic['exo']['var']:
					exosp = exosp + [[x,[0,1]]]
				consp = []
				for x in self.vardic['con']['var']:
					consp = consp + [[x,[0,1]]]

				alldic = {}
				alldic.update(self.paramdic)
				alldic.update(self.sstate)

				spvdic = {}
				spvdic['endo'] = endosp
				spvdic['exo'] = exosp
				spvdic['con'] = consp

				spvdic2 = COP.deepcopy(spvdic)

				endoli = [x[0].split('(')[0].strip() for x in self.vardic['endo']['var']]
				exoli = [x[0].split('(')[0].strip() for x in self.vardic['exo']['var']]
				conli = [x[0].split('(')[0].strip() for x in self.vardic['con']['var']]


				patup = ('{-100,100}|None','all','{-100,100}')
				count = 0
				for i1,line in enumerate(list_tmp2):
					iterob = self.vreg(patup,line,True,'max')
					if iterob:
						iterob = list(iterob)
					else:
						continue
					iterob.reverse()
					for vma in iterob:
						vtype = vma[1][0]
						vartime = vma[2][2]
						vari = vma[2][1]
						pos = vma[3][0]
						poe = vma[3][1]
						if vtype == 'endo':
							indx = endoli.index(vari)
						elif vtype == 'con':
							indx = conli.index(vari)
						elif vtype == 'exo':
							indx = exoli.index(vari)
						else:
							continue

						# Check for endo
						if vtype == 'endo' and int(vartime) > spvdic['endo'][indx][1][1]:
							spvdic2['endo'][indx][1][1] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_F'+tind+'(t)'
							newname =  vari+'_F'+str(abs(int(vartime)))
							list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp2[i1][poe:]
							for i2 in range(int(vartime)):
								tind = (5-len(str(abs(i2+1)-1)))*'0'+str(abs(i2+1)-1)
								newvar = vari+'_F'+tind+'(t)'
								newname =  vari+'_F'+str(abs(i2))
								if [newvar,newname] not in self.vardic['con']['var']:
									self.vardic['con']['var'].append([newvar,newname])
									self.vardic['con']['mod'].append(self.vardic['endo']['mod'][indx])
									self.audic['con']['var'].append([newvar,newname])
									self.audic['con']['mod'].append(self.vardic['endo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						elif vtype == 'endo' and int(vartime) < spvdic['endo'][indx][1][0]:
							spvdic2['endo'][indx][1][0] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_B'+tind+'(t)'
							newname =  vari+'_B'+str(abs(int(vartime)))
							list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp2[i1][poe:]
							for i2 in range(1,abs(int(vartime))):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_B'+tind+'(t)'
								newname =  vari+'_B'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['endo']['var']:
									self.vardic['endo']['var'].append([newvar,newname])
									self.vardic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
									self.audic['endo']['var'].append([newvar,newname])
									self.audic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						# Check for exo
						if vtype == 'exo' and int(vartime) > spvdic['exo'][indx][1][1]:
							spvdic2['exo'][indx][1][1] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_F'+tind+'(t)'
							newname =  vari+'_F'+str(abs(int(vartime)-1))
							list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp2[i1][poe:]
							for i2 in range(int(vartime)):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_F'+tind+'(t)'
								newname =  vari+'_F'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['con']['var']:
									self.vardic['con']['var'].append([newvar,newname])
									self.vardic['con']['mod'].append(self.vardic['exo']['mod'][indx])
									self.audic['con']['var'].append([newvar,newname])
									self.audic['con']['mod'].append(self.vardic['exo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						elif vtype == 'exo' and int(vartime) < spvdic['exo'][indx][1][0]:
							spvdic2['exo'][indx][1][0] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_B'+tind+'(t)'
							newname =  vari+'_B'+str(abs(int(vartime)))
							list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp2[i1][poe:]
							for i2 in range(abs(int(vartime))):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_B'+tind+'(t)'
								newname =  vari+'_B'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['endo']['var']:
									self.vardic['endo']['var'].append([newvar,newname])
									self.vardic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
									self.audic['endo']['var'].append([newvar,newname])
									self.audic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						# Check for con
						if vtype == 'con' and int(vartime) > spvdic['con'][indx][1][1]:
							spvdic2['con'][indx][1][1] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_F'+tind+'(t)'
							newname =  vari+'_F'+str(abs(int(vartime)-1))
							list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp2[i1][poe:]
							for i2 in range(int(vartime)):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_F'+tind+'(t)'
								newname =  vari+'_F'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['con']['var']:
									self.vardic['con']['var'].append([newvar,newname])
									self.vardic['con']['mod'].append(self.vardic['con']['mod'][indx])
									self.audic['con']['var'].append([newvar,newname])
									self.audic['con']['mod'].append(self.vardic['con']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						elif vtype == 'con' and int(vartime) < spvdic['con'][indx][1][0]:
							spvdic2['con'][indx][1][0] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_B'+tind+'(t)'
							newname =  vari+'_B'+str(abs(int(vartime)))
							list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp2[i1][poe:]
							for i2 in range(abs(int(vartime))):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_B'+tind+'(t)'
								newname =  vari+'_B'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['endo']['var']:
									self.vardic['endo']['var'].append([newvar,newname])
									self.vardic['endo']['mod'].append(self.vardic['con']['mod'][indx])
									self.audic['endo']['var'].append([newvar,newname])
									self.audic['endo']['mod'].append(self.vardic['con']['mod'][indx])
									self.sstate[newvar.split('(')[0]+'_bar'] = alldic[vari+'_bar']
								continue

				# Now change the system to include possible augmented variables
				endo_r = filter(lambda x: x[1] != [-1,0], spvdic2['endo']) 
				if endo_r:
					endo_r = [[x[0],[abs(x[1][0]+1),x[1][1]]] for x in endo_r ]
					# Create lags and forwards equations
					for vari in endo_r:
						for lag in range(abs(vari[1][0])):
							tind = (5-len(str(lag+2)))*'0'+str(lag+2)
							if lag == 0:
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t-1)')
							else:
								tind1 = (5-len(str(lag+1)))*'0'+str(lag+1)
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
						for lead in range(vari[1][1]):
							tind = (5-len(str(lead)))*'0'+str(lead)
							if lead == 0:
								if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0] not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0])
							else:
								tind1 = (5-len(str(lead-1)))*'0'+str(lead-1)
								if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
				exo_r = filter(lambda x: x[1] != [0,1], spvdic2['exo'])
				if exo_r:
					exo_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in exo_r ]
					# Create lags and forwards equations
					for vari in exo_r:
						for lag in range(abs(vari[1][0])):
							tind = (5-len(str(lag+1)))*'0'+str(lag+1)
							if lag == 0:
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)')
							else:
								tind1 = (5-len(str(lag)))*'0'+str(lag)
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
						if vari[1][1] > 0:
							for lead in range(vari[1][1]+1):
								tind = (5-len(str(lead+1)))*'0'+str(lead+1)
								if lead == 0:
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
								else:
									tind1 = (5-len(str(lead)))*'0'+str(lead)
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
				con_r = filter(lambda x: x[1] != [0,1], spvdic2['con']) 
				if con_r:
					con_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in con_r ]
					# Create lags and forwards equations
					for vari in con_r:
						for lag in range(abs(vari[1][0])):
							tind = (5-len(str(lag+1)))*'0'+str(lag+1)
							if lag == 0:
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)')
							else:
								tind1 = (5-len(str(lag)))*'0'+str(lag)
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1: 
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
						if vari[1][1] > 0:
							for lead in range(int(vari[1][1])+1):
								tind = (5-len(str(lead+1)))*'0'+str(lead+1)
								if lead == 0:
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
								else:
									tind1 = (5-len(str(lead)))*'0'+str(lead)
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
						
						

				return (list_tmp1,list_tmp2)

			def mkaug2(insys):			
				# Determine the lengths of augmented vars

				list_tmp1 = COP.deepcopy(insys)

				endosp = []
				for x in self.vardic['endo']['var']:
					endosp = endosp + [[x,[-1,0]]]
				exosp = []
				for x in self.vardic['exo']['var']:
					exosp = exosp + [[x,[0,1]]]
				consp = []
				for x in self.vardic['con']['var']:
					consp = consp + [[x,[0,1]]]

				alldic = {}
				alldic.update(self.paramdic)
				alldic.update(self.sstate)

				spvdic = {}
				spvdic['endo'] = endosp
				spvdic['exo'] = exosp
				spvdic['con'] = consp

				spvdic2 = COP.deepcopy(spvdic)

				endoli = [x[0].split('(')[0].strip() for x in self.vardic['endo']['var']]
				exoli = [x[0].split('(')[0].strip() for x in self.vardic['exo']['var']]
				conli = [x[0].split('(')[0].strip() for x in self.vardic['con']['var']]


				patup = ('{-100,100}|None','all','{-100,100}')
				count = 0
				for i1,line in enumerate(list_tmp1):
					iterob = self.vreg(patup,line,True,'max')
					if iterob:
						iterob = list(iterob)
					else:
						continue
					iterob.reverse()
					for vma in iterob:
						vtype = vma[1][0]
						vartime = vma[2][2]
						vari = vma[2][1]
						pos = vma[3][0]
						poe = vma[3][1]
						if vtype == 'endo':
							indx = endoli.index(vari)
						elif vtype == 'con':
							indx = conli.index(vari)
						elif vtype == 'exo':
							indx = exoli.index(vari)
						else:
							continue

						# Check for endo
						if vtype == 'endo' and int(vartime) > spvdic['endo'][indx][1][1]:
							spvdic2['endo'][indx][1][1] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_F'+tind+'(t)'
							newname =  vari+'_F'+str(abs(int(vartime)))
							list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp1[i1][poe:]
							for i2 in range(int(vartime)):
								tind = (5-len(str(abs(i2+1)-1)))*'0'+str(abs(i2+1)-1)
								newvar = vari+'_F'+tind+'(t)'
								newname =  vari+'_F'+str(abs(i2))
								if [newvar,newname] not in self.vardic['con']['var']:
									self.vardic['con']['var'].append([newvar,newname])
									self.vardic['con']['mod'].append(self.vardic['endo']['mod'][indx])
									self.audic['con']['var'].append([newvar,newname])
									self.audic['con']['mod'].append(self.vardic['endo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						elif vtype == 'endo' and int(vartime) < spvdic['endo'][indx][1][0]:
							spvdic2['endo'][indx][1][0] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_B'+tind+'(t)'
							newname =  vari+'_B'+str(abs(int(vartime)))
							list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp1[i1][poe:]
							for i2 in range(1,abs(int(vartime))):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_B'+tind+'(t)'
								newname =  vari+'_B'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['endo']['var']:
									self.vardic['endo']['var'].append([newvar,newname])
									self.vardic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
									self.audic['endo']['var'].append([newvar,newname])
									self.audic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						# Check for exo
						if vtype == 'exo' and int(vartime) > spvdic['exo'][indx][1][1]:
							spvdic2['exo'][indx][1][1] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_F'+tind+'(t)'
							newname =  vari+'_F'+str(abs(int(vartime)-1))
							list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp1[i1][poe:]
							for i2 in range(int(vartime)):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_F'+tind+'(t)'
								newname =  vari+'_F'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['con']['var']:
									self.vardic['con']['var'].append([newvar,newname])
									self.vardic['con']['mod'].append(self.vardic['exo']['mod'][indx])
									self.audic['con']['var'].append([newvar,newname])
									self.audic['con']['mod'].append(self.vardic['exo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						elif vtype == 'exo' and int(vartime) < spvdic['exo'][indx][1][0]:
							spvdic2['exo'][indx][1][0] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_B'+tind+'(t)'
							newname =  vari+'_B'+str(abs(int(vartime)))
							list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp1[i1][poe:]
							for i2 in range(abs(int(vartime))):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_B'+tind+'(t)'
								newname =  vari+'_B'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['endo']['var']:
									self.vardic['endo']['var'].append([newvar,newname])
									self.vardic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
									self.audic['endo']['var'].append([newvar,newname])
									self.audic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						# Check for con
						if vtype == 'con' and int(vartime) > spvdic['con'][indx][1][1]:
							spvdic2['con'][indx][1][1] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_F'+tind+'(t)'
							newname =  vari+'_F'+str(abs(int(vartime)-1))
							list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp1[i1][poe:]
							for i2 in range(int(vartime)):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_F'+tind+'(t)'
								newname =  vari+'_F'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['con']['var']:
									self.vardic['con']['var'].append([newvar,newname])
									self.vardic['con']['mod'].append(self.vardic['con']['mod'][indx])
									self.audic['con']['var'].append([newvar,newname])
									self.audic['con']['mod'].append(self.vardic['con']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue
						elif vtype == 'con' and int(vartime) < spvdic['con'][indx][1][0]:
							spvdic2['con'][indx][1][0] = int(vartime)
							tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
							newvar = vari+'_B'+tind+'(t)'
							newname =  vari+'_B'+str(abs(int(vartime)))
							list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp1[i1][poe:]
							for i2 in range(abs(int(vartime))):
								tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
								newvar = vari+'_B'+tind+'(t)'
								newname =  vari+'_B'+str(abs(i2+1))
								if [newvar,newname] not in self.vardic['endo']['var']:
									self.vardic['endo']['var'].append([newvar,newname])
									self.vardic['endo']['mod'].append(self.vardic['con']['mod'][indx])
									self.audic['endo']['var'].append([newvar,newname])
									self.audic['endo']['mod'].append(self.vardic['con']['mod'][indx])
									if self.sstate.has_key(vari+'_bar'):
										self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
									else:
										self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
								continue


				# Now change the system to include possible augmented variables
				endo_r = filter(lambda x: x[1] != [-1,0], spvdic2['endo']) 
				if endo_r:
					endo_r = [[x[0],[abs(x[1][0]+1),x[1][1]]] for x in endo_r ]
					# Create lags and forwards equations
					for vari in endo_r:
						for lag in range(abs(vari[1][0])):
							tind = (5-len(str(lag+2)))*'0'+str(lag+2)
							if lag == 0:
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t-1)')
							else:
								tind1 = (5-len(str(lag+1)))*'0'+str(lag+1)
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
						for lead in range(vari[1][1]):
							tind = (5-len(str(lead)))*'0'+str(lead)
							if lead == 0:
								if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0] not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0])
							else:
								tind1 = (5-len(str(lead-1)))*'0'+str(lead-1)
								if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
				exo_r = filter(lambda x: x[1] != [0,1], spvdic2['exo'])
				if exo_r:
					exo_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in exo_r ]
					# Create lags and forwards equations
					for vari in exo_r:
						for lag in range(abs(vari[1][0])):
							tind = (5-len(str(lag+1)))*'0'+str(lag+1)
							if lag == 0:
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)')
							else:
								tind1 = (5-len(str(lag)))*'0'+str(lag)
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
						if vari[1][1] > 0:
							for lead in range(vari[1][1]+1):
								tind = (5-len(str(lead+1)))*'0'+str(lead+1)
								if lead == 0:
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
								else:
									tind1 = (5-len(str(lead)))*'0'+str(lead)
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
				con_r = filter(lambda x: x[1] != [0,1], spvdic2['con']) 
				if con_r:
					con_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in con_r ]
					# Create lags and forwards equations
					for vari in con_r:
						for lag in range(abs(vari[1][0])):
							tind = (5-len(str(lag+1)))*'0'+str(lag+1)
							if lag == 0:
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)')
							else:
								tind1 = (5-len(str(lag)))*'0'+str(lag)
								if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1: 
									list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
						if vari[1][1] > 0:
							for lead in range(int(vari[1][1])+1):
								tind = (5-len(str(lead+1)))*'0'+str(lead+1)
								if lead == 0:
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
								else:
									tind1 = (5-len(str(lead)))*'0'+str(lead)
									if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
										list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')

				return list_tmp1



			list_tmp1 = mkaug2(list_tmp1)
			list_tmp3 = [x[1] for x in list_tmp2]
			outtup = mkaug1(list_tmp3,list_tmp1)
				
			list_tmp3 = outtup[1]
			list_tmp1 = outtup[0]
			for i1,x in enumerate(list_tmp3):
				list_tmp2[i1][1] = list_tmp3[i1]

			self.nlsubs = dict(list_tmp2)
			nlsubs2 = {}
			for x in [x[0] for x in self.vardic['other']['var']]:
				nlsubs2[x] = self.nlsubs['@'+x]
			self.nlsubs2 = nlsubs2

			# Create ordered nlsubsys
			if self.vardic['other']['var']:
				nlsubsys = []
				varother = self.vardic['other']['var']
				for vari in [x[0] for x in varother]:
					nlsubsys.append(nlsubs2[vari])
				self.nlsubsys = nlsubsys

			# Count the number of distinct forward expectations
			# AND the no of equations that containt them
			ffli = []
			count = 0
			patup = ('0','endo|con|exo','{1,10}')
			for line in list_tmp1[:-self.nexo]:
				outli = self.vreg(patup,line,True,'min')
				if outli:
					count = count + 1
					for elem in [x[0] for x in outli]:
						if elem not in ffli:
							ffli.append(elem)
			self.ffen = count
			self.ffli = ffli
			self.fflin = len(ffli)

			# Now create simplified vdic with possibly new augmented vars
			vdic = dict([[x,self.vardic[x]['var']] for x in self.vardic.keys()])
			self.vdic = vdic

			# Finally, recomputed nendo, ncon, etc.
			self.nendo = len(self.vardic['endo']['var'])
			self.nexo = len(self.vardic['exo']['var'])
			self.ncon = len(self.vardic['con']['var'])
			self.nother = len(self.vardic['other']['var'])
			self.nstat = self.nendo+self.nexo
			self.nall = self.nstat+self.ncon


			return list_tmp1
		if sum([nreg.search(x)!=None for x in input.secs['focs'][0]]) == 0:
			self.nlsys_list = mknlinsys()

		# Creates the log-linear system
		def mkloglinsys1():
			list_tmp = []
			i1=0
			counter=0
			for x in input.secs['modeq'][0]:
				if '...' in x:
					counter = counter + 1
				elif '...' not in x:
					if counter == 0:
						list_tmp.append(x.replace(';','').split(']')[1].strip())
					elif counter > 0:
						str_tmp = ''
						for y in input.secs['modeq'][0][i1-counter:i1+1]:
							str_tmp = str_tmp + y.replace('...','').replace(';','').strip()
						list_tmp.append(str_tmp.split(']')[1].strip())
						counter = 0	
				i1=i1+1
			return list_tmp

		def mkloglinsys2(inlist):
			"""
			This function simply takes the left-hand side
			of each equation and takes it over to the
			right-hand side with a minus sign in front.
			Now we have equations of the form g(x)=0.
			Also takes care of fact that the term may
			already be negative. Also takes into account
			that some equations have already been written
			as g(x)=0!
			"""
			_revar='(E\(t-{0,1}\+{0,1}\d*\)\|){0,1}[a-z]*(\(t-{0,1}\+{0,1}\d*\)){1,1}'
			re_var = RE.compile(_revar)
			list_tmp1 = COP.deepcopy(inlist)
			str_tmp2 = ''
			i1 = 0
			while i1 < len(list_tmp1):

				if list_tmp1[i1].split('=')[0].strip().strip(';') == '0':
					list_tmp1[i1] = list_tmp1[i1].split('=')[1].strip().strip(';')
					i1 = i1 + 1
					continue
				elif list_tmp1[i1].split('=')[1].strip().strip(';') == '0':
					list_tmp1[i1] = list_tmp1[i1].split('=')[0].strip().strip(';')
					i1 = i1 + 1
					continue

				list_tmp1[i1] = list_tmp1[i1].strip(';')[:]
				str_tmp2 = list_tmp1[i1].split('=')[0].strip()
				list_tmp1[i1] = list_tmp1[i1].split('=')[1].strip()

				while re_var.search(str_tmp2):
					ma1 = re_var.search(str_tmp2)
					pos1 = ma1.span()[0]
					poe1 = ma1.span()[1]
					# for coeff*var and coeff not begins with '-', but '+'
					if pos1 > 1 and str_tmp2[0] != '-' and str_tmp2[0] == '+':
						coeff = '-'+str_tmp2[1:pos1]
						vari = ma1.group(0)
						str_tmp2=str_tmp2[poe1:]
						if str_tmp2 == '':
							str_tmp2 = '@'
					# for coeff*var and coeff not begins with '+', but '-'		
					elif pos1 > 1 and str_tmp2[0] != '+' and str_tmp2[0] == '-':
						coeff = '+'+str_tmp2[1:pos1]
						vari = ma1.group(0)
						str_tmp2=str_tmp2[poe1:]
						if str_tmp2 == '':
							str_tmp2 = '@'
					# for coeff*var and coeff not begins with '+' and not begins with '-'
					elif pos1 > 1 and str_tmp2[0] != '+' and str_tmp2[0] != '-':
						coeff = '-'+str_tmp2[:pos1]
						vari = ma1.group(0)
						str_tmp2=str_tmp2[poe1:]
						if str_tmp2 == '':
							str_tmp2 = '@'
					# for coeff == '-'
					elif pos1 == 1 and str_tmp2[0] == '-':
						coeff = '+'
						vari = ma1.group(0)
						str_tmp2=str_tmp2[poe1:]
					# for coeff == '+'
					elif pos1 == 1 and str_tmp2[0] == '+':
						coeff = '-'
						vari = ma1.group(0)
						str_tmp2=str_tmp2[poe1:]
					# for var*coeff
					elif pos1 == 0 and len(re_var.findall(str_tmp2)) > 1\
					     and str_tmp2[poe1] != '-' and str_tmp2[poe1] != '+':
						ma2 = re_var.search(str_tmp2[poe1:])
						pos2 = ma2.span()[0]+poe1
						coeff = '-'+str_tmp2[poe1:pos2]
						vari = ma1.group(0)
						str_tmp2=str_tmp2[pos2:]
					# for last bit
					elif pos1 == 0 and len(re_var.findall(str_tmp2)) == 1:
						coeff = '-'+str_tmp2[1:][poe1-1:]
						vari = ma1.group(0)
						str_tmp2=''
					# for coeff == ''
					elif pos1 == 0 and len(re_var.findall(str_tmp2)) > 1 and str_tmp2[poe1] == '-':
						coeff = '-'
						vari = ma1.group(0)
						str_tmp2=str_tmp2[poe1:]
					# for coeff == ''
					elif pos1 == 0 and len(re_var.findall(str_tmp2)) > 1 and str_tmp2[poe1] == '+':
						coeff = '-'
						vari = ma1.group(0)
						str_tmp2=str_tmp2[poe1:]
					# for last bit
					elif pos1 == 0 and len(re_var.findall(str_tmp2)) == 1:
						coeff = '-'
						vari = ma1.group(0)
						str_tmp2=''

					if coeff[-1] != '*' and coeff != '-':
						coeff = coeff+'*'
					list_tmp1[i1] = list_tmp1[i1]+coeff+vari
				str_tmp2 = ''	
				i1 = i1 + 1
			return list_tmp1

		if sum([nreg.search(x)!=None for x in input.secs['modeq'][0]]) == 0:
			self.llsys_list = mkloglinsys1()
			self.llsys_list = mkloglinsys2(self.llsys_list)

		# Make symbolic system and numeric as well
		def mksymsys():
			func = []
			subli = []
			patup = ('{-10,10}|None','all','{-10,10}')
			for x in self.sstate.keys()+self.paramdic.keys():
				locals()[x] = eval("SP.Symbol('"+x+"')")
			for y in self.llsys_list:
				str_tmp = y[:]
				vali = [x[0] for x in self.vreg(patup,y,True,'min')]
				vali2 = [[x,'sub'+str(u)] for x,u in zip(vali,range(0,len(vali),1))]
				vali3 = [(x[0],x[3]) for x in self.vreg(patup,str_tmp,True,'max')]
				vali3.reverse()
				valdic = dict(vali2)
				for x in vali3:
					str_tmp = str_tmp[:x[1][0]] + valdic[x[0]] + str_tmp[x[1][1]:]
				subli.append(valdic)
				for x in valdic.values():
					locals()[x] = eval("SP.Symbol('"+x+"')")
				func.append(eval(str_tmp))
			diffli = []
			i1 = 0
			for x in func:
				diffdic = {}
				for y in subli[i1].items():
					diffdic[y[0]] = eval('SP.diff(func[i1],'+y[1]+')')
				diffli.append(diffdic)
				i1=i1+1
			self.diffli1 = diffli
			def evaldi():
				diffli2 = []
				locals().update(self.sstate)
				locals().update(self.paramdic)
				for x in self.diffli1:
					tmpdic={}
					for y in x.keys():
						tmpdic[y] = eval(x[y].tostr())
					diffli2.append(tmpdic) 
				self.diffli2 = diffli2
			evaldi()
		if sum([nreg.search(x)!=None for x in input.secs['modeq'][0]]) == 0:
			mksymsys()
		# Create a Variance-Covariance matrix
		def mksigmat():
			str_tmp1 = []
			str_tmp1 = input.secs['vcvm'][0][0].split('[')[1].rstrip()
			if str_tmp1[-2:] == '];':
				str_tmp1=str_tmp1[:-2]
			if len(input.secs['vcvm']) > 1:
				for x in input.secs['vcvm'][0][1:]:
					if x[-2:] != '];':
						str_tmp1=str_tmp1+x.lstrip().rstrip()[:]
					elif x[-2:] == '];':
						str_tmp1=str_tmp1+x.lstrip().rstrip()[:-2]
			locals().update(self.paramdic)
			locals().update(self.sstate)
			list_tmp1 = str_tmp1.split(';')
			len1 = len(list_tmp1[0].split())
			mat1=MAT.zeros((len1,len1))
			i1=0
			for x in list_tmp1:
				i2=0
				for y in x.split():
					mat1[i1,i2] = float(eval(y))
					i2=i2+1
				i1=i1+1
			return mat1	
		if sum([nreg.search(x)!=None for x in input.secs['vcvm'][0]]) == 0 and 'sstate' in dir(self):
			self.sigma = mksigmat()

		def mkeqtype():
			lsys = self.llsys_list
			tup1 = ('{-1,1}|None','iid|exo','{-1,1}')
			tup2 = ('{-1,1}|None','all','{-1,1}')
			tup3 = ('{-1,1}','all','{-1,1}')
			err_indx = []
			det_indx = []
			exp_indx = []
			for x,y in zip(lsys,range(len(lsys))):
				if self.vreg(('{-1,1}','all','{-1,1}'),x,False,'max'):
					exp_indx.append(y)
				elif self.vreg((None,'exo|iid','{-1,1}'),x,False,'max'):
					if len(self.vreg((None,'exo|iid','{-1,1}'),x,True,'max'))==\
					   len(self.vreg(('{-1,1}|None','all','{-1,1}'),x,True,'max')):
						err_indx.append(y)
					else:
						det_indx.append(y)
				else:
					det_indx.append(y)


			self.eqindx = {}
			self.eqindx['err'] = err_indx
			self.eqindx['det'] = det_indx
			self.eqindx['exp'] = exp_indx
		if sum([nreg.search(x)!=None for x in input.secs['modeq'][0]]) == 0:
			mkeqtype()			

	# html info of model opened with webbrowser
	def info(self):
		tmplist = glob.glob('tempz*.html')
		for x in tmplist:
			OPS.remove(x)
		modname = self.modname
		secs = self.txtpars.secs
		direc = OPS.getcwd()
		fd,fpath = TMF.mkstemp(prefix='tempz',suffix='.html',dir=direc)
		htmlfile = OPS.fdopen(fd,'w+b')
		htmlfile.write('<HTML><BODY BGCOLOR="white">\n')
		htmlfile.write('<H2>%s</H2>\n'%modname)
		htmlfile.write('\n\n')
		for x1 in secs['mod'][1]:
			htmlfile.write('<P>'+x1+'\n')
		htmlfile.write('<P>'+'<H4>Model Parameters</H4>\n')
		for x1 in secs['para'][1]:
			htmlfile.write('<P>'+x1+'\n')
		htmlfile.write('<P>'+'<H4>Nonlinear First-Order Conditions</H4>\n')
		for x1 in secs['focs'][1]:
			if x1[0] == '[':
				htmlfile.write('<P>'+x1+'\n')
			elif x1[0] != '[':
				htmlfile.write('<P>'+4*'&nbsp;'+ x1+'\n')
		htmlfile.write('<P>'+'<H4>Log-Linearized Model Equations</H4>\n')
		for x1 in secs['modeq'][1]:
			if x1[0] == '[':
				htmlfile.write('<P>'+x1+'\n')
			elif x1[0] != '[':
				htmlfile.write('<P>'+4*'&nbsp;'+ x1+'\n')
		htmlfile.write('</BODY></HTML>\n')
		htmlfile.write('<P>'+'<H2>Solution Section</H2>\n')
		htmlfile.write('<P>'+'<H4>Uhlig Toolkit Output</H4>\n')
		htmlfile.close()
		cmd = browserpath+' '+fpath
		OPS.system(cmd)
		tmplist = glob.glob('tempz*.html')
		for x in tmplist:
			OPS.remove(x)
		return 'Model Website opened!'
	# Set the author of the current model
	def setauthor(self,author=None):
		if author == None:
			self.author = 'No author'
		else:
			self.author = author
	# pdflatex the model tex file and open with pdf viewer
	def pdf(self):
		modfile = self.modfile
		modfname = modfile.split('.')[0]
		if modfname+'.pdf' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.pdf'))
		if modfname+'.log' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.log'))
		if modfname+'.aux' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.aux'))
		if modfname+'.dvi' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.dvi'))
		if modfname+'.log' in OPS.listdir(OPS.getcwd()):
			OPS.remove(OPS.path.join(OPS.getcwd(),modfname+'.log'))

		# Does the model tex file even exist?
		if modfname+'.tex' not in OPS.listdir(modfpath):
			print 'Error: The model tex file does not exist!'
			print 'Use model.texed() to create a new one!'
			return

		# First check for Chktex syntax errors
		if 'texer.log' in OPS.listdir(OPS.getcwd()):
			OPS.remove(OPS.path.join(OPS.getcwd(),'texer.log'))
		args = '--quiet -v2 -n3 -n25 -n12 -n35'
		cmd = 'chktex '+args+' '+OPS.path.join(modfpath,modfname+'.tex'+' > texer.log')
		OPS.system(cmd)
		file = open(OPS.path.join(OPS.getcwd(),'texer.log'),'rU')
		errlist = file.readlines()
		if len(errlist)== 0:
			file.close()
			OPS.remove(OPS.path.join(OPS.getcwd(),'texer.log'))
		else:
			print 'There were errors in the model tex file! Abort!\n'

			for x in errlist:
				print x
			file.close()
			OPS.remove(OPS.path.join(OPS.getcwd(),'texer.log'))
			return

		# PDFLatex the tex file and open
		cmd = 'pdflatex '+'-output-directory='+modfpath+' '+OPS.path.join(modfpath,modfname+'.tex'+' > out.log')
		OPS.system(cmd)
		cmd2 = pdfpath + ' ' + OPS.path.join(modfpath,modfname+'.pdf')
		OPS.system(cmd2)
		if modfname+'.pdf' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.pdf'))
		if modfname+'.log' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.log'))
		if modfname+'.aux' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.aux'))
		if modfname+'.dvi' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.dvi'))
		if modfname+'.log' in OPS.listdir(OPS.getcwd()):
			OPS.remove(OPS.path.join(OPS.getcwd(),modfname+'.log'))
	# Opens the model text file for editing
	def txted(self):
		modfile = self.modfile
		modfname = modfile.split('.')[0]
		if modfname+'.log' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.log'))
		if modfname+'.log' in OPS.listdir(OPS.getcwd()):
			OPS.remove(OPS.path.join(OPS.getcwd(),modfname+'.log'))
		cmd = txtedpath+' '+OPS.path.join(modfpath,self.modfile+' > out.log')
		timestamp0 = OPS.stat(OPS.path.join(modfpath,self.modfile))[8]
		OPS.system(cmd)
		timestamp1 = OPS.stat(OPS.path.join(modfpath,self.modfile))[8]
		if timestamp0 != timestamp1:
			self.__init__(self.modfile,self.dbase)
			self.init2()
		if modfname+'.log' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.log'))
		if modfname+'.log' in OPS.listdir(OPS.getcwd()):
			OPS.remove(OPS.path.join(OPS.getcwd(),modfname+'.log'))
	# Opens the model latex file for editing
	def texed(self):
		modfile = self.modfile
		modfname = modfile.split('.')[0]
		if modfname+'.log' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.log'))
		if modfname+'.log' in OPS.listdir(OPS.getcwd()):
			OPS.remove(OPS.path.join(OPS.getcwd(),modfname+'.log'))
		if modfname+'.tex' not in OPS.listdir(modfpath):
			file = open(OPS.path.join(modfpath,modfname+'.tex'),'wb')
			file.write('\\documentclass[a4paper,11pt]{article}\n')
			file.write('\\title{%s}\n'% self.modname)
			file.write('\\author{%s}\n'% self.author)
			file.write('\\begin{document}\n')
			file.write('\\maketitle\n\n')
			file.write('Start writing your model in Latex here!')
			file.write(8*'\n')
			file.write('\\end{document}\n')
			file.close()
		cmd = texedpath+' '+OPS.path.join(modfpath,self.modfile[:-3]+'tex'+' > out.log')
		OPS.system(cmd)
		if modfname+'.log' in OPS.listdir(modfpath):
			OPS.remove(OPS.path.join(modfpath,modfname+'.log'))
		if modfname+'.log' in OPS.listdir(OPS.getcwd()):
			OPS.remove(OPS.path.join(OPS.getcwd(),modfname+'.log'))
	# Choose to delete the current model's latex file. Your are asked again whether yes/no
	def deltex(self):
		answ = 0
		while answ not in ['Y','y','N','n']:
			answ = raw_input("REALLY delete this model's tex file ? (Y/N)")
		if answ in ['Y','y']:
			modfile = self.modfile
			modfname = modfile.split('.')[0]
			if modfname+'.tex' in OPS.listdir(modfpath):
				OPS.remove(OPS.path.join(modfpath,modfname+'.tex'))
		elif answ in ['N','n']:
			return
	# Change the 'current view' of model solution algorithms
	def ccv(self,instring=None):
		if instring == None:
			return 'Error: Your have not specified a current view name!'
		elif instring not in dir(self.modsolvers):
			return 'Error: Your chosen view does not exist for this model!'
		else:
			self.cv = eval('self.modsolvers.'+instring)
	# Get data method, in case model has not been loaded with database object
	def getdata(self,datafile=None,dbase=None):
		self.data = {}
		if datafile != None:
			input = open(OPS.path.join(OPS.getcwd(),datafile), 'r')
			output = open(OPS.path.join(OPS.getcwd(),'tmp_001.csv'), 'w')
			wholefile = input.read()
			lines = wholefile.splitlines()
			varnames = lines[0]
			varnlist = varnames.split(',')[:]
			varnlist = [x.replace('"','') for x in varnlist][:]
			str_tmp=''
			for x1 in lines[1:]:
				str_tmp = str_tmp + x1+'\n'
			output.write(str_tmp)
			output.close()
			rawdata = P.load('tmp_001.csv',delimiter=',')
			OPS.remove(OPS.getcwd()+'/'+'tmp_001.csv')
			self.rawdata = rawdata
			for x1 in self.dataprop.items():
				if type(x1[1]) != type('x'):
					exec(x1[0]+'='+str(x1[1]))
				elif type(x1[1]) == type('x'):
					exec(x1[0]+'='+"'"+x1[1]+"'")

			if N.shape(rawdata)[1] != len(self.varnames):
				print 'DAT IMPORT ERR: Len(varnames) != Len(data)'
				data = TimeSeries(rawdata,start_date=TSS.Date(freq=freq,year=year,month=month))
			else:
				data = TimeSeries(rawdata,start_date=TSS.Date(freq=freq,year=year,month=month))
			self.data = TimeSeriesCol(varnlist,data)
		elif dbase != None:
			varnames = []
			for x in self.vardic.keys():
				for y in self.vardic[x]['var']:
					varnames.append(y[1])
			for x1 in varnames:
				if dbase.modif.has_key(x1):
					self.data[x1] = dbase.modif[x1]
			self.VAR = dbase.VAR
		elif datafile == None and dbase == None:
			self.data = None

	# The regex function
	def vreg(self,paratuple=(None,'all','0'),cinstring='',iter=False,info='min'):
		'''
		The regex function for variable detection out of strings
		'''
		vars = {}
		vars['con'] = []
		vars['endo'] = []
		vars['exo'] = []
		vars['other'] = []
		vars['iid'] = []
		vars['state'] = []
		vars['all'] = []
		instring = COP.deepcopy(cinstring)

		for x in self.vardic['con']['var']:
			vars['con'].append(x[0].split('(')[0])
		for x in self.vardic['endo']['var']:
			vars['endo'].append(x[0].split('(')[0])
		for x in self.vardic['exo']['var']:
			vars['exo'].append(x[0].split('(')[0])
			vars['iid'].append(x[2].split('(')[0])
		for x in self.vardic['other']['var']:
			if '|' in x[0]:
				vars['other'].append(x[0].split('|')[1].split('(')[0])
			else:
				vars['other'].append(x[0].split('(')[0])
		vars['state'] = vars['endo'] + vars['exo']
		vars['all'] = vars['state'] + vars['con'] + vars['iid'] + vars['other']

		# Create regex part for expectation info set
		if paratuple[0] != None and '{' not in paratuple[0]:
			if paratuple[0] == '0':
				iti = '(E\(t(?P<itime>.{0,0})\)\|)'
			elif paratuple[0][0] == '-':
				iti = '(E\(t(?P<itime>-'+paratuple[0][1:]+')\)\|)'
			else:
				iti = '(E\(t(?P<itime>\+'+paratuple[0]+')\)\|)'
		elif paratuple[0] != None and '{' in paratuple[0]:
			if '|' not in paratuple[0]:
				indx1 = paratuple[0].split(',')[0][1:]
				indx2 = paratuple[0].split(',')[1][:-1]
				indx2 = str(eval(indx2+'+1'))
			elif '|' in paratuple[0]:
				strtmp = paratuple[0].split('|')[0]
				indx1 = strtmp.split(',')[0][1:]
				indx2 = strtmp.split(',')[1][:-1]
				indx2 = str(eval(indx2+'+1'))
			itrange = eval('range('+indx1+','+indx2+',1)')
			zswitch = 0
			if 0 in itrange:
				zswitch = 1
				itrange.pop(itrange.index(0))		
			str_tmp=''
			for x in itrange:
				if str(x)[0] == '-':
					str_tmp = str_tmp + '|'+str(x)
				elif str(x) == '0':
					str_tmp = str_tmp + '|'+str(x)
				else:
					str_tmp = str_tmp + '|\+'+str(x)
			str_tmp = str_tmp[1:]+zswitch*'|.{0,0}'
			orswitch = False
			if '|' in paratuple[0]:
				orswitch = True
			iti = '(E\(t(?P<itime>'+str_tmp[:]+')\)\|)'+(1-orswitch)*'{1,1}'+orswitch*'{0,1}'
		elif paratuple[0] == None:
			iti = ''


		# Create regex part for variable type
		if '|' not in paratuple[1]:
			varp = ''
			for x in vars[paratuple[1]]:
				varp = varp + '|'+x
		elif '|' in paratuple[1]:
			typeli = paratuple[1].split('|')
			varp = ''
			for x in typeli:
				for y in vars[x]:
					varp = varp + '|'+y
		varp = varp[1:]
		varp = '(?P<var>'+varp+')'


		# Create the time part
		str_tmp = ''
		if '{' not in paratuple[2]:
			if paratuple[2][0] == '-':
				ti = '\(t(?P<time>'+paratuple[2]+')\)'
			elif paratuple[2] == '0':
				ti = '\(t(?P<time>.{0,0})\)'
			else:
				ti = '\(t(?P<time>\+'+paratuple[2]+')\)'
		elif '{' in paratuple[2]:
			indx1 = paratuple[2].split(',')[0][1:]
			indx2 = paratuple[2].split(',')[1][:-1]
			indx2 = str(eval(indx2+'+1'))
			trange = eval('range('+indx1+','+indx2+',1)')
			zswitch = False
			if 0 in trange:
				zswitch = True
				trange.pop(trange.index(0))
			str_tmp=''
			for x in trange:
				if str(x)[0] == '-':
					str_tmp = str_tmp + '|'+str(x)
				elif str(x) == '0':
					str_tmp = str_tmp + '|'+str(x)
				else:
					str_tmp = str_tmp + '|\+'+str(x)
			str_tmp = str_tmp[1:]+zswitch*'|.{0,0}'
			ti = '(\(t(?P<time>'+str_tmp[:]+')\))'

		reva = RE.compile('(?P<pre>^[^a-zA-Z^_]{0,1}|[^a-zA-Z^_]{1,1})('+iti+varp+ti+')')

		if reva.search(instring) and iter == False:
			ma = reva.search(instring)
			if ma.group('pre') != None:
				pre = ma.group('pre')
				prel = len(pre)
			else:
				prel = 0
			alcap = ma.group()[prel:]
			var = ma.group('var')
			time = ma.group('time')
			if time != '':
				time = ma.group('time')
			else:
				time = '0'
			span = (ma.span()[0]+prel,ma.span()[1])
			vtype = (var in vars['endo'])*'endo'+\
			      (var in vars['exo'])*'exo'+\
			      (var in vars['con'])*'con'+\
			      (var in vars['iid'])*'iid'+\
			      (var in vars['other'])*'other'
			if vtype == 'endo':
				posx = [x[0].split('(')[0] for x in self.vardic['endo']['var']].index(var)
			elif vtype == 'exo':
				posx = [x[0].split('(')[0] for x in self.vardic['exo']['var']].index(var)
			elif vtype == 'con':
				posx = [x[0].split('(')[0] for x in self.vardic['con']['var']].index(var)
			elif vtype == 'iid':
				posx = [x[2].split('(')[0] for x in self.vardic['exo']['var']].index(var)
			elif vtype == 'other':
				posx = [x[0] for x in self.vardic['other']['var']].index(alcap)

			if ma.groupdict().has_key('itime'):
				itime = ma.group('itime')
				if itime != '':
					itime = ma.group('itime')
				else:
					itime = '0'
			else:
				itime = None

			if info == 'max':
				return (alcap,(vtype,posx),(itime,var,time),span)
			elif info == 'min':
				return (alcap,(vtype,posx),(itime,var,time))
		elif reva.search(instring) and iter == True:
			itobj = reva.finditer(instring)
			relist=[]
			for ma in itobj:
				if ma.group('pre') != None:
					pre = ma.group('pre')
					prel = len(pre)
				else:
					prel = 0
				alcap = ma.group()[prel:]
				var = ma.group('var')
				time = ma.group('time')
				if time != '':
					time = ma.group('time')
				else:
					time = '0'
				span = (ma.span()[0]+prel,ma.span()[1])
				vtype = (var in vars['endo'])*'endo'+\
				      (var in vars['exo'])*'exo'+\
				      (var in vars['con'])*'con'+\
				      (var in vars['iid'])*'iid'+\
				      (var in vars['other'])*'other'
				if vtype == 'endo':
					posx = [x[0].split('(')[0] for x in self.vardic['endo']['var']].index(var)
				elif vtype == 'exo':
					posx = [x[0].split('(')[0] for x in self.vardic['exo']['var']].index(var)
				elif vtype == 'con':
					posx = [x[0].split('(')[0] for x in self.vardic['con']['var']].index(var)
				elif vtype == 'iid':
					posx = [x[2].split('(')[0] for x in self.vardic['exo']['var']].index(var)
				elif vtype == 'other':
					posx = [x[0] for x in self.vardic['other']['var']].index(alcap)
				if ma.groupdict().has_key('itime'):
					itime = ma.group('itime')
					if itime != '':
						itime = ma.group('itime')
					else:
						itime = '0'
				else:
					itime = None

				relist.append((alcap,(vtype,posx),(itime,var,time),span))
			if info == 'min':
				relist2 = [(x[0],x[1],x[2]) for x in relist]
				for x in relist2:
					while relist2.count(x) > 1:
						indx = relist2.index(x)
						relist2.pop(indx)
				return relist2
			elif info == 'max':
				return relist
		else:
			return False
###########ANALYTIC AND NUMERICAL JACOBIAN AND HESSIAN METHODS############
	# The unparallelized Jacobian, Hessian computation method
	def mkjahe(self):
		'''
		This produces the Jacobian for the non-linear system of equations
		around the steady state, and evaluated at the steady state.
		Also produces the Hessian!
		'''
		exo_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['exo']['var']]
		endo_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['endo']['var']]
		con_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['con']['var']]
		exo_0 = [x[0] for x in self.vardic['exo']['var']]
		endo_0 = [x[0].split('(')[0]+'(t-1)' for x in self.vardic['endo']['var']]
		con_0 = [x[0] for x in self.vardic['con']['var']]
		inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
		intup=tuple(inlist)

		patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
		jcols = len(intup)
		jrows = len(self.nlsys_list)
		nlsys = COP.deepcopy(self.nlsys_list)

		if 'nlsubsys' in dir(self):
			lsubs = len(self.nlsubsys)
			jrows = jrows + lsubs
			nlsys = nlsys + COP.deepcopy(self.nlsubsys)

		# Create substitution var list and dictionary
		tmpli = []
		for i,x in enumerate(intup):
			nolen = len(str(i))
			inds = (5-nolen)*'0'+str(i)
			tmpli.append([x,'SUB'+inds])
			dicli = dict(tmpli)
			dicli2 = dict([[x[1],x[0]] for x in tmpli])

		func = []
		subli = []
		symdic = {}
		for x in dicli.values():
			symdic[x] = SPC.Symbol(x)
		locals().update(symdic)
		for x in self.paramdic.keys():
			locals()[x] = SPC.Symbol(x)
		for x in self.sstate.keys():
			locals()[x] = SPC.Symbol(x)
		for x in nlsys:
			str_tmp = x[:]
			list1 = self.vreg(patup,x,True,'max')
			list1.reverse()
			for y in list1:
				pos = y[3][0]
				poe = y[3][1]
				vari = y[0]
				str_tmp = str_tmp[:pos] + dicli[vari] +str_tmp[poe:]
			list2 = self.vreg(('{-10,10}|None','iid','{-10,10}'),str_tmp,True,'max')
			if list2:
				list2.reverse()
				for y in list2:
					pos = y[3][0]
					poe = y[3][1]
					vari = y[0]
					str_tmp = str_tmp[:pos]+'0'+str_tmp[poe:]
			# Now substitute out exp and log in terms of sympycore expressions
			elog = RE.compile('LOG\(')
			while elog.search(str_tmp):
				ma = elog.search(str_tmp)
				pos = ma.span()[0]
				poe = ma.span()[1]
				str_tmp = str_tmp[:pos]+'SPC.log('+str_tmp[poe:]
			eexp = RE.compile('EXP\(')
			while eexp.search(str_tmp):
				ma = eexp.search(str_tmp)
				pos = ma.span()[0]
				poe = ma.span()[1]
				str_tmp = str_tmp[:pos]+'SPC.exp('+str_tmp[poe:]
			func.append(eval(str_tmp))
		self.func1 = func

		# Only exp() when variable needs to be put into ln !
		list1 = [x[1] for x in self.vardic['endo']['var']]
		list1 = list1 + [x[1] for x in self.vardic['con']['var']]
		list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
		list2 = [x for x in self.vardic['endo']['mod']]
		list2 = list2 + [x for x in self.vardic['con']['mod']]
		list2 = list2 + [x for x in self.vardic['exo']['mod']]
		if self.vardic['other']['var']:
			list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
			list2 = list2 +  [x for x in self.vardic['other']['mod']]
		moddic = {}
		for x,y in zip(list1,list2):
			moddic[x] = y
			moddic[x] = y

		lookli = []
		for key in self.vardic.keys():
			lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
		lookdic = dict(lookli)
		_mreg = 'SUB\d{5,5}'
		mreg = RE.compile(_mreg)
		func2 = []
		for i,x in enumerate(func):
			func2.append(func[i])
			liobj = list(mreg.finditer(str(func2[i])))
			doneli = []
			for ma in reversed(liobj):
				suba = ma.group()
				if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
				   and suba not in doneli:
					func2[i] = func2[i].subs(SPC.Symbol(suba),SPC.exp(SPC.Symbol(suba)))
					doneli.append(suba)
		self.func2 = func2

		# Also take into account log or not for other equations !
		if 'nlsubsys' in dir(self):
			tmpvarli0 = [x[0] for x in self.vardic['other']['var']]
			for i1,line in enumerate(self.func2[-lsubs:]):
				vari0 = tmpvarli0[i1]
				if 'log' in self.vardic['other']['mod'][i1]:
					pass
				else:
					vma = self.vreg(patup,vari0,False,'min')
					vmavar = vma[2][1]+'_bar'
					self.func2[i1] = SPC.Symbol(vmavar)*(self.func2[i1])

		# Create list with log of steady states
		evalli = []
		alldic = COP.deepcopy(self.paramdic)
		alldic.update(self.sstate)
		for x in intup:
			vma = self.vreg(patup, x, False, 'min')
			vari = vma[2][1]
			if 'log' in moddic[lookdic[vari]]:
				evalli.append(N.log(alldic[vari+'_bar']))
			else:
				evalli.append(alldic[vari+'_bar'])
		evaldic = {}
		for i,x in enumerate(tmpli):
			evaldic[x[1]] = evalli[i]

		# Make 2D symbolic and numeric Jacobian
		def mkjac(jrows=jrows,jcols=jcols):
			rdic = dict([[x,'0'] for x in range(jcols)])
			jdic = dict([[x,COP.deepcopy(rdic)] for x in range(jrows)])
			jcols = len(jdic[0])
			jrows = len(jdic)
			numj = MAT.zeros((jrows,jcols))
			alldic = {}
			alldic.update(self.paramdic)
			alldic.update(self.sstate)
			alldic.update(evaldic)
			locals().update(alldic)
			for x in range(jrows):
				for y in range(jcols):
					jdic[x][y] = func2[x].diff(symdic[tmpli[y][1]])
					numj[x,y] = eval(str(jdic[x][y].evalf()))
			return numj,jdic

		# Now make 3D symbolic and numeric Hessian
		def mkhes(jrows=jrows,jcols=jcols):
			rdic = dict([[x,'0'] for x in range(jcols)])
			rdic1 = dict([[x,COP.deepcopy(rdic)] for x in range(jcols)])
			hdic = dict([[x,COP.deepcopy(rdic1)] for x in range(jrows)])
			hcols = len(hdic[0])
			hrows = len(hdic[0])*len(hdic)
			numh = MAT.zeros((hrows,hcols))
			jdic = self.jdic
			alldic = {}
			alldic.update(self.paramdic)
			alldic.update(self.sstate)
			alldic.update(evaldic)
			locals().update(alldic)
			count = 0
			for x in range(jrows):
				for y in range(jcols):
					for z in range(jcols):
						hdic[x][y][z] = jdic[x][y].diff(symdic[tmpli[z][1]])
						numh[count,z] = eval(str(hdic[x][y][z].evalf()))
					count = count + 1
			return numh,hdic

		self.numj,self.jdic = mkjac()
		numj = self.numj
		if mk_hessian:
			self.numh,self.hdic = mkhes()
			numh = self.numh

		if 'nlsubsys' in dir(self):
			numjs = numj[:-lsubs,:]
			numjl = numj[-lsubs:,:]
			self.numj = numjs
			self.numjl = numjl
			if mk_hessian:
				numhs = numh[:-lsubs*jcols,:]
				numhl = numh[-lsubs*jcols:,:]
				self.numhl = numhl
				self.numh = numhs
		else:
			self.numj = numj
			if mk_hessian:
				self.numh = numh

		self.jAA = self.numj[:,:int(len(intup)/2)]
		self.jBB = -self.numj[:,int(len(intup)/2):]
	# The parallelized mkjahe version using parallel python
	def mkjahepp(self):
		'''
		This produces the Jacobian for the non-linear system of equations
		around the steady state, and evaluated at the steady state.
		Also produces the Hessian!
		'''
		# import local sympycore
		import sympycore

		exo_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['exo']['var']]
		endo_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['endo']['var']]
		con_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['con']['var']]
		exo_0 = [x[0] for x in self.vardic['exo']['var']]
		endo_0 = [x[0].split('(')[0]+'(t-1)' for x in self.vardic['endo']['var']]
		con_0 = [x[0] for x in self.vardic['con']['var']]
		inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
		intup=tuple(inlist)

		patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
		jcols = len(intup)
		jrows = len(self.nlsys_list)
		nlsys = COP.deepcopy(self.nlsys_list)

		if 'nlsubsys' in dir(self):
			lsubs = len(self.nlsubsys)
			jrows = jrows + lsubs
			nlsys = nlsys + COP.deepcopy(self.nlsubsys)

		# Create substitution var list and dictionary
		tmpli = []
		for i,x in enumerate(intup):
			nolen = len(str(i))
			inds = (5-nolen)*'0'+str(i)
			tmpli.append([x,'SUB'+inds])
			dicli = dict(tmpli)
			dicli2 = dict([[x[1],x[0]] for x in tmpli])

		func = []
		subli = []
		symdic = {}
		for x in dicli.values():
			symdic[x] = sympycore.Symbol(x)
		locals().update(symdic)
		for x in self.paramdic.keys():
			locals()[x] = sympycore.Symbol(x)
		for x in self.sstate.keys():
			locals()[x] = sympycore.Symbol(x)
		for x in nlsys:
			str_tmp = x[:]
			list1 = self.vreg(patup,x,True,'max')
			list1.reverse()
			for y in list1:
				pos = y[3][0]
				poe = y[3][1]
				vari = y[0]
				str_tmp = str_tmp[:pos] + dicli[vari] +str_tmp[poe:]
			list2 = self.vreg(('{-10,10}|None','iid','{-10,10}'),str_tmp,True,'max')
			if list2:
				list2.reverse()
				for y in list2:
					pos = y[3][0]
					poe = y[3][1]
					vari = y[0]
					str_tmp = str_tmp[:pos]+'0'+str_tmp[poe:]
			# Now substitute out exp and log in terms of sympycore expressions
			elog = RE.compile('LOG\(')
			while elog.search(str_tmp):
				ma = elog.search(str_tmp)
				pos = ma.span()[0]
				poe = ma.span()[1]
				str_tmp = str_tmp[:pos]+'sympycore.log('+str_tmp[poe:]
			eexp = RE.compile('EXP\(')
			while eexp.search(str_tmp):
				ma = eexp.search(str_tmp)
				pos = ma.span()[0]
				poe = ma.span()[1]
				str_tmp = str_tmp[:pos]+'sympycore.exp('+str_tmp[poe:]
			func.append(eval(str_tmp))
		self.func1 = func

		# Only exp() when variable needs to be put into ln !
		list1 = [x[1] for x in self.vardic['endo']['var']]
		list1 = list1 + [x[1] for x in self.vardic['con']['var']]
		list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
		list2 = [x for x in self.vardic['endo']['mod']]
		list2 = list2 + [x for x in self.vardic['con']['mod']]
		list2 = list2 + [x for x in self.vardic['exo']['mod']]
		if self.vardic['other']['var']:
			list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
			list2 = list2 +  [x for x in self.vardic['other']['mod']]
		moddic = {}
		for x,y in zip(list1,list2):
			moddic[x] = y
			moddic[x] = y

		lookli = []
		for key in self.vardic.keys():
			lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
		lookdic = dict(lookli)
		_mreg = 'SUB\d{5,5}'
		mreg = RE.compile(_mreg)
		func2 = []
		for i,x in enumerate(func):
			func2.append(func[i])
			liobj = list(mreg.finditer(str(func2[i])))
			doneli = []
			for ma in reversed(liobj):
				suba = ma.group()
				if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
				   and suba not in doneli:
					func2[i] = func2[i].subs(sympycore.Symbol(suba),sympycore.exp(sympycore.Symbol(suba)))
					doneli.append(suba)
		self.func2 = func2

		# Also take into account log or not for other equations !
		if 'nlsubsys' in dir(self):
			tmpvarli0 = [x[0] for x in self.vardic['other']['var']]
			for i1,line in enumerate(self.func2[-lsubs:]):
				vari0 = tmpvarli0[i1]
				if 'log' in self.vardic['other']['mod'][i1]:
					pass
				else:
					vma = self.vreg(patup,vari0,False,'min')
					vmavar = vma[2][1]+'_bar'
					self.func2[i1] = sympycore.Symbol(vmavar)*(self.func2[i1])

		# Create list with log of steady states
		evalli = []
		alldic = {}
		alldic.update(self.paramdic)
		alldic.update(self.sstate)
		for x in intup:
			vma = self.vreg(patup, x, False, 'min')
			vari = vma[2][1]
			if 'log' in moddic[lookdic[vari]]:
				evalli.append(N.log(alldic[vari+'_bar']))
			else:
				evalli.append(alldic[vari+'_bar'])
		evaldic = {}
		for i,x in enumerate(tmpli):
			evaldic[x[1]] = evalli[i]

		# Now make the 2D and 3D symbolic and numeric Jacobian and Hessian
		def mkjaheseq(lcount,func2,jcols,symdic,tmpli,paramdic,sstate,evaldic,mk_hessian,exp,log):

			jdic = dict([[x,'0'] for x in range(jcols)])
			numj = numpy.matlib.zeros((1,jcols))
			if mk_hessian:
				rdic = dict([[x,'0'] for x in range(jcols)])
				hdic = dict([[x,copy.deepcopy(rdic)] for x in range(jcols)])
				numh = numpy.matlib.zeros((jcols,jcols))

			alldic = {}
			alldic.update(paramdic)
			alldic.update(sstate)
			alldic.update(evaldic)
			locals().update(alldic)
			globals().update(paramdic)
			globals().update(sstate)
			globals().update(evaldic)
			count = 0
			for y in range(jcols):
				jdic[y] = func2[lcount].diff(symdic[tmpli[y][1]])
				numj[0,y] = eval(str(jdic[y].evalf()))
				if mk_hessian:
					for z in range(jcols):
						hdic[y][z] = jdic[y].diff(symdic[tmpli[z][1]])
						numh[count,z] = eval(str(hdic[y][z].evalf()))
					count = count + 1
			if mk_hessian:
				return (numj,jdic,numh,hdic)
			else:
				return (numj,jdic)

		# Start parallel Python job server
		ppservers = ()
		inputs = [x for x in xrange(len(self.func2))]
		job_server = PP.Server(ncpus,ppservers=ppservers)
		from sympycore import exp
		from sympycore import log
		imports = ('numpy','copy','numpy.matlib',)
		jobs = [job_server.submit(mkjaheseq,(input,self.func2,jcols,symdic,tmpli,self.paramdic,self.sstate,evaldic,mk_hessian,exp,log),(),imports) for input in inputs]
		if mk_hessian:
			jdic = {}
			hdic = {}
			job_0 = jobs[0]
			numj = job_0()[0]
			jdic[0] = job_0()[1]
			numh = job_0()[2]
			hdic[0] = job_0()[3]
			for i1,job in enumerate(jobs[1:len(jobs)]):
				numj = MAT.vstack((numj,job()[0]))
				jdic[i1+1] = job()[1]
				numh = MAT.vstack((numh,job()[2]))
				hdic[i1+1] = job()[3]
			self.numj = numj
			self.jdic = jdic
			self.numh = numh
			self.hdic = hdic
		else:
			jdic = {}
			job_0 = jobs[0]
			numj = job_0()[0]
			jdic[0] = job_0()[1]
			for i1,job in enumerate(jobs[1:len(jobs)]):
				numj = MAT.vstack((numj,job()[0]))
				jdic[i1+1] = job()[1]
			self.numj = numj
			self.jdic = jdic

		if 'nlsubsys' in dir(self):
			numjs = numj[:-lsubs,:]
			numjl = numj[-lsubs:,:]
			self.numj = numjs
			self.numjl = numjl

			numhs = numh[:-lsubs*jcols,:]
			numhl = numh[-lsubs*jcols:,:]
			self.numhl = numhl
			self.numh = numhs
		else:
			self.numj = numj
			self.numh = numh

		self.jAA = self.numj[:,:int(len(intup)/2)]
		self.jBB = -self.numj[:,int(len(intup)/2):]
	# The numerical (Paul Klein) Jacobian and Hessian computation method (uses matlab)
	def mkjahenmat(self,msess=sess1):

		exo_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['exo']['var']]
		endo_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['endo']['var']]
		con_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['con']['var']]
		exo_0 = [x[0] for x in self.vardic['exo']['var']]
		endo_0 = [x[0].split('(')[0]+'(t-1)' for x in self.vardic['endo']['var']]
		con_0 = [x[0] for x in self.vardic['con']['var']]
		inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
		intup=tuple(inlist)

		patup = ('{-10,10}|None','endo|con|exo','{-10,10}')
		jcols = len(intup)
		jrows = len(self.nlsys_list)
		nlsys = COP.deepcopy(self.nlsys_list)

		if 'nlsubsys' in dir(self):
			lsubs = len(self.nlsubsys)
			jrows = jrows + lsubs
			nlsys = nlsys + COP.deepcopy(self.nlsubsys)

		tmpli = []
		for i,x in enumerate(intup):
			tmpli.append([x,'invar('+str(i+1)+')'])
			dicli = dict(tmpli)
			dicli2 = dict([[x[1],x[0]] for x in tmpli])

		func = []
		subli = []
		for x in nlsys:
			str_tmp = x[:]
			list1 = self.vreg(patup,x,True,'max')
			list1.reverse()
			for y in list1:
				pos = y[3][0]
				poe = y[3][1]
				vari = y[0]
				str_tmp = str_tmp[:pos] + dicli[vari] +str_tmp[poe:]
			list2 = self.vreg(('{-10,10}|None','iid','{-10,10}'),str_tmp,True,'max')
			if list2:
				list2.reverse()
				for y in list2:
					pos = y[3][0]
					poe = y[3][1]
					vari = y[0]
					str_tmp = str_tmp[:pos] + '0' +str_tmp[poe:]
			func.append(str_tmp)

			# Now substitute out exp and log in terms of matlab expressions
			elog = RE.compile('LOG\(')
			while elog.search(str_tmp):
				ma = elog.search(str_tmp)
				pos = ma.span()[0]
				poe = ma.span()[1]
				str_tmp = str_tmp[:pos]+'log('+str_tmp[poe:]
			eexp = RE.compile('EXP\(')
			while eexp.search(str_tmp):
				ma = eexp.search(str_tmp)
				pos = ma.span()[0]
				poe = ma.span()[1]
				str_tmp = str_tmp[:pos]+'exp('+str_tmp[poe:]
			func.append(eval(str_tmp))
		self.func1 = func


		# Only exp() when variable needs to be put into ln !
		list1 = [x[1] for x in self.vardic['endo']['var']]
		list1 = list1 + [x[1] for x in self.vardic['con']['var']]
		list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
		list2 = [x for x in self.vardic['endo']['mod']]
		list2 = list2 + [x for x in self.vardic['con']['mod']]
		list2 = list2 + [x for x in self.vardic['exo']['mod']]
		if self.vardic['other']['var']:
			list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
			list2 = list2 +  [x for x in self.vardic['other']['mod']]
		moddic = {}
		for x,y in zip(list1,list2):
			moddic[x] = y
			moddic[x] = y

		lookli = []
		for key in self.vardic.keys():
			lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
		lookdic = dict(lookli)
		_mreg = 'invar\(\d+\)'
		mreg = RE.compile(_mreg)
		func2 = []
		for i,x in enumerate(func):
			func2.append(func[i])
			liobj = list(mreg.finditer(str(func2[i])))
			doneli = []
			for ma in reversed(liobj):
				suba = ma.group()
				if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
				   and suba not in doneli:
					func2[i] = func2[i].subs(SPC.Symbol(suba),SPC.exp(SPC.Symbol(suba)))
					doneli.append(suba)
		self.func2 = func2

		# Also take into account log or not for other equations !
		if 'nlsubsys' in dir(self):
			tmpvarli0 = [x[0] for x in self.vardic['other']['var']]
			for i1,line in enumerate(self.func2[-lsubs:]):
				vari0 = tmpvarli0[i1]
				if 'log' in self.vardic['other']['mod'][i1]:
					pass
				else:
					vma = self.vreg(patup,vari0,False,'min')
					vmavar = vma[2][1]+'_bar'
					self.func2[i1] = SPC.Symbol(vmavar)*(self.func2[i1])

					vma = self.vreg(patup,vari0,False,'min')
					vmavar = vma[2][1]+'_bar'
					self.func2[i1] = SPC.Symbol(vmavar)*(self.func2[i1])


		# Transorm python exponentiation into matlab exponentiation
		_mreg='\*{2,2}'
		mreg=RE.compile(_mreg)
		for i,x in enumerate(self.func2):
			self.func2[i] = mreg.sub('^',self.func2[i])

		# Create matlab function
		if 'mfunc.m' in OPS.listdir(OPS.path.join(mlabpath,'Klein')):
			OPS.remove(OPS.path.join(mlabpath,'Klein/mfunc.m'))

		mfunc = open(OPS.path.join(mlabpath,'Klein/mfunc.m'),'w')
		mfunc.write('function vecreturn = mfunc(invar);\n\n')
		mfunc.write('%Parameters\n')
		for x in self.paramdic.items():
			mfunc.write(x[0]+' = '+str(x[1])+';\n')
		mfunc.write('\n\n')
		mfunc.write('%Steady States\n')
		for x in self.sstate.items():
			mfunc.write(x[0]+' = '+str(x[1])+';\n')
		mfunc.write('\n\n')
		mfunc.write('vecreturn = zeros('+str(len(self.func2))+',1);\n')
		mfunc.write('\n\n')
		for i,x in enumerate(self.func2):
			mfunc.write('vecreturn('+str(i+1)+') = '+x+';\n')
		mfunc.close()

		#Prepare ln inmatrix (vector)
		inmat = MAT.zeros((len(tmpli),1))
		alldic={}
		alldic.update(self.paramdic)
		alldic.update(self.sstate)
		for i,x in enumerate(tmpli):
			vma = self.vreg(patup, x, False, 'min')
			vari = vma[2][1]
			inmat[i,0] = alldic[vari+'_bar']
			if 'log' in moddic[lookdic[vari]]:
				inmat[i,0] = N.log(inmat[i,0])


		#Make Jacobian and Hessian
		sess1 = msess
		directory = OPS.path.join(mlabpath,'Klein')
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+directory)
		mlabraw.put(sess1,'x0',inmat)
		mlabraw.eval(sess1,"jacob = centgrad('mfunc',x0);")
		self.numj = mlabraw.get(sess1,'jacob')
		mlabraw.eval(sess1,"hessi = centhess('mfunc',x0);")
		self.numh = mlabraw.get(sess1,'hessi')

		numj = self.numj
		numh = self.numh

		if 'nlsubsys' in dir(self):
			numjs = numj[:-lsubs,:]
			numjl = numj[-lsubs:,:]
			self.numj = numjs
			self.numjl = numjl
			numhs = numh[:-lsubs*jcols,:]

			numhl = numh[-lsubs*jcols:,:]
			self.numhl = numhl
			self.numh = numhs
		else:
			self.numj = numj
			self.numh = numh

		self.jAA = N.matrix(self.numj[:,:int(len(intup)/2)])
		self.jBB = N.matrix(-self.numj[:,int(len(intup)/2):])

		del self.func1
		del self.func2

		OPS.remove(OPS.path.join(mlabpath,'Klein/mfunc.m'))

###NUMERICAL ONLY JACOBIAN AND HESSIAN METHODS - FOR UPDATING###############
	# The unparallelized Jacobian, Hessian computation method
	def mkjahen(self):
		'''
		This produces the Jacobian for the non-linear system of equations
		around the steady state, and evaluated at the steady state.
		Also produces the Hessian!
		'''
		exo_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['exo']['var']]
		endo_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['endo']['var']]
		con_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['con']['var']]
		exo_0 = [x[0] for x in self.vardic['exo']['var']]
		endo_0 = [x[0].split('(')[0]+'(t-1)' for x in self.vardic['endo']['var']]
		con_0 = [x[0] for x in self.vardic['con']['var']]
		inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
		intup=tuple(inlist)

		patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
		jcols = len(intup)
		jrows = len(self.nlsys_list)
		nlsys = COP.deepcopy(self.nlsys_list)

		if 'nlsubsys' in dir(self):
			lsubs = len(self.nlsubsys)
			jrows = jrows + lsubs
			nlsys = nlsys + COP.deepcopy(self.nlsubsys)

		# Create substitution var list and dictionary
		tmpli = []
		for i,x in enumerate(intup):
			nolen = len(str(i))
			inds = (5-nolen)*'0'+str(i)
			tmpli.append([x,'SUB'+inds])
			dicli = dict(tmpli)
			dicli2 = dict([[x[1],x[0]] for x in tmpli])
			
		lookli = []
		for key in self.vardic.keys():
			lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
		lookdic = dict(lookli)

		# Only exp() when variable needs to be put into ln !
		list1 = [x[1] for x in self.vardic['endo']['var']]
		list1 = list1 + [x[1] for x in self.vardic['con']['var']]
		list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
		list2 = [x for x in self.vardic['endo']['mod']]
		list2 = list2 + [x for x in self.vardic['con']['mod']]
		list2 = list2 + [x for x in self.vardic['exo']['mod']]
		if self.vardic['other']['var']:
			list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
			list2 = list2 +  [x for x in self.vardic['other']['mod']]
		moddic = {}
		for x,y in zip(list1,list2):
			moddic[x] = y
			moddic[x] = y

		# Create list with log of steady states
		evalli = []
		alldic = {}
		alldic.update(self.paramdic)
		alldic.update(self.sstate)
		for x in intup:
			vma = self.vreg(patup, x, False, 'min')
			vari = vma[2][1]
			if 'log' in moddic[lookdic[vari]]:
				evalli.append(N.log(alldic[vari+'_bar']))
			else:
				evalli.append(alldic[vari+'_bar'])
		evaldic = {}
		for i,x in enumerate(tmpli):
			evaldic[x[1]] = evalli[i]

		# Make 2D symbolic and numeric Jacobian
		def mkjac(jrows=jrows,jcols=jcols):
			jdic = self.jdic
			numj = MAT.zeros((jrows,jcols))
			alldic = {}
			alldic.update(self.paramdic)
			alldic.update(self.sstate)
			alldic.update(evaldic)
			locals().update(alldic)
			for x in range(jrows):
				for y in range(jcols):
					numj[x,y] = eval(str(jdic[x][y].evalf()))
			return numj

		# Now make 3D symbolic and numeric Hessian
		def mkhes(jrows=jrows,jcols=jcols):
			hdic = self.hdic
			hrows = jrows*jcols
			numh = MAT.zeros((hrows,jcols))
			alldic = {}
			alldic.update(self.paramdic)
			alldic.update(self.sstate)
			alldic.update(evaldic)
			locals().update(alldic)
			count = 0
			for x in range(jrows):
				for y in range(jcols):
					for z in range(jcols):
						numh[count,z] = eval(str(hdic[x][y][z].evalf()))
					count = count + 1
			return numh

		self.numj = mkjac()
		numj = self.numj
		if mk_hessian:
			self.numh = mkhes()
			numh = self.numh

		if 'nlsubsys' in dir(self):
			numjs = numj[:-lsubs,:]
			numjl = numj[-lsubs:,:]
			self.numj = numjs
			self.numjl = numjl
			if mk_hessian:
				numhs = numh[:-lsubs*jcols,:]
				numhl = numh[-lsubs*jcols:,:]
				self.numhl = numhl
				self.numh = numhs
		else:
			self.numj = numj
			if mk_hessian:
				self.numh = numh

		self.jAA = self.numj[:,:int(len(intup)/2)]
		self.jBB = -self.numj[:,int(len(intup)/2):]
	# The parallelized mkjahe version using parallel python
	def mkjaheppn(self):
		'''
		This produces the Jacobian for the non-linear system of equations
		around the steady state, and evaluated at the steady state.
		Also produces the Hessian!
		'''
		# import local sympycore
		import sympycore

		exo_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['exo']['var']]
		endo_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['endo']['var']]
		con_1 = ['E(t)|'+x[0].split('(')[0]+'(t+1)' for x in self.vardic['con']['var']]
		exo_0 = [x[0] for x in self.vardic['exo']['var']]
		endo_0 = [x[0].split('(')[0]+'(t-1)' for x in self.vardic['endo']['var']]
		con_0 = [x[0] for x in self.vardic['con']['var']]
		inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
		intup=tuple(inlist)

		patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
		jcols = len(intup)
		jrows = len(self.nlsys_list)
		nlsys = COP.deepcopy(self.nlsys_list)

		if 'nlsubsys' in dir(self):
			lsubs = len(self.nlsubsys)
			jrows = jrows + lsubs
			nlsys = nlsys + COP.deepcopy(self.nlsubsys)

		# Create substitution var list and dictionary
		tmpli = []
		for i,x in enumerate(intup):
			nolen = len(str(i))
			inds = (5-nolen)*'0'+str(i)
			tmpli.append([x,'SUB'+inds])
			dicli = dict(tmpli)
			dicli2 = dict([[x[1],x[0]] for x in tmpli])

		lookli = []
		for key in self.vardic.keys():
			lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
		lookdic = dict(lookli)

		# Only exp() when variable needs to be put into ln !
		list1 = [x[1] for x in self.vardic['endo']['var']]
		list1 = list1 + [x[1] for x in self.vardic['con']['var']]
		list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
		list2 = [x for x in self.vardic['endo']['mod']]
		list2 = list2 + [x for x in self.vardic['con']['mod']]
		list2 = list2 + [x for x in self.vardic['exo']['mod']]
		if self.vardic['other']['var']:
			list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
			list2 = list2 +  [x for x in self.vardic['other']['mod']]
		moddic = {}
		for x,y in zip(list1,list2):
			moddic[x] = y
			moddic[x] = y

		# Create list with log of steady states
		evalli = []
		alldic = {}
		alldic.update(self.paramdic)
		alldic.update(self.sstate)
		for x in intup:
			vma = self.vreg(patup, x, False, 'min')
			vari = vma[2][1]
			if 'log' in moddic[lookdic[vari]]:
				evalli.append(N.log(alldic[vari+'_bar']))
			else:
				evalli.append(alldic[vari+'_bar'])
		evaldic = {}
		for i,x in enumerate(tmpli):
			evaldic[x[1]] = evalli[i]

		# Now make the 2D and 3D symbolic and numeric Jacobian and Hessian
		def mkjaheseq(lcount,func2,jcols,tmpli,paramdic,sstate,evaldic,mk_hessian,exp,log,jdic,hdic):
			numj = numpy.matlib.zeros((1,jcols))
			if mk_hessian:
				numh = numpy.matlib.zeros((jcols,jcols))

			alldic = {}
			alldic.update(paramdic)
			alldic.update(sstate)
			alldic.update(evaldic)
			locals().update(alldic)
			globals().update(paramdic)
			globals().update(sstate)
			globals().update(evaldic)
			count = 0
			for y in range(jcols):
				numj[0,y] = eval(str(jdic[lcount][y].evalf()))
				if mk_hessian:
					for z in range(jcols):
						numh[count,z] = eval(str(hdic[lcount][y][z].evalf()))
					count = count + 1
			if mk_hessian:
				return (numj,numh)
			else:
				return numj

		# Start parallel Python job server
		ppservers = ()
		inputs = [x for x in xrange(len(self.func2))]
		job_server = PP.Server(ncpus,ppservers=ppservers)
		from sympycore import exp
		from sympycore import log
		imports = ('numpy','copy','numpy.matlib',)
		jobs = [job_server.submit(mkjaheseq,(input,self.func2,jcols,tmpli,self.paramdic,self.sstate,evaldic,mk_hessian,exp,log,self.jdic,self.hdic),(),imports) for input in inputs]
		if mk_hessian:
			job_0 = jobs[0]
			numj = job_0()[0]
			numh = job_0()[1]
			for i1,job in enumerate(jobs[1:len(jobs)]):
				numj = MAT.vstack((numj,job()[0]))
				numh = MAT.vstack((numh,job()[1]))
			self.numj = numj
			self.numh = numh
		else:
			job_0 = jobs[0]
			numj = job_0()
			for i1,job in enumerate(jobs[1:len(jobs)]):
				numj = MAT.vstack((numj,job()))
			self.numj = numj

		if 'nlsubsys' in dir(self):
			numjs = numj[:-lsubs,:]
			numjl = numj[-lsubs:,:]
			self.numj = numjs
			self.numjl = numjl

			numhs = numh[:-lsubs*jcols,:]
			numhl = numh[-lsubs*jcols:,:]
			self.numhl = numhl
			self.numh = numhs
		else:
			self.numj = numj
			self.numh = numh

		self.jAA = self.numj[:,:int(len(intup)/2)]
		self.jBB = -self.numj[:,int(len(intup)/2):]
##########################################################	
	# Method updating model IF model file has been changed externally !
	def updf(self):
		self.txtpars.__init__(self.modfile)

	# Method updating model after anything has been changed manually, like paramvalue, etc !
	def updm(self):
		# Create None tester regular expression
		_nreg = '^\s*None\s*$'
		nreg = RE.compile(_nreg)
################## STEADY STATE CALCULATIONS !!! ####################
		self.sssolvers = SSsolvers()
		# Solve for steady-state using fsolve
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['ssm'][0]]) == 0:
			intup = (self.ssys_list,self.ssidic,self.paramdic)
			self.sssolvers.fsolve = Fsolve(intup)
			self.sssolvers.fsolve.solve()
			if self.sssolvers.fsolve.ier == 1:
				self.sstate = self.sssolvers.fsolve.fsout
				self.ssidic = self.sssolvers.fsolve.fsout
				self.numssdic = self.sssolvers.fsolve.fsout
				self.switches['ss_suc'] = ['1','1']
			else:
				self.switches['ss_suc'] = ['1','0']
		# Solve for steady-state using manss
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['sss'][0]]) == 0:
			intup = (self.manss_sys,self.paramdic)
			self.sssolvers.manss = Manss(intup)
			self.sssolvers.manss.solve()
			self.sstate = self.sssolvers.manss.sstate
		if initlev == 1: return

		# No populate more with stuff that needs steady state
		self.pop2(self.txtpars)

		# Open the model solution tree branch
		self.modsolvers = MODsolvers()
		######################## LINEAR METHODS !!! ############################
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['modeq'][0]]) == 0:
			# Open the matlab Uhlig object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1,
				 self.vardic)
			self.modsolvers.matuhlig = MatUhlig(intup)
			# Open the native Uhlig object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1)
			self.modsolvers.pyuhlig = PyUhlig(intup)
			# Open the matlab Klein object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1)
			self.modsolvers.matklein = MatKlein(intup)
			# Open the Fortran Klein object
			intup = ((self.nendo,self.ncon,self.nexo),
				 self.eqindx,
				 self.vreg,
				 self.llsys_list,
				 self.diffli1,
				 self.diffli2,
				 sess1)
			self.modsolvers.forklein = ForKlein(intup)
		################## 1ST-ORDER NON-LINEAR METHODS !!! ##################
		if sum([nreg.search(x)!=None for x in self.txtpars.secs['focs'][0]]) == 0:

			# First, create the Jacobian and (possibly-->mk_hessian==True?) Hessian
			if use_anaderiv:
				if ncpus > 1 and mk_hessian:
					self.mkjahepp()
				elif ncpus > 1 and not mk_hessian:
					self.mkjahepp()
				else:
					self.mkjahe()
			else:
				self.mkjahenmat()

			# Open the MatWood object
			intup = (self.jAA,self.jBB,
				 self.nexo,self.ncon,
				 self.nendo,sess1)
			self.modsolvers.matwood = MatWood(intup)
			# Open the Fortran KleinD object
			if 'nlsubsys' in dir(self):
				intup = (self.numj,
					 self.nendo,self.nexo,
					 self.ncon,self.sigma,
					 self.jAA,self.jBB,
					 self.vardic,self.vdic,
					 self.modname,
					 self.numjl,
					 self.nother)
			else:
				intup = (self.numj,
					 self.nendo,self.nexo,
					 self.ncon,self.sigma,
					 self.jAA,self.jBB,
					 self.vardic,self.vdic,
					 self.modname)
			self.modsolvers.forkleind = ForKleinD(intup)
		################## 2ND-ORDER NON-LINEAR METHODS !!! ##################
			if sum([nreg.search(x)!=None for x in self.txtpars.secs['vcvm'][0]]) == 0 and\
			   'numh' in dir(self):
				# Open the MatKlein2D object
				if 'nlsubsys' in dir(self):
					intup = (self.numj,self.numh,
						 self.nendo,self.nexo,
						 self.ncon,self.sigma,
						 self.jAA,self.jBB,
						 self.vardic,self.vdic,
						 self.modname,
						 self.numjl,self.numhl,
						 self.nother,sess1)
				else:
					intup = (self.numj,self.numh,
						 self.nendo,self.nexo,
						 self.ncon,self.sigma,
						 self.jAA,self.jBB,
						 self.vardic,self.vdic,
						 self.modname,
						 sess1)
				self.modsolvers.matklein2d = MatKlein2D(intup)
				# Open the PyKlein2D object
				intup = intup[:-1]
				self.modsolvers.pyklein2d = PyKlein2D(intup)

		if 'jAA' in dir(self):
			self.mkeigv()

	def mkeigv(self):
		AA = self.jAA
		BB = self.jBB
		try:
			CC = BB.I*AA
			self.eigvals = LIN.eig(CC)[0]
		except:
			pass
"""***********************************************************"""
################THE GENERAL TEXTPARSER (WORKS)#######################
class TXTparser:

	def __init__(self,ffile=None,type='type1'):
		self.type = type
		self.ffile = ffile
		self.filename = ffile
		self.secs = {}
		input = open(OPS.path.join(modfpath,ffile), 'r')
		wholefile = input.read()
		self.filestring = wholefile
		self.lines = self.filestring.splitlines()
		self.numlines = len(self.filestring.splitlines())

		secnames = [['%Model Description','mod','MOD_loc'],
			    ['%Model Information','info','INF_loc'],
			    ['%Parameters','para','PA_loc'],
			    ['%Variable Vectors','varvec','VV_loc'],
			    ['%Boundary conditions','bocond','BC_loc'],
			    ['%Variable Substitution Non-Linear System','vsfocs','VSFO_loc'],
			    ['%Non-linear first-order conditions','focs','FO_loc'],
			    ['%Steady States[closed-form]','sss','SS_loc'],
			    ['%Manual entry of sstate non-linear system','ssm','SSM_loc'],
			    ['%Log-Linearized Model Equations','modeq','ME_loc'],
			    ['%Variance-Covariance Matrix','vcvm','VCM_loc'],
			    ['%Minford Model Evaluation','mme','MME_loc']]

		if type == 'type1':
			self.type1(secnames)
		else:
			pass

	def type1(self,secnames):
		lines = self.lines
		numlines = self.numlines
		secs = self.secs		
		locdic = locate(lines,secnames)

		def readerfu(locname,locat,secname):
			list_tmp1 = []
			list_tmp2 = []
			row_iter=locat+1
			for x in lines[row_iter:]:
				x = RE.sub("\s+", "", x)
				if '#' not in x and x != '' and x[0] != '%':
					list_tmp1.append(lines[row_iter])
					list_tmp2.append(lines[row_iter])
					row_iter=row_iter+1
				elif '#' in x:
					list_tmp2.append(lines[row_iter])
					row_iter=row_iter+1
				elif x == '':
					row_iter=row_iter+1
				elif x[0] == '%':
					break
			while '' in list_tmp1:
				list_tmp1.remove('')
			self.secs[secname] = (list_tmp1,list_tmp2)

		for x in secnames:
			readerfu(x[1],locdic[x[1]],x[1])
"""***********************************************************"""
###########THE MODSOLVER CLASS AND IT'S SUBCLASSES (WORKS)###############
class MODsolvers:

	def __init__(self):
		pass
#----------------------------------------------------------------------------------------------------------------------
class PyUhlig:

	def __init__(self,intup):
		self.ntup = intup[0]
		self.nendo = self.ntup[0]
		self.ncon = self.ntup[1]
		self.nexo = self.ntup[2]

		self.eqindx = intup[1]
		self.vreg = intup[2]
		self.llsys_list = intup[3]
		self.diffli1 = intup[4]
		self.diffli2 = intup[5]
		diffli1 = self.diffli1
		diffli2 = self.diffli2

		deteq = []
		for x in self.eqindx['det']:
			deteq.append(self.llsys_list[x])
		self.deteq = deteq
		expeq = []
		for x in self.eqindx['exp']:
			expeq.append(self.llsys_list[x])
		self.expeq = expeq
		erreq = []
		for x in self.eqindx['err']:
			erreq.append(self.llsys_list[x])
		self.erreq = erreq

		detdif1 = []
		detdif2 = []
		for x in self.eqindx['det']:
			detdif1.append(diffli1[x])
			detdif2.append(diffli2[x])
		expdif1 = []
		expdif2 = []
		for x in self.eqindx['exp']:
			expdif1.append(diffli1[x])
			expdif2.append(diffli2[x])
		errdif1 = []
		errdif2 = []
		for x in self.eqindx['err']:
			errdif1.append(diffli1[x])
			errdif2.append(diffli2[x])
		self.detdif1 = detdif1
		self.detdif2 = detdif2
		self.expdif1 = expdif1
		self.expdif2 = expdif2
		self.errdif1 = errdif1
		self.errdif2 = errdif2

		self.mkmats()

	def mkmats(self):
		deteq = self.deteq
		expeq = self.expeq
		erreq = self.erreq
		nendo = self.nendo
		ncon = self.ncon
		nexo = self.nexo		
		vreg = self.vreg
		detdif1 = self.detdif1
		detdif2 = self.detdif2
		expdif1 = self.expdif1
		expdif2 = self.expdif2
		errdif1 = self.errdif1
		errdif2 = self.errdif2

		# Create argument list for matrix creation
		argli = [['AA',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','0')],
			 ['BB',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','-1')],
			 ['CC',(detdif2,detdif1),deteq,(len(deteq),ncon),(None,'con','0')],
			 ['DD',(detdif2,detdif1),deteq,(len(deteq),nexo),(None,'exo','0')],
			 ['FF',(expdif2,expdif1),expeq,(len(expeq),nendo),('0','endo','1')],
			 ['GG',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','0')],
			 ['HH',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','-1')],
			 ['JJ',(expdif2,expdif1),expeq,(len(expeq),ncon),('0','con','1')],
			 ['KK',(expdif2,expdif1),expeq,(len(expeq),ncon),(None,'con','0')],
			 ['LL',(expdif2,expdif1),expeq,(len(expeq),nexo),('0','exo','1')],
			 ['MM',(expdif2,expdif1),expeq,(len(expeq),nexo),(None,'exo','0')],
			 ['NN',(errdif2,errdif1),erreq,(len(erreq),nexo),(None,'exo','0')]]

		# Create all matrices in a loop over argli
		for argx in argli:
			XX = MAT.zeros(argx[3])
			dic1 = dict([[x,'0'] for x in range(argx[3][1])])
			sXX = dict([[x,COP.deepcopy(dic1)] for x in range(len(argx[2]))])
			for y,x in enumerate(argx[2]):
				if vreg(argx[4],x,False,'max'):
					cont=[[z[0],z[1][1]] for z in vreg(argx[4],x,True,'min')]
					for z in cont:
						XX[y,z[1]] = argx[1][0][y][z[0]]
						sXX[y][z[1]] = argx[1][1][y][z[0]]
			exec('self.'+argx[0]+' = XX')
			exec('self.'+'s'+argx[0]+' = sXX')

	def solve(self,sol_method='do_QZ',Xi_select='all'):
		self.output = {}
		self.Xi_select = Xi_select
		self.sol_method = sol_method
		if self.sol_method == 'do_PP':
			self.do_PP(self.Xi_select)
		elif self.sol_method == 'do_QZ':
			self.do_QZ(self.Xi_select,Tol=1.0000e-06)
		self.do_QRS()

	def do_QZ(self,Xi_select='all',Tol=1.0000e-06):
		# Make uhlig matrices locally available for computations
		AA = self.AA
		BB = self.BB
		CC = self.CC
		DD = self.DD
		FF = self.FF
		GG = self.GG
		HH = self.HH
		JJ = self.JJ
		KK = self.KK
		LL = self.LL
		MM = self.MM
		NN = self.NN
		q_expectational_equ = N.shape(FF)[0]
		m_states = N.shape(FF)[1]
		l_equ = N.shape(CC)[0]
		n_endog = N.shape(CC)[1]
		k_exog = min(N.shape(NN))

		if HLP.rank(CC) < n_endog:
			raise uhlerr, 'uhlerror: Rank(CC) needs to be at least\n\
			      equal to the number of endogenous variables'
		if l_equ == 0:
			CC_plus = LIN.pinv(CC)
			CC_0 = (HLP.null(CC.T)).T
			Psi_mat = FF
			Gamma_mat = -GG
			Theta_mat = -HH
			Xi_mat = \
			       N.concatenate((Gamma_mat,Theta_mat,MAT.eye(m_states),\
					      MAT.zeros((m_states,m_states))),1)
			Delta_mat = \
				  N.concatenate((Psi_mat,MAT.zeros((m_states,m_states)),\
						 MAT.zeros((m_states,m_states)),\
						 MAT.zeros((m_states,m_states)),MAT.eye(m_states),1))


		else:
			CC_plus = LIN.pinv(CC)
			CC_0 = HLP.null(CC.T)
			if l_equ > n_endog:
				Psi_mat = \
					N.concatenate((MAT.zeros((l_equ-n_endog,m_states)),FF-JJ*CC_plus*AA),1)
			elif l_equ == n_endog:
				Psi_mat = N.concatenate((FF-JJ*CC_plus*AA),1)
			Gamma_mat = N.concatenate((CC_0*AA, JJ*CC_plus*BB-GG+KK*CC_plus*AA))
			Theta_mat = N.concatenate((CC_0*BB, KK*CC_plus*BB-HH))
			Xi_mat = \
			       N.concatenate((\
				       N.concatenate((Gamma_mat,Theta_mat),1),\
				       N.concatenate((MAT.identity(m_states),MAT.zeros((m_states,m_states))),1)\
			       ))
			Psi_mat = Psi_mat.reshape(m_states,m_states)
			Delta_mat = \
				  N.concatenate((\
					  N.concatenate((Psi_mat,MAT.zeros((m_states,m_states))),1),\
					  N.concatenate((MAT.zeros((m_states,m_states)),MAT.identity(m_states)),1)\
				  ))

		(Delta_up,Xi_up,UUU,VVV)=\
		 HLP.qz(Delta_mat,Xi_mat,\
			mode='complex',\
			lapackname=lapackname,\
			lapackpath=lapackpath)

		d_Delta_up = MAT.diag(Delta_up)
		d_Xi_up = MAT.diag(Xi_up)
		Xi_eigval = MAT.zeros(N.shape(Delta_up))
		for i1 in range(0,N.shape(Delta_up)[0],1):
			Xi_eigval[i1,i1] = d_Xi_up[i1]/d_Delta_up[i1]
		d_Xi_eigval = N.diag(Xi_eigval)
		mat_tmp = MAT.zeros((N.shape(Xi_eigval)[0],3))
		i1=0
		for x in d_Xi_eigval:
			mat_tmp[i1,0] = d_Xi_eigval[i1]
			mat_tmp[i1,1] = abs(d_Xi_eigval)[i1]
			mat_tmp[i1,2] = i1
			i1=i1+1

		Xi_sortval = HLP.sortrows(mat_tmp,1)[:,0]
		# Need to do an argsort() on index column to turn to integers (Booleans) !!
		Xi_sortindex = HLP.sortrows(mat_tmp,1)[:,2].argsort(axis=0)
		Xi_sortabs = HLP.sortrows(mat_tmp,1)[:,1]

		# Root selection branch with Xi_select
		if Xi_select == 'all':
			Xi_select = N.arange(0,m_states,1)

		stake = max(abs(Xi_sortval[Xi_select])) + Tol
		(Delta_up,Xi_up,UUU,VVV) = self.qzdiv(stake,Delta_up,Xi_up,UUU,VVV)
		Lamda_mat = N.diag(Xi_sortval[Xi_select])
		trVVV = VVV.conjugate().T
		VVV_2_1 = trVVV[m_states:2*m_states,0:m_states]
		VVV_2_2 = trVVV[m_states:2*m_states,m_states:2*m_states]
		UUU_2_1 = UUU[m_states:2*m_states,0:m_states]

		PP = -(VVV_2_1.I*VVV_2_2)
		PP = N.real(PP)

		self.PP = PP
		self.CC_plus = CC_plus

	def do_PP(self,Xi_select='all'):
		# Make uhlig matrices locally available for computations
		AA = self.AA
		BB = self.BB
		CC = self.CC
		DD = self.DD
		FF = self.FF
		GG = self.GG
		HH = self.HH
		JJ = self.JJ
		KK = self.KK
		LL = self.LL
		MM = self.MM
		NN = self.NN
		q_expectational_equ = N.shape(FF)[0]
		m_states = N.shape(FF)[1]
		l_equ = N.shape(CC)[0]
		n_endog = N.shape(CC)[1]
		k_exog = min(N.shape(NN))

		if HLP.rank(CC) < n_endog:
			raise uhlerr, 'uhlerror: Rank(CC) needs to be at least\n\
			      equal to the number of endogenous variables'
		if l_equ == 0:
			CC_plus = LIN.pinv(CC)
			CC_0 = (HLP.null(CC.T)).T
			Psi_mat = FF
			Gamma_mat = -GG
			Theta_mat = -HH
			Xi_mat = \
			       N.concatenate((Gamma_mat,Theta_mat,MAT.eye(m_states),\
					      MAT.zeros((m_states,m_states))),1)
			Delta_mat = \
				  N.concatenate((Psi_mat,MAT.zeros((m_states,m_states)),\
						 MAT.zeros((m_states,m_states)),\
						 MAT.zeros((m_states,m_states)),MAT.eye(m_states),1))


		else:
			CC_plus = LIN.pinv(CC)
			CC_0 = HLP.null(CC.T)
			if l_equ > n_endog:
				Psi_mat = \
					N.concatenate((MAT.zeros((l_equ-n_endog,m_states)),FF-JJ*CC_plus*AA),1)
			elif l_equ == n_endog:
				Psi_mat = N.concatenate((FF-JJ*CC_plus*AA),1)
			Gamma_mat = N.concatenate((CC_0*AA, JJ*CC_plus*BB-GG+KK*CC_plus*AA))
			Theta_mat = N.concatenate((CC_0*BB, KK*CC_plus*BB-HH))
			Xi_mat = \
			       N.concatenate((\
				       N.concatenate((Gamma_mat,Theta_mat),1),\
				       N.concatenate((MAT.eye(m_states),MAT.zeros((m_states,m_states))),1)\
			       ))
			Delta_mat = \
				  N.concatenate((\
					  N.concatenate((Psi_mat,MAT.zeros((m_states,m_states))),1),\
					  N.concatenate((MAT.zeros((m_states,m_states)),MAT.eye(m_states)),1)\
				  ))

		(Xi_eigvec,Xi_eigval) = HLP.eig(Xi_mat,Delta_mat)
		tmp_mat = MAT.eye(N.rank(Xi_eigvec))
		for i1 in range(0,N.rank(Xi_eigvec),1):
			tmp_mat[i1,i1] = float(Xi_eigval[0,i1])
		Xi_eigval = tmp_mat


		if HLP.rank(Xi_eigvec) < m_states:
			raise uhlerr, 'uhlerror: Xi_mat is not diagonalizable!\n\
			      Cannot solve for PP. Maybe you should try the Schur-Decomposition\n\
			      method instead, use do_QZ()!!'


		d_Xi_eigval = N.diag(Xi_eigval)
		mat_tmp = MAT.zeros((N.rank(Xi_eigval),3))
		i1=0
		for x in d_Xi_eigval:
			mat_tmp[i1,0] = d_Xi_eigval[i1]
			mat_tmp[i1,1] = abs(d_Xi_eigval[i1])
			mat_tmp[i1,2] = i1
			i1=i1+1
		Xi_sortval = HLP.sortrows(mat_tmp,1)[:,0]
		# Need to do an argsort() on index column to turn to integers (Booleans) !!
		Xi_sortindex = HLP.sortrows(mat_tmp,1)[:,2].argsort(axis=0)
		Xi_sortabs = HLP.sortrows(mat_tmp,1)[:,1]
		Xi_sortvec = Xi_eigvec[0:2*m_states,:].take(Xi_sortindex.T.A, axis=1)

		# Root selection branch with Xi_select
		if Xi_select == 'all':
			Xi_select = N.arange(0,m_states,1)

		#Create Lambda_mat and Omega_mat
		Lambda_mat = MAT.zeros((len(Xi_select),len(Xi_select)))
		for i1 in range(0,len(Xi_select),1):
			Lambda_mat[i1,i1] = Xi_sortval[Xi_select][i1]
		Omega_mat = Xi_sortvec[m_states:2*m_states,:].take(Xi_select, axis=1)

		if HLP.rank(Omega_mat) < m_states:
			raise uhlerr, 'uhlerror: Omega_mat is not invertible!!\n\
			      Therefore, solution for PP is not available'

		PP = Omega_mat*HLP.ctr((HLP.ctr(Omega_mat).I*HLP.ctr(Lambda_mat)))
		self.PP = PP
		self.CC_plus = CC_plus

	def do_QRS(self):
		# Make uhlig matrices locally available for computations
		AA = self.AA
		BB = self.BB
		CC = self.CC
		DD = self.DD
		FF = self.FF
		GG = self.GG
		HH = self.HH
		JJ = self.JJ
		KK = self.KK
		LL = self.LL
		MM = self.MM
		NN = self.NN
		PP = self.PP
		CC_plus = self.CC_plus
		(l_equ,m_states) = MAT.shape(AA)
		(l_equ,n_endog) = MAT.shape(CC)
		(l_equ,k_exog) = MAT.shape(DD)
		if l_equ == 0:
			RR = MAT.zeros((0,m_states))
			VV1 = MAT.kron(MAT.transpose(NN),FF)+MAT.kron(MAT.identity(k_exog),FF*PP+GG)
			VV2 = MAT.kron(MAT.transpose(NN),JJ)+MAT.kron(MAT.identity(k_exog),KK)
			VV = MAT.hstack((VV1,VV2))
		else:
			RR = -CC_plus*(AA*PP+BB)
			VV1 = MAT.kron(MAT.identity(k_exog),AA)
			VV2 = MAT.kron(MAT.identity(k_exog),CC)
			VV3 = MAT.kron(MAT.transpose(NN),FF)+MAT.kron(MAT.identity(k_exog),FF*PP+JJ*RR+GG)
			VV4 = MAT.kron(MAT.transpose(NN),JJ)+MAT.kron(MAT.identity(k_exog),KK)
			VV = MAT.vstack((MAT.hstack((VV1,VV2)),MAT.hstack((VV3,VV4))))
		self.RR = RR

		LLNN_plus_MM = LL*NN+MM
		NON = MAT.hstack((DD.flatten(1),LLNN_plus_MM.flatten(1)))
		try:
			QQSS_vec = -(VV.I*MAT.transpose(NON))
		except MyErr:
			print 'Error: Matrix VV is not invertible!'
		QQ = QQSS_vec[0:m_states*k_exog,:].reshape(-1,m_states).transpose()
		SS = QQSS_vec[(m_states*k_exog):((m_states+n_endog)*k_exog),:].reshape(-1,n_endog).transpose()
		WW1 = MAT.hstack((MAT.identity(m_states),MAT.zeros((m_states,k_exog))))
		WW2 = MAT.hstack((RR*LIN.pinv(PP),SS-RR*LIN.pinv(PP)*QQ))
		WW3 = MAT.hstack((MAT.zeros((k_exog,m_states)),MAT.identity(k_exog)))
		WW = MAT.vstack((WW1,WW2,WW3))
		self.WW = WW
		self.QQ = QQ
		self.SS = SS
		del self.CC_plus

	def qzdiv(self,stake,A,B,Q,Z):
		n = N.shape(A)[0]
		root = N.hstack((abs(N.mat(N.diag(A)).T),abs(N.mat(N.diag(B)).T)))
		index_mat = (root[:,0]<1e-13).choose(root[:,0],1)
		index_mat = (index_mat>1e-13).choose(root[:,0],0)
		root[:,0] = root[:,0]-MAT.multiply(index_mat,(root[:,0]+root[:,1]))
		root[:,1] = root[:,1]/root[:,0]
		for i1 in range(n-1,-1,-1):
			m='none'
			for i2 in range(i1,-1,-1):
				if root[i2,1] > stake or root[i2,1] < -0.1:
					m=i2
					break
			if m == 'none':
				return (A,B,Q,Z)
			else:
				for i3 in range(m,i1,1):
					(A,B,Q,Z) = self.qzswitch(i3,A,B,Q,Z)
					tmp = COP.deepcopy(root[i3,1])
					root[i3,1] = root[i3+1,1]
					root[i3+1,1] = tmp

	def qzswitch(self,i,A,B,Q,Z):
		a = A[i,i]
		d = B[i,i]
		b = A[i,i+1]
		e = B[i,i+1]
		c = A[i+1,i+1]
		f = B[i+1,i+1]
		wz = N.mat(N.hstack(((c*e-f*b),(c*d-f*a).conjugate().T)))
		xy = N.mat(N.hstack(((b*d-e*a).conjugate().T,(c*d-f*a).conjugate().T)))
		n = N.mat(N.sqrt(wz*wz.conjugate().T))
		m = N.mat(N.sqrt(xy*xy.conjugate().T))
		if n.all() == 0:
			return
		else:
			wz = N.mat(wz/n)
			xy = N.mat(xy/m)
			wz = N.vstack((wz,N.hstack((-wz[:,1].T,wz[:,0].T))))
			xy = N.vstack((xy,N.hstack((-xy[:,1].T,xy[:,0].T))))
			A[i:i+2,:] = xy*A[i:i+2,:]
			B[i:i+2,:] = xy*B[i:i+2,:]
			A[:,i:i+2] = A[:,i:i+2]*wz
			B[:,i:i+2] = B[:,i:i+2]*wz
			Z[:,i:i+2] = Z[:,i:i+2]*wz
			Q[i:i+2,:] = xy*Q[i:i+2,:]
		return (A,B,Q,Z)
#----------------------------------------------------------------------------------------------------------------------
class MatUhlig:

	def __init__(self,intup):
		self.ntup = intup[0]
		self.nendo = self.ntup[0]
		self.ncon = self.ntup[1]
		self.nexo = self.ntup[2]

		self.eqindx = intup[1]
		self.vreg = intup[2]
		self.llsys_list = intup[3]
		self.diffli1 = intup[4]
		self.diffli2 = intup[5]
		self.sess1 = intup[6]
		self.vardic = intup[7]
		diffli1 = self.diffli1
		diffli2 = self.diffli2

		deteq = []
		for x in self.eqindx['det']:
			deteq.append(self.llsys_list[x])
		self.deteq = deteq
		expeq = []
		for x in self.eqindx['exp']:
			expeq.append(self.llsys_list[x])
		self.expeq = expeq
		erreq = []
		for x in self.eqindx['err']:
			erreq.append(self.llsys_list[x])
		self.erreq = erreq

		detdif1 = []
		detdif2 = []
		for x in self.eqindx['det']:
			detdif1.append(diffli1[x])
			detdif2.append(diffli2[x])
		expdif1 = []
		expdif2 = []
		for x in self.eqindx['exp']:
			expdif1.append(diffli1[x])
			expdif2.append(diffli2[x])
		errdif1 = []
		errdif2 = []
		for x in self.eqindx['err']:
			errdif1.append(diffli1[x])
			errdif2.append(diffli2[x])
		self.detdif1 = detdif1
		self.detdif2 = detdif2
		self.expdif1 = expdif1
		self.expdif2 = expdif2
		self.errdif1 = errdif1
		self.errdif2 = errdif2

		self.mkmats()

	def mkmats(self):
		deteq = self.deteq
		expeq = self.expeq
		erreq = self.erreq
		nendo = self.nendo
		ncon = self.ncon
		nexo = self.nexo		
		vreg = self.vreg
		detdif1 = self.detdif1
		detdif2 = self.detdif2
		expdif1 = self.expdif1
		expdif2 = self.expdif2
		errdif1 = self.errdif1
		errdif2 = self.errdif2

		# Create argument list for matrix creation
		argli = [['AA',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','0')],
			 ['BB',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','-1')],
			 ['CC',(detdif2,detdif1),deteq,(len(deteq),ncon),(None,'con','0')],
			 ['DD',(detdif2,detdif1),deteq,(len(deteq),nexo),(None,'exo','0')],
			 ['FF',(expdif2,expdif1),expeq,(len(expeq),nendo),('0','endo','1')],
			 ['GG',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','0')],
			 ['HH',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','-1')],
			 ['JJ',(expdif2,expdif1),expeq,(len(expeq),ncon),('0','con','1')],
			 ['KK',(expdif2,expdif1),expeq,(len(expeq),ncon),(None,'con','0')],
			 ['LL',(expdif2,expdif1),expeq,(len(expeq),nexo),('0','exo','1')],
			 ['MM',(expdif2,expdif1),expeq,(len(expeq),nexo),(None,'exo','0')],
			 ['NN',(errdif2,errdif1),erreq,(len(erreq),nexo),(None,'exo','0')]]

		# Create all matrices in a loop over argli
		for argx in argli:
			XX = MAT.zeros(argx[3])
			dic1 = dict([[x,'0'] for x in range(argx[3][1])])
			sXX = dict([[x,COP.deepcopy(dic1)] for x in range(len(argx[2]))])
			for y,x in enumerate(argx[2]):
				if vreg(argx[4],x,False,'max'):
					cont=[[z[0],z[1][1]] for z in vreg(argx[4],x,True,'min')]
					for z in cont:
						XX[y,z[1]] = argx[1][0][y][z[0]]
						sXX[y][z[1]] = argx[1][1][y][z[0]]
			exec('self.'+argx[0]+' = XX')
			exec('self.'+'s'+argx[0]+' = sXX')

	def solve(self,sol_method='do_QZ',Xi_select='all'):
		# Create varnames with correct string length
		tmp_list =[x[1] for x in self.vardic['endo']['var']]\
			 +[x[1] for x in self.vardic['con']['var']]\
			 +[x[1] for x in self.vardic['exo']['var']]
		tmp_list2 = [[len(x),x] for x in tmp_list]
		tmp_list2.sort()
		max_len = tmp_list2[-1][0]
		i1=0
		for x in tmp_list:
			tmp_list[i1]=x+(max_len-len(x))*' '
			i1=i1+1
		varnstring = 'VARNAMES = ['
		for x in tmp_list:
			varnstring = varnstring + "'" +x[1]+ "'" + ','
		varnstring = varnstring[0:-1]+'];'

		# Start matlab session and calculate
		sess1 = self.sess1
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,varnstring)
		mlabraw.put(sess1,'AA',self.AA)
		mlabraw.put(sess1,'BB',self.BB)
		mlabraw.put(sess1,'CC',self.CC)
		mlabraw.put(sess1,'DD',self.DD)
		mlabraw.put(sess1,'FF',self.FF)
		mlabraw.put(sess1,'GG',self.GG)
		mlabraw.put(sess1,'HH',self.HH)
		mlabraw.put(sess1,'JJ',self.JJ)
		mlabraw.put(sess1,'KK',self.KK)
		mlabraw.put(sess1,'LL',self.LL)
		mlabraw.put(sess1,'MM',self.MM)
		mlabraw.put(sess1,'NN',self.NN)
		mlabraw.eval(sess1,'cd '+mlabpath)
		mlabraw.eval(sess1,'cd Toolkit41')
		message = ' '*70
		eval_list = ['message='+"'"+message+"'",
			     'warnings = [];',
			     'options;',
			     'solve;']
		try:
			for x in eval_list:
				mlabraw.eval(sess1,x)
			self.PP = N.matrix(mlabraw.get(sess1,'PP'))
			self.QQ = N.matrix(mlabraw.get(sess1,'QQ'))
			self.RR = N.matrix(mlabraw.get(sess1,'RR'))
			self.SS = N.matrix(mlabraw.get(sess1,'SS'))
			self.WW = N.matrix(mlabraw.get(sess1,'WW'))
		except MyErr:
			pass
		finally:
			return
#----------------------------------------------------------------------------------------------------------------------
class MatKlein:

	def __init__(self,intup):
		self.sess1 = intup[-1]
		self.uhlig = PyUhlig(intup[:-1])
		self.uhlig.mkmats()
		self.AA = self.uhlig.AA
		self.BB = self.uhlig.BB
		self.CC = self.uhlig.CC
		self.DD = self.uhlig.DD
		self.FF = self.uhlig.FF
		self.GG = self.uhlig.GG
		self.HH = self.uhlig.HH
		self.JJ = self.uhlig.JJ
		self.KK = self.uhlig.KK
		self.LL = self.uhlig.LL
		self.MM = self.uhlig.MM
		self.NN = self.uhlig.NN
		self.mkKleinMats()
		self.outdic = {}

	def mkKleinMats(self):
		# Make uhlig matrices locally available for computations
		AA = self.AA
		BB = self.BB
		CC = self.CC
		DD = self.DD
		FF = self.FF
		GG = self.GG
		HH = self.HH
		JJ = self.JJ
		KK = self.KK
		LL = self.LL
		MM = self.MM
		NN = self.NN

		# Determine size of states, endogenous
		exo_st = MAT.shape(NN)[1]
		endo_st = MAT.shape(BB)[1]
		endo_cn = MAT.shape(CC)[1]
		n_deteq = MAT.shape(AA)[0]
		n_expeq = MAT.shape(JJ)[0]
		tot_st = exo_st+endo_st
		self.tstates = tot_st
		tot_var = tot_st+endo_cn

		klein_A_rtwo = MAT.hstack((LL,GG))
		klein_A_rtwo = MAT.hstack((klein_A_rtwo,JJ))
		klein_A_rtwo_rows = MAT.shape(klein_A_rtwo)[0]
		klein_A_rtwo_cols = MAT.shape(klein_A_rtwo)[1]

		klein_A_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))
		klein_A_rone = MAT.hstack((MAT.identity(exo_st),klein_A_rone[:,exo_st:]))

		klein_A_rthree = MAT.hstack((MAT.zeros((n_deteq,exo_st)),AA))
		klein_A_rthree = MAT.hstack((klein_A_rthree,MAT.zeros((n_deteq,endo_cn))))

		klein_A = MAT.vstack((klein_A_rone,klein_A_rtwo))
		klein_A = MAT.vstack((klein_A,klein_A_rthree))

		klein_B_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))

		klein_B_rone = MAT.hstack((NN,klein_B_rone[:,exo_st:]))
		klein_B_rtwo = MAT.hstack((-MM,-HH))
		klein_B_rtwo = MAT.hstack((klein_B_rtwo,-KK))
		klein_B_rthree = MAT.hstack((-DD,-BB))
		klein_B_rthree = MAT.hstack((klein_B_rthree,-CC))

		klein_B = MAT.vstack((klein_B_rone,klein_B_rtwo))
		klein_B = MAT.vstack((klein_B,klein_B_rthree))	

		self.A = klein_A
		self.B = klein_B

	def solve(self):
		A = self.A
		B = self.B
		tstates = self.tstates
		sess1 = self.sess1
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+mlabpath)
		mlabraw.eval(sess1,'cd Klein')
		mlabraw.put(sess1,'AA',A)
		mlabraw.put(sess1,'BB',B)
		mlabraw.put(sess1,'tstates',tstates)
		try:
			mlabraw.eval(sess1,'[F,P,Z11]=solab(AA,BB,tstates)')
			self.F = N.matrix(mlabraw.get(sess1,'F'))
			self.P = N.matrix(mlabraw.get(sess1,'P'))
			self.Z11 = N.matrix(mlabraw.get(sess1,'Z11'))
			self.outdic['P'] = self.P
			self.outdic['F'] = self.F
		except MyErr:
			pass
		finally:
			return
#----------------------------------------------------------------------------------------------------------------------
class MatKleinD:

	def __init__(self,intup):
		self.interm3 = intup[0]
		self.param = intup[1]
		self.sstate_list = intup[2]
		self.varvecs = intup[3]
		self.vardic = intup[4]
		self.vardic2 = intup[5]
		self.shockdic = intup[6]
		self.shockdic2 = intup[7]
		self.varns = intup[8]
		self.re_var = intup[9]
		self.re_var2 = intup[10]
		self.ishockdic = intup[11]
		self.ishockdic2 = intup[12]
		self.sess1 = intup[13]
		self.sstate = intup[14]
		self.mkeqs()
		self.mkvarl()
		self.mkeqs2()
		self.mkgrad()
		self.mkAB()
		self.solab()
		self.mkhess()
		self.solab2()

	def mkeqs(self):
		list_tmp = COP.deepcopy(self.interm3)
		reg2 = RE.compile('\*{2,2}')
		reva2 = self.re_var2
		i1=0
		for x in list_tmp:
			list_tmp[i1]=reg2.sub('^',x.split('=')[0]).strip()
			for y in self.subvars(tuple(self.varvecs['e'])):
				while y in list_tmp[i1]:
					list_tmp[i1] = list_tmp[i1].replace(y,'1')
			i1=i1+1
		self.interm4 = list_tmp

	def mkvarl(self):
		varlist=[]
		nexost = len(self.varvecs['z'])
		nendost = len(self.varvecs['k'])
		nendocon = len(self.varvecs['y'])
		exost_f = makeForward(self.varvecs['z'],'1','0')
		for x in exost_f:
			varlist.append(x)
		for x in self.varvecs['k']:
			varlist.append(x)
		endocons_f = makeForward(self.varvecs['y'],'1','0')
		for x in endocons_f:
			varlist.append(x)
		exost_l = self.varvecs['z']
		for x in exost_l:
			varlist.append(x)
		endost_l = makeLags(self.varvecs['k'],'1')
		for x in endost_l:
			varlist.append(x)
		for x in self.varvecs['y']:
			varlist.append(x)
		varlist2 = list(self.subvars(tuple(varlist)))
		self.vlist = varlist2

		varlist=[]
		for x in self.varvecs['z']:
			varlist.append(x)
		for x in self.varvecs['k']:
			varlist.append(x)
		endocons_f = makeForward(self.varvecs['y'],'1','0')
		for x in endocons_f:
			varlist.append(x)
		exost_l = makeLags(self.varvecs['z'],'1')
		for x in exost_l:
			varlist.append(x)
		endost_l = makeLags(self.varvecs['k'],'1')
		for x in endost_l:
			varlist.append(x)
		for x in self.varvecs['y']:
			varlist.append(x)
		varlist2 = list(self.subvars(tuple(varlist)))
		self.vlist2 = varlist2

	def mkeqs2(self):
		nexost = len(self.varvecs['z'])
		nendost = len(self.varvecs['k'])
		nendocon = len(self.varvecs['y'])
		list_tmp1 = COP.deepcopy(self.interm4)
		eq_len = len(self.interm4)
		sub_list2 = self.vlist2
		sub_list = self.vlist
		reva2 = self.re_var2
		i1=0
		for x in list_tmp1:
			str_tmp1 = x[:]
			i2=0
			for y1,y2 in zip(sub_list,sub_list2):
				if i1 < eq_len-nexost:
					while reva2.search(str_tmp1) != None:
						ma1 = reva2.search(str_tmp1)
						indx1 = sub_list.index(ma1.group(0))
						str_tmp1 = str_tmp1.replace(ma1.group(0),'exp(x0('+str(indx1+1)+'))')[:]
					i2 = i2 + 1
				else:
					while reva2.search(str_tmp1) != None:
						ma2 = reva2.search(str_tmp1)
						indx2 = sub_list2.index(ma2.group(0))
						str_tmp1 = str_tmp1.replace(ma2.group(0),'exp(x0('+str(indx2+1)+'))')[:]
					i2 = i2 + 1
			list_tmp1[i1] = str_tmp1
			i1=i1+1

		self.interm5 = list_tmp1		

	def mkfunc(self):
		rege = RE.compile('\*{2,2}')
		inde = ' '*2
		eq_len=len(self.interm5)
		mfile=open(path_opts['MlabPath']['path']+'/Klein/'+'mfunc.m','w')
		mfile.write('function y = mfunc(x0)\n')
		mfile.write(inde+'%Parameter Values\n')
		for x in self.param.items():
			mfile.write(inde+x[0]+'='+str(x[1])+';\n')
		mfile.write('\n')
		mfile.write(inde+'%Steady State Values\n')
		for x in self.sstate_list:
			x[0]=rege.sub('^',x[0])
			x[1]=rege.sub('^',x[1])
			mfile.write(inde+x[0]+'='+str(x[1])+';\n')
		mfile.write('\n')
		mfile.write(inde+'y=zeros('+str(eq_len)+',1);\n')
		mfile.write('\n')
		for i1,i2 in zip(range(1,eq_len+1,1),self.interm5):
			mfile.write(inde+'y('+str(i1)+')='+i2+';'+'\n')
		mfile.close()

	def mkgrad(self):
		self.mkfunc()
		sess1 = self.sess1
		sstate = self.sstate
		ssdic = sstate[0]
		ssdic.update(sstate[1])
		ssdic2 = {}
		for x in ssdic.keys():
			if x[-4:] == '_bar':
				ssdic2[x] = ssdic[x]
		for x in ssdic2.keys():
			ssdic2[x] = N.log(ssdic2[x])
		vlist2 = []
		for x in self.vlist:
			ma = self.re_var2.search(x)
			tvar = ma.group('var')
			tvar = tvar.upper()
			tvar = tvar+'_bar'
			vlist2.append(tvar)
		invec = []
		for x in vlist2:
			invec.append(ssdic2[x])
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+path_opts['MlabPath']['path'])
		mlabraw.eval(sess1,'cd Klein')
		mlabraw.put(sess1,'x0',invec)
		mlabraw.eval(sess1,"grdm=centgrad('mfunc',x0)")
		gradm = N.matrix(mlabraw.get(sess1,'grdm'))
		self.gradm = gradm

	def mkAB(self):
		nexost = len(self.varvecs['z'])
		nendost = len(self.varvecs['k'])
		nendocon = len(self.varvecs['y'])
		A = self.gradm[:,:nexost+nendost+nendocon]
		B = -self.gradm[:,nexost+nendost+nendocon:]
		self.AA = A
		self.BB = B

	def solab(self):
		A = self.AA
		B = self.BB
		tstates = len(self.varvecs['k'])+len(self.varvecs['z'])
		(F,P,retcon) = isolab(A,B,tstates,MAT.shape(A)[0])
		if MAT.sum(P.reshape(-1,1)) == 0.0:
			return
		else:
			self.P = N.matrix(P)
			self.F = N.matrix(F)

	def mkhess(self):
		self.mkfunc()
		sess1 = self.sess1
		sstate = self.sstate
		ssdic = sstate[0]
		ssdic.update(sstate[1])
		ssdic2 = {}
		for x in ssdic.keys():
			if x[-4:] == '_bar':
				ssdic2[x] = ssdic[x]
		for x in ssdic2.keys():
			ssdic2[x] = N.log(ssdic2[x])
		vlist2 = []
		for x in self.vlist:
			ma = self.re_var2.search(x)
			tvar = ma.group('var')
			tvar = tvar.upper()
			tvar = tvar+'_bar'
			vlist2.append(tvar)
		invec = []
		for x in vlist2:
			invec.append(ssdic2[x])
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+path_opts['MlabPath']['path'])
		mlabraw.eval(sess1,'cd Klein')
		mlabraw.put(sess1,'x0',invec)
		mlabraw.eval(sess1,"hessm=centhess('mfunc',x0)")
		hessm = N.matrix(mlabraw.get(sess1,'hessm'))
		self.hessm = hessm

	def solab2(self):
		sess1 = self.sess1
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+path_opts['MlabPath']['path'])
		mlabraw.eval(sess1,'cd Klein')
		mlabraw.put(sess1,'gradm',self.gradm)
		mlabraw.put(sess1,'hessm',self.hessm)
#----------------------------------------------------------------------------------------------------------------------
class MatWood:

	def __init__(self,intup):
		self.jAA = intup[0]
		self.jBB = intup[1]
		self.nexo = intup[2]
		self.ncon = intup[3]
		self.nendo = intup[4]
		self.sess1 = intup[5]
		self.NY = self.ncon+self.nendo
		self.NK = self.nendo
		self.NX = self.nexo
		self.mkwoodm()

	def solve(self):
		self.reds()
		self.solds()

	def mkwoodm(self):
		AA = self.jAA
		BB = self.jBB
		nendo = self.nendo
		nexo = self.nexo
		ncon = self.ncon
		nstates = nexo+nendo

		wAA = MAT.hstack((AA[:-nexo,nstates:],AA[:-nexo,nexo:nstates]))
		#wAA = MAT.hstack((AA[:,nstates:],AA[:,:nstates]))
		wBB = MAT.hstack((BB[:-nexo,nstates:],BB[:-nexo,nexo:nstates]))
		#wBB = MAT.hstack((BB[:,nstates:],BB[:,:nstates]))
		wCC = BB[:-nexo,:nexo]
		#wCC = MAT.zeros((ncon+nstates,1))


		self.wAA = wAA
		self.wBB = wBB
		self.wCC = wCC

	def reds(self):
		A = self.wAA
		B = self.wBB
		C = self.wCC
		NX = self.NX
		NK = self.NK
		NY = self.NY
		sess1 = self.sess1
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+mlabpath)
		mlabraw.eval(sess1,'cd Woodford')
		mlabraw.put(sess1,'wAA',A)
		mlabraw.put(sess1,'wBB',B)
		mlabraw.put(sess1,'wCC',C)
		mlabraw.put(sess1,'NX',NX)
		mlabraw.put(sess1,'NY',NY)
		mlabraw.put(sess1,'NK',NK)
		mlabraw.eval(sess1,'[Br,Cr,Lr,NF] = redsf(wAA,wBB,wCC,NY,NX,NK)')
		self.Br = N.matrix(mlabraw.get(sess1,'Br'))
		self.Cr = N.matrix(mlabraw.get(sess1,'Cr'))
		self.Lr = N.matrix(mlabraw.get(sess1,'Lr'))
		self.NF = N.matrix(mlabraw.get(sess1,'NF'))

	def solds(self):
		Br = self.Br
		Cr = self.Cr
		Lr = self.Lr
		NF = self.NF
		NX = self.NX
		NK = self.NK
		NY = self.NY
		sess1 = self.sess1
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+mlabpath)
		mlabraw.eval(sess1,'cd Woodford')
		mlabraw.put(sess1,'Br',Br)
		mlabraw.put(sess1,'Cr',Cr)
		mlabraw.put(sess1,'Lr',Lr)
		mlabraw.put(sess1,'NF',NF)
		mlabraw.put(sess1,'NX',NX)
		mlabraw.put(sess1,'NY',NY)
		mlabraw.put(sess1,'NK',NK)
		mlabraw.eval(sess1,'[D,F,G,H] = soldsf(Br,Cr,Lr,NY,NX,NK,NF)')
		self.D = N.matrix(mlabraw.get(sess1,'D'))
		self.F = N.matrix(mlabraw.get(sess1,'F'))
		self.G = N.matrix(mlabraw.get(sess1,'G'))
		self.H = N.matrix(mlabraw.get(sess1,'H'))
#----------------------------------------------------------------------------------------------------------------------
class ForKlein:

	def __init__(self,intup):
		self.uhlig = PyUhlig(intup)
		self.uhlig.mkmats()
		self.AA = self.uhlig.AA
		self.BB = self.uhlig.BB
		self.CC = self.uhlig.CC
		self.DD = self.uhlig.DD
		self.FF = self.uhlig.FF
		self.GG = self.uhlig.GG
		self.HH = self.uhlig.HH
		self.JJ = self.uhlig.JJ
		self.KK = self.uhlig.KK
		self.LL = self.uhlig.LL
		self.MM = self.uhlig.MM
		self.NN = self.uhlig.NN
		self.mkKleinMats()
		self.outdic = {}

	def mkKleinMats(self):
		# Make uhlig matrices locally available for computations
		AA = self.AA
		BB = self.BB
		CC = self.CC
		DD = self.DD
		FF = self.FF
		GG = self.GG
		HH = self.HH
		JJ = self.JJ
		KK = self.KK
		LL = self.LL
		MM = self.MM
		NN = self.NN

		# Determine size of states, endogenous
		exo_st = MAT.shape(NN)[1]
		endo_st = MAT.shape(BB)[1]
		endo_cn = MAT.shape(CC)[1]
		n_deteq = MAT.shape(AA)[0]
		n_expeq = MAT.shape(JJ)[0]
		tot_st = exo_st+endo_st
		self.tstates = tot_st
		tot_var = tot_st+endo_cn

		klein_A_rtwo = MAT.hstack((LL,GG))
		klein_A_rtwo = MAT.hstack((klein_A_rtwo,JJ))
		klein_A_rtwo_rows = MAT.shape(klein_A_rtwo)[0]
		klein_A_rtwo_cols = MAT.shape(klein_A_rtwo)[1]

		klein_A_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))
		klein_A_rone = MAT.hstack((MAT.identity(exo_st),klein_A_rone[:,exo_st:]))

		klein_A_rthree = MAT.hstack((MAT.zeros((n_deteq,exo_st)),AA))
		klein_A_rthree = MAT.hstack((klein_A_rthree,MAT.zeros((n_deteq,endo_cn))))

		klein_A = MAT.vstack((klein_A_rone,klein_A_rtwo))
		klein_A = MAT.vstack((klein_A,klein_A_rthree))

		klein_B_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))

		klein_B_rone = MAT.hstack((NN,klein_B_rone[:,exo_st:]))
		klein_B_rtwo = MAT.hstack((-MM,-HH))
		klein_B_rtwo = MAT.hstack((klein_B_rtwo,-KK))
		klein_B_rthree = MAT.hstack((-DD,-BB))
		klein_B_rthree = MAT.hstack((klein_B_rthree,-CC))

		klein_B = MAT.vstack((klein_B_rone,klein_B_rtwo))
		klein_B = MAT.vstack((klein_B,klein_B_rthree))	

		self.A = klein_A
		self.B = klein_B

	def solve(self):
		A = self.A
		B = self.B
		tstates = MAT.shape(self.AA)[1] + MAT.shape(self.DD)[1]
		(F,P,retcon) = isolab(A,B,tstates,MAT.shape(A)[0])
		if MAT.sum(P.reshape(-1,1)) == 0.0:
			return
		else:
			self.P = N.matrix(P)
			self.F = N.matrix(F)
#----------------------------------------------------------------------------------------------------------------------
class PyKlein2D:

	def __init__(self,intup):
		self.gra = intup[0]
		self.hes = intup[1]
		self.nendo = intup[2]
		self.nexo = intup[3]
		self.ncon = intup[4]
		self.sigma = intup[5]
		self.A = intup[6]
		self.B = intup[7]
		self.vardic = intup[8]
		self.vdic = intup[9]
		self.modname = intup[10]
		self.audic = intup[11]
		self.nacon = len(self.audic['con']['var'])
		self.naendo = len(self.audic['endo']['var'])
		if self.vardic['other']['var']:
			self.numjl = intup[12]
			self.numhl = intup[13]
			self.nother = intup[14]
			self.oswitch = 1
			kintup = (self.gra,self.nendo,
				  self.nexo,self.ncon,
				  self.sigma,self.A,self.B,
				  self.vardic,self.vdic,
				  self.modname,self.audic,
				  self.numjl,self.nother)
		else:
			self.oswitch = 0
			kintup = (self.gra,self.nendo,
				  self.nexo,self.ncon,
				  self.sigma,self.A,self.B,
				  self.vardic,self.vdic,
				  self.modname,self.audic)
		self.tstates = self.nendo+self.nexo
		self.forkleind = ForKleinD(kintup)
		self.ssigma = self.mkssigma(self.tstates,self.nexo,self.sigma)

	def mkssigma(self,tstates,nexo,sigma):
		ssigma = MAT.zeros((tstates,tstates))
		for i in range(nexo):
			ssigma[i,i] = sigma[i,i]
		return ssigma

	def solab2(self):

		# Tracem function
		def tracem(xmat):
			n = MAT.shape(xmat)[1]
			m = MAT.shape(xmat)[0]/n
			ymat = MAT.zeros((m,1))
			for i1 in xrange(int(m)):
				ymat[i1,0]=MAT.trace(xmat[n*i1:i1*n+1,:n])
			return ymat

		ssigma = self.ssigma
		pp = self.P
		ff = self.F
		gra = self.gra
		hes = self.hes
		m = self.ncon+self.nendo+self.nexo
		ny = self.ncon
		nx = self.tstates

		f1 = gra[:,:nx]
		f2 = gra[:,nx:m]
		f4 = gra[:,m+nx:2*m]

		mm = MAT.vstack((pp,ff*pp,MAT.eye(nx),ff))
		aa1 = MAT.kron(MAT.eye(m),mm.T)*hes*mm
		bb1 = MAT.kron(f1,MAT.eye(nx))
		bb2 = MAT.kron(f2,MAT.eye(nx))
		bb4 = MAT.kron(f4,MAT.eye(nx))
		cc1 = MAT.kron(MAT.eye(ny),pp.T)
		cc2 = MAT.kron(ff,MAT.eye(nx))

		aa_1 = MAT.kron(MAT.eye(nx),bb4)+MAT.kron(pp.T,bb2*cc1)
		aa_2 = MAT.kron(MAT.eye(nx),bb1+bb2*cc2)
		aa = MAT.hstack((aa_1,aa_2))
		sol = N.linalg.solve(-aa,HLP.ctr(aa1.flatten(1)))
		ee = HLP.ctr(sol[:nx**2*ny].reshape(-1,nx*ny))
		gg = HLP.ctr(sol[nx**2*ny:].reshape(-1,nx**2))
		ma = 2.0*MAT.hstack((f1+f2*ff,f2+f4))
		eyeff = MAT.vstack((MAT.eye(nx),ff,MAT.zeros((m,nx))))
		ve_1 = f2*tracem(MAT.kron(MAT.eye(ny),ssigma)*ee)
		ve_2 = tracem(MAT.kron(MAT.eye(m),HLP.ctr(eyeff))*hes*eyeff*ssigma)
		ve = ve_1 + ve_2
		kxy = -(ma.I)*ve
		kx = kxy[:nx]
		ky = kxy[nx:]
		self.E = ee
		self.G = gg
		self.KX = kx
		self.KY = ky

	def solve(self):
		self.forkleind.solve()
		self.P = self.forkleind.P
		self.F = self.forkleind.F
		self.solab2()

	def sim(self,tlen,sntup=None):
		# Add 1000 observations, as they will be thrown away
		# Add one more observation to start first-order vector
		exoli = [x[1] for x in self.vardic['exo']['var']]
		indx=[]
		if sntup == None:
			indx = range(len(exoli))
		else:
			for name in sntup:
				if name not in exoli:
					return 'Error: '+name+' is not a valid exo variable for this model!'
				else:
					indx.append(exoli.index(name))
		indx.sort()
		ncon = self.ncon
		nexo = self.nexo
		nendo = self.nendo
		tstates = self.tstates
		tlena = 1000+tlen
		sigma = self.sigma
		ssigma = self.ssigma
		if self.oswitch:
			numjl = self.numjl
			numhl = self.numhl
		kx = self.KX
		ky = self.KY
		pp = self.P
		ff = self.F
		ee = self.E
		gg = self.G
		count = 0
		for varia in MAT.diag(sigma):
			if locals().has_key('ranvec'):
				if count in indx:
					ranvec = MAT.vstack((ranvec,N.sqrt(varia)*MAT.matrix(N.random.standard_normal(tlena))))
				else:
					ranvec = MAT.vstack((ranvec,MAT.zeros((1,tlena))))
			else:
				if count in indx:
					ranvec = N.sqrt(varia)*MAT.matrix(N.random.standard_normal(tlena))
				else:
					ranvec = MAT.zeros((1,tlena))
			count = count + 1

		ranvec = MAT.vstack((ranvec,MAT.zeros((nendo,tlena))))

		x_one_m1 = ranvec[:,0]
		x_one_0 = pp*x_one_m1
		x_two_m1 = kx + ranvec[:,0]
		y_one_0 = ff*x_one_m1
		y_one_m1 = MAT.zeros(y_one_0.shape)
		y_two_0 = ky+0.5*MAT.kron(MAT.eye(ncon),x_one_m1.T)*ee*x_one_m1
		if self.oswitch:
			nother = self.nother
			o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))
			o_two_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))+\
				0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1)).T)*\
				numhl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

		x_one_c = COP.deepcopy(x_one_m1)
		y_one_c = COP.deepcopy(y_one_0)
		x_two_c = COP.deepcopy(x_two_m1)
		y_two_c = COP.deepcopy(y_two_0)
		if self.oswitch:
			o_one_c = COP.deepcopy(o_one_0)
			o_two_c = COP.deepcopy(o_two_0)

		x_one = COP.deepcopy(x_one_c)
		y_one = COP.deepcopy(y_one_c)
		x_two = COP.deepcopy(x_two_c)
		y_two = COP.deepcopy(y_two_c)
		if self.oswitch:
			o_one = COP.deepcopy(o_one_c)
			o_two = COP.deepcopy(o_two_c)

		for i1 in range(1,ranvec.shape[1],1):

			x_one_n = pp*x_one_c+ranvec[:,i1]
			y_one_n = ff*x_one_c

			x_two_n = kx+pp*x_two_c+\
				0.5*MAT.kron(MAT.eye(tstates),x_one_c.T)*\
				gg*x_one_c+ranvec[:,i1]

			y_two_n = ky+ff*x_two_c+\
				0.5*MAT.kron(MAT.eye(ncon),x_one_c.T)*\
				ee*x_one_c

			if self.oswitch:
				o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))
				o_two_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))+\
					0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c)).T)*\
					numhl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

			x_one = MAT.hstack((x_one,x_one_c))
			y_one = MAT.hstack((y_one,y_one_n))
			x_two = MAT.hstack((x_two,x_two_c))
			y_two = MAT.hstack((y_two,y_two_n))
			if self.oswitch:
				o_one = MAT.hstack((o_one,o_one_n))
				o_two = MAT.hstack((o_two,o_two_n))

			x_one_c = x_one_n
			y_one_c = y_one_n
			x_two_c = x_one_n
			y_two_c = y_two_n
			if self.oswitch:
				o_one_c = o_one_n
				o_two_c = o_two_n

		# Throw away first 1000 obs
		x_one = x_one[:,1000:]
		y_one = y_one[:,1000:]
		y_two = y_two[:,1000:]
		x_two = x_two[:,1000:]
		if self.oswitch:
			o_one = o_one[:,1000:]
			o_two = o_two[:,1000:]

		self.sim_y_two = y_two
		self.sim_y_one = y_one
		self.sim_x_two = x_two
		self.sim_x_one = x_one
		self.insim = [self.sim_x_two,self.sim_y_two]
		if self.oswitch:
			self.sim_o_one = o_one
			self.sim_o_two = o_two
			self.insim = self.insim + [self.sim_o_two,]

	def mkcvar(self,vname=False,insim='insim'):
		# Now also produce unconditional variances, and normal. variances with v
		# Check if simul data is available
		if insim not in dir(self):
			return "Error: Can't produce uncond. variance table without simulation data!"

		exoli = [x[1] for x in self.vardic['exo']['var']]
		endoli = [x[1] for x in self.vardic['endo']['var']]
		stateli = exoli + endoli
		conli = [x[1] for x in self.vardic['con']['var']]
		alli = stateli+conli
		if self.oswitch:
			otherli = [x[1] for x in self.vardic['other']['var']]
			alli = alli + otherli

		# Check if vname is valid
		if vname not in alli:
			return 'Error: '+vname+' is not a valid variable for this model!'
		insim = eval('self.'+insim)
		osim_x = COP.deepcopy(insim[0])
		osim_y = COP.deepcopy(insim[1])
		nendo = self.nendo
		nexo = self.nexo
		ncon = self.ncon

		# Now hp filter the simulations before graphing according to filtup
		for i1 in xrange(osim_y.shape[0]):
			yy = osim_y[i1,:].__array__().T
			woy = N.zeros((osim_y.shape[1],3))
			lam = 1600
			yyf = MAT.matrix(hpfilt(yy,woy,osim_y.shape[1],1600,0))
			osim_y[i1,:] = yyf
		# Now filter the state variables!
		for i1 in xrange(osim_x.shape[0]):
			xx = osim_x[i1,:].__array__().T
			wox = N.zeros((osim_x.shape[1],3))
			lam = 1600
			xxf = MAT.matrix(hpfilt(xx,wox,osim_x.shape[1],1600,0))
			osim_x[i1,:] = xxf

		if self.oswitch:
			osim_o = COP.deepcopy(insim[2])
			# Now hp filter the other variables!
			for i1 in xrange(osim_o.shape[0]):
				oo = osim_o[i1,:].__array__().T
				woo = N.zeros((osim_o.shape[1],3))
				lam = 1600
				oof = MAT.matrix(hpfilt(oo,woo,osim_o.shape[1],1600,0))
				osim_o[i1,:] = oof

		# Stack all matrices into one
		allmat = MAT.vstack((osim_x,osim_y))
		if self.oswitch:
			allmat = MAT.vstack((allmat,osim_o))
		varmat = MAT.zeros((allmat.shape[0],1))
		for i1 in range(allmat.shape[0]):
			varmat[i1,0] = S.stats.var(allmat[i1,:].A.flatten(1))

		if vname:
			relmat = MAT.zeros((allmat.shape[0],1))
			indx = alli.index(vname)
			varvari = S.stats.var(allmat[indx,:].A.flatten(1))
			for i1 in range(allmat.shape[0]):
				relmat[i1,0] = varmat[i1,0]/varvari
			self.cvarm = MAT.hstack((varmat,relmat))
		else:
			self.cvarm = varmat





	def mkact(self,vname,intup=None,insim='insim'):
		# Now also produce autocorrelation table
		# Check if simul data is available
		if insim not in dir(self):
			return "Error: Can't produce autocorrelation table without simulation data!"

		exoli = [x[1] for x in self.vardic['exo']['var']]
		endoli = [x[1] for x in self.vardic['endo']['var']]
		stateli = exoli + endoli
		conli = [x[1] for x in self.vardic['con']['var']]
		alli = stateli+conli
		if self.oswitch:
			otherli = [x[1] for x in self.vardic['other']['var']]
			alli = alli + otherli
		# Check if vname is valid
		if vname not in alli:
			return 'Error: '+vname+' is not a valid variable for this model!'

		if intup == None:
			lags = 3
			leads = 3
		else:
			lags = intup[0]
			leads = intup[1]

		insim = eval('self.'+insim)
		osim_x = COP.deepcopy(insim[0])
		osim_y = COP.deepcopy(insim[1])

		# Now hp filter the simulations before graphing according to filtup
		for i1 in xrange(osim_y.shape[0]):
			yy = osim_y[i1,:].__array__().T
			woy = N.zeros((osim_y.shape[1],3))
			lam = 1600
			yyf = MAT.matrix(hpfilt(yy,woy,osim_y.shape[1],1600,0))
			osim_y[i1,:] = yyf
		# Now filter the state variables!
		for i1 in xrange(osim_x.shape[0]):
			xx = osim_x[i1,:].__array__().T
			wox = N.zeros((osim_x.shape[1],3))
			lam = 1600
			xxf = MAT.matrix(hpfilt(xx,wox,osim_x.shape[1],1600,0))
			osim_x[i1,:] = xxf

		sim_x = osim_x[:,lags:-leads]
		sim_xf = osim_x[:,leads+1:]
		sim_yf = osim_y[:,leads+1:]
		sim_y = osim_y[:,lags:-leads]
		sim_xl = osim_x[:,:-lags-1]
		sim_yl = osim_y[:,:-lags-1]

		if self.oswitch:
			osim_o = COP.deepcopy(insim[2])
			# Now hp filter the other variables!
			for i1 in xrange(osim_o.shape[0]):
				oo = osim_o[i1,:].__array__().T
				woo = N.zeros((osim_o.shape[1],3))
				lam = 1600
				oof = MAT.matrix(hpfilt(oo,woo,osim_o.shape[1],1600,0))
				osim_o[i1,:] = oof
			sim_o = osim_o[:,lags:-leads]
			sim_of = osim_o[:,leads+1:]
			sim_ol = osim_o[:,:-lags-1]
			alli = alli + otherli

		posa = alli.index(vname)
		actm = MAT.zeros((len(alli),lags+1+leads))
		sima = MAT.vstack((sim_x,sim_y))
		sima_f = MAT.vstack((sim_xf,sim_yf))
		sima_l = MAT.vstack((sim_xl,sim_yl))
		if self.oswitch:
			sima = MAT.vstack((sima,sim_o))
			sima_f = MAT.vstack((sima_f,sim_of))
			sima_l = MAT.vstack((sima_l,sim_ol))
		for i1 in range(sima.shape[0]):
			# Do lags
			for i2 in range(lags):
				actm[i1,lags-i2-1] = N.round(STA.corrcoef(sima[posa,:].A.flatten(1),sima_l[i1,i2:sima_l.shape[1]-lags+i2+1].A.flatten(1))[1,0],3)
			# Do current
			actm[i1,lags] = N.round(STA.corrcoef(sima[posa,:].A.flatten(1),sima[i1,:].A.flatten(1))[1,0],3)
			# Do leads
			for i2 in range(leads):
				actm[i1,lags+1+i2] = N.round(STA.corrcoef(sima[posa,:].A.flatten(1),sima_f[i1,i2:sima_f.shape[1]-leads+i2+1].A.flatten(1))[1,0],3)

		self.actm = actm
		self.actname = vname
		if intup:
			self.actin = intup

	def show_act(self):
		if 'actname' not in dir(self):
			return 'Error: You have not produced any autocorrelation table yet!'
		nexo = self.nexo
		nendo = self.nendo
		ncon = self.ncon
		naendo = self.naendo
		nacon = self.nacon
		nother = len(self.vardic['other']['var'])
		respva = COP.deepcopy(self.actname)
		if 'actin' in dir(self):
			lags = COP.deepcopy(int(self.actin[0]))
			leads = COP.deepcopy(int(self.actin[1]))
		else:
			lags = 3
			leads = 3
		actm = COP.deepcopy(self.actm)
		vdic = self.vdic
		allvari = [x[1] for x in vdic['exo']]
		allvari = allvari + [x[1] for x in vdic['endo']]
		allvari = allvari + [x[1] for x in vdic['con']]
		if self.oswitch:
			allvari = allvari + [x[1] for x in vdic['other']]
		allvari = allvari[:nexo+nendo-naendo]+allvari[nexo+nendo:nexo+nendo+ncon-nacon]
		actm = MAT.vstack((actm[:nexo+nendo-naendo,:],actm[nexo+nendo:nexo+nendo+ncon-nacon,:]))

		# determine longest varname
		vnmax = 0
		for vari in allvari:
			if len(vari) > vnmax:
				vnmax = len(vari)
		print ''
		print 'Autocorrelation table, current '+respva
		print '='*65
		for i1,vari in enumerate(allvari):
			print vari+(vnmax-len(vari))*' '+'  |'+str(actm[i1,:])[2:-2]

	def show_sim(self,intup,filtup=None,insim='insim'):
		# Check if simulations have been carried out
		if insim not in dir(self):
			return 'Error: You have not produced any simulations yet! Nothing to graph!'
		mname = self.modname
		vardic = self.vardic
		insim = eval('self.'+insim)
		sim_x = COP.deepcopy(insim[0])
		sim_y = COP.deepcopy(insim[1])
		tlen = sim_x.shape[1]
		if self.oswitch:
			sim_o = COP.deepcopy(insim[2])
			nother = self.nother
			otherli = [x[1] for x in vardic['other']['var']]
		nexo = self.nexo
		nendo = self.nendo
		ncon = self.ncon
		conli = [x[1] for x in vardic['con']['var']]
		endoli = [x[1] for x in vardic['endo']['var']]
		exoli = [x[1] for x in vardic['exo']['var']]

		# If filtup is None, create filtup from vardic
		if filtup == None:
			filli = [0]*len(intup)
			for i1,x in enumerate(intup):
				for y in vardic.keys():
					if x in [z[1] for z in vardic[y]['var']]:
						indx = [z[1] for z in vardic[y]['var']].index(x)
						if 'hp' in [z for z in vardic[y]['mod']][indx]:
							filli[i1] = 1
			filtup = tuple(filli)


		stateli = exoli+endoli
		alli = stateli+conli
		if self.oswitch:
			alli = alli+otherli
		# Check if all name in intup are available to graph
		for name in intup:
			if name not in alli:
				return 'Error: '+name+' is not a valid variable name for this model!'
		# Create x, y and z indeces
		indx = []
		indy = []
		indo = []
		for name in intup:
			if name in stateli:
				indx.append(stateli.index(name))
			elif name in conli:
				indy.append(conli.index(name))
			elif self.oswitch and name in otherli:
				indo.append(otherli.index(name))
		leg = []
		indx.sort()
		indy.sort()
		indo.sort()		

		# Now hp filter the simulations before graphing according to filtup
		for i1 in xrange(sim_y.shape[0]):
			if indy and (i1 in indy) and filtup[list(intup).index(conli[i1])]:
				conli[i1] = conli[i1]+'(hpf)'
				yy = sim_y[i1,:].__array__().T
				woy = N.zeros((tlen,3))
				lam = 1600
				yyf = MAT.matrix(hpfilt(yy,woy,tlen,1600,0))
				sim_y[i1,:] = yyf
		# Now filter the state variables!
		for i1 in xrange(sim_x.shape[0]):
			if indx and (i1 in indx) and filtup[list(intup).index(stateli[i1])]:
				stateli[i1] = stateli[i1]+'(hpf)'
				xx = sim_x[i1,:].__array__().T
				wox = N.zeros((tlen,3))
				lam = 1600
				xxf = MAT.matrix(hpfilt(xx,wox,tlen,1600,0))
				sim_x[i1,:] = xxf
			# Now hp filter the other variables!
		if self.oswitch:
			for i1 in xrange(sim_o.shape[0]):
				if indo and (i1 in indo) and filtup[list(intup).index(otherli[i1])]:
					otherli[i1] = otherli[i1]+'(hpf)'
					oo = sim_o[i1,:].__array__().T
					woo = N.zeros((tlen,3))
					lam = 1600
					oof = MAT.matrix(hpfilt(oo,woo,tlen,1600,0))
					sim_o[i1,:] = oof

		if indx and indy and indo:
			for x in indx:
				leg.append(stateli[x])
			for y in indy:
				leg.append(conli[y])
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((sim_x.T[:,indx],sim_y.T[:,indy],sim_o.T[:,indo])).A)
			P.title(str(tlen)+' simulated periods, '+mname)
			P.xlabel('Time')
			P.legend(leg)
			P.show()
		elif not indx and indy and indo:
			for y in indy:
				leg.append(conli[y])
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((sim_y.T[:,indy],sim_o.T[:,indo])).A)
			P.title(str(tlen)+' simulated periods, '+mname)
			P.xlabel('Time')
			P.legend(leg)
			P.show()
		elif indx and not indy and indo:
			for x in indx:
				leg.append(stateli[x])
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((sim_x.T[:,indx],sim_o.T[:,indo])).A)
			P.title(str(tlen)+' simulated periods, '+mname)
			P.xlabel('Time')
			P.legend(leg)
			P.show()
		elif indx and indy and not indo:
			for x in indx:
				leg.append(stateli[x])
			for y in indy:
				leg.append(conli[y])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((sim_x.T[:,indx],sim_y.T[:,indy])).A)
			P.title(str(tlen)+' simulated periods, '+mname)
			P.xlabel('Time')
			P.legend(leg)
			P.show()
		elif indx and not indy and not indo:
			for x in indx:
				leg.append(stateli[x])
			leg = tuple(leg)
			P.figure()
			P.plot(sim_x.T[:,indx].A)
			P.title(str(tlen)+' simulated periods, '+mname)
			P.xlabel('Time')
			P.legend(leg)
			P.show()
		elif not indx and indy and not indo:
			for y in indy:
				leg.append(conli[y])
			leg = tuple(leg)
			P.figure()
			P.plot(sim_y.T[:,indy].A)
			P.title(str(tlen)+' simulated periods, '+mname)
			P.xlabel('Time')
			P.legend(leg)
			P.show()
		elif not indx and not indy and indo:
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(sim_o.T[:,indo].A)
			P.title(str(tlen)+' simulated periods, '+mname)
			P.xlabel('Time')
			P.legend(leg)
			P.show()

	def irf(self,tlen,sntup):
		tlen = tlen + 1
		ncon = self.ncon
		nexo = self.nexo
		nendo = self.nendo
		tstates = self.tstates
		sigma = self.sigma
		ssigma = self.ssigma
		if self.oswitch:
			numjl = self.numjl
			numhl = self.numhl
		kx = self.KX
		ky = self.KY
		pp = self.P
		ff = self.F
		ee = self.E
		gg = self.G

		sposli=[]
		exoli = [x[1] for x in self.vardic['exo']['var']]
		# Check if names are valid
		for name in sntup:
			if name not in exoli:
				return 'Error: '+name+' is not a valid exoshock name for this model!'
		for name in sntup:
			sposli.append(exoli.index(name))
		# Expose spos choice to self for show_irf
		self.spos = (sntup,sposli)

		shock = MAT.zeros((nexo,1))
		sendo = MAT.zeros((nendo,1))

		for spos in sposli:
			shock[spos,0] = 1.0
		shock = MAT.vstack((shock,sendo))

		x_one_m1 = shock
		x_one_0 = pp*x_one_m1
		x_two_m1 = shock
		y_one_0 = ff*x_one_m1
		y_one_m1 = MAT.zeros(y_one_0.shape)
		y_two_0 = 0.5*MAT.kron(MAT.eye(ncon),x_one_m1.T)*ee*x_one_m1
		if self.oswitch:
			nother = self.nother
			o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))
			o_two_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))+\
				0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1)).T)*\
				numhl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

		x_one_c = COP.deepcopy(x_one_m1)
		y_one_c = COP.deepcopy(y_one_0)
		x_two_c = COP.deepcopy(x_two_m1)
		y_two_c = COP.deepcopy(y_two_0)
		if self.oswitch:
			o_one_c = COP.deepcopy(o_one_0)
			o_two_c = COP.deepcopy(o_two_0)

		x_one = COP.deepcopy(x_one_c)
		y_one = COP.deepcopy(y_one_c)
		x_two = COP.deepcopy(x_two_c)
		y_two = COP.deepcopy(y_two_c)
		if self.oswitch:
			o_one = COP.deepcopy(o_one_c)
			o_two = COP.deepcopy(o_two_c)

		for i1 in range(tlen):
			x_one_n = pp*x_one_c
			y_one_n = ff*x_one_c
			x_two_n = pp*x_two_c+0.5*MAT.kron(MAT.eye(tstates),x_one_c.T)*gg*x_one_c
			y_two_n = ff*x_two_c+0.5*MAT.kron(MAT.eye(ncon),x_one_c.T)*ee*x_one_c
			if self.oswitch:
				o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))
				o_two_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))+\
					0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c)).T)*\
					numhl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

			x_one = MAT.hstack((x_one,x_one_c))
			y_one = MAT.hstack((y_one,y_one_n))
			x_two = MAT.hstack((x_two,x_two_c))
			y_two = MAT.hstack((y_two,y_two_n))
			if self.oswitch:
				o_one = MAT.hstack((o_one,o_one_n))
				o_two = MAT.hstack((o_two,o_two_n))

			x_one_c = COP.deepcopy(x_one_n)
			y_one_c = COP.deepcopy(y_one_n)
			x_two_c = COP.deepcopy(x_one_n)
			y_two_c = COP.deepcopy(y_two_n)

		self.irf_x_two = x_two[:,1:-1]
		self.irf_y_two = y_two[:,1:-1]
		self.inirf = [self.irf_x_two,self.irf_y_two]
		if self.oswitch:
			self.irf_o_one = o_one[:,2:]
			self.irf_o_two = o_two[:,2:]
			self.inirf = self.inirf + [self.irf_o_two,]

	def show_irf(self,intup,inirf='inirf'):
		# Check if simulations have been carried out
		if inirf not in dir(self):
			return 'Error: You have not produced any IRFs yet! Nothing to graph!'
		inirf = eval('self.'+inirf)
		irf_x = COP.deepcopy(inirf[0])
		irf_y = COP.deepcopy(inirf[1])
		mname = self.modname
		vardic = self.vardic
		tlen = irf_x.shape[1]
		if self.oswitch:
			irf_o = COP.deepcopy(inirf[2])
			nother = self.nother
			varother = self.vardic['other']['var']
			otherli = [x[1] for x in varother]
		nexo = self.nexo
		nendo = self.nendo
		ncon = self.ncon
		conli = [x[1] for x in vardic['con']['var']]
		endoli = [x[1] for x in vardic['endo']['var']]
		exoli = [x[1] for x in vardic['exo']['var']]

		stateli = exoli+endoli
		alli = stateli+conli
		if self.oswitch:
			alli = alli+otherli
		# Check if all name in intup are available to graph
		for name in intup:
			if name not in alli:
				return 'Error: '+name+' is not a valid variable name for this model!'
		# Create x, y and z indeces
		indx = []
		indy = []
		indo = []
		for name in intup:
			if name in stateli:
				indx.append(stateli.index(name))
			elif name in conli:
				indy.append(conli.index(name))
			elif self.oswitch and name in otherli:
				indo.append(otherli.index(name))
		leg = []
		indx.sort()
		indy.sort()
		indo.sort()

		if indx and indy and indo:
			for x in indx:
				leg.append(stateli[x])
			for y in indy:
				leg.append(conli[y])
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((irf_x.T[:,indx],irf_y.T[:,indy],irf_o.T[:,indo])).A)
			P.title(str(tlen)+' simulated IRF periods, '+mname)
			P.xlabel('Time')
			P.ylabel('Log-Dev from SS')
			P.legend(leg)
			P.show()
		elif not indx and indy and indo:
			for y in indy:
				leg.append(conli[y])
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((irf_y.T[:,indy],irf_o.T[:,indo])).A)
			P.title(str(tlen)+' simulated IRF periods, '+mname)
			P.xlabel('Time')
			P.ylabel('Log-Dev from SS')
			P.legend(leg)
			P.show()
		elif indx and not indy and indo:
			for x in indx:
				leg.append(stateli[x])
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((irf_x.T[:,indx],irf_o.T[:,indo])).A)
			P.title(str(tlen)+' simulated IRF periods, '+mname)
			P.xlabel('Time')
			P.ylabel('Log-Dev from SS')
			P.legend(leg)
			P.show()
		elif indx and indy and not indo:
			for x in indx:
				leg.append(stateli[x])
			for y in indy:
				leg.append(conli[y])
			leg = tuple(leg)
			P.figure()
			P.plot(MAT.hstack((irf_x.T[:,indx],irf_y.T[:,indy])).A)
			P.title(str(tlen)+' simulated IRF periods, '+mname)
			P.xlabel('Time')
			P.ylabel('Log-Dev from SS')
			P.legend(leg)
			P.show()
		elif indx and not indy and not indo:
			for x in indx:
				leg.append(stateli[x])
			leg = tuple(leg)
			P.figure()
			P.plot(irf_x.T[:,indx].A)
			P.title(str(tlen)+' simulated IRF periods, '+mname)
			P.xlabel('Time')
			P.ylabel('Log-Dev from SS')
			P.legend(leg)
			P.show()
		elif not indx and indy and not indo:
			for y in indy:
				leg.append(conli[y])
			leg = tuple(leg)
			P.figure()
			P.plot(irf_y.T[:,indy].A)
			P.title(str(tlen)+' simulated IRF periods, '+mname)
			P.xlabel('Time')
			P.ylabel('Log-Dev from SS, hp-filtered')
			P.legend(leg)
			P.show()
		elif not indx and not indy and indo:
			for o in indo:
				leg.append(otherli[o])
			leg = tuple(leg)
			P.figure()
			P.plot(irf_o.T[:,indo].A)
			P.title(str(tlen)+' simulated IRF periods, '+mname)
			P.xlabel('Time')
			P.ylabel('Log-Dev from SS, hp-filtered')
			P.legend(leg)
			P.show()

	def tester(self):
		import sys
		from wx.lib.mixins.inspection import InspectableApp
		app = InspectableApp(False)
		frame = TestFrame(None)
		frame.Show(True)
		app.MainLoop()
#----------------------------------------------------------------------------------------------------------------------
class MatKlein2D(PyKlein2D):

	def __init__(self,intup):
		self.sess1 = intup[-1]
		intup2 = intup[:-1]
		# Get all PyKlein2D attributes
		PyKlein2D.__init__(self, intup2)

	def solve(self):
		tstates = self.tstates
		sess1 = self.sess1
		ssigma = self.ssigma
		mlabraw.eval(sess1,'clear all;')
		mlabraw.eval(sess1,'cd '+mlabpath)
		mlabraw.eval(sess1,'cd Klein')
		mlabraw.put(sess1,'jac',self.gra)
		mlabraw.put(sess1,'hess',self.hes)
		mlabraw.put(sess1,'nstates',self.tstates)
		mlabraw.put(sess1,'ssigma',self.ssigma)
		mlabraw.eval(sess1,'[ff,pp,ee,gg,kx,ky] = solab2(jac,hess,ssigma,nstates)')
		self.F = N.matrix(mlabraw.get(sess1,'ff'))
		self.P = N.matrix(mlabraw.get(sess1,'pp'))
		self.E = N.matrix(mlabraw.get(sess1,'ee'))
		self.G = N.matrix(mlabraw.get(sess1,'gg'))
		self.KX = N.matrix(mlabraw.get(sess1,'kx'))
		self.KY = N.matrix(mlabraw.get(sess1,'ky'))
#----------------------------------------------------------------------------------------------------------------------
class ForKleinD(PyKlein2D):

	def __init__(self,intup):
		self.gra = intup[0]
		self.nendo = intup[1]
		self.nexo = intup[2]
		self.ncon = intup[3]
		self.sigma = intup[4]
		self.A = intup[5]
		self.B = intup[6]
		self.vardic = intup[7]
		self.vdic = intup[8]
		self.modname = intup[9]
		self.audic = intup[10]
		if self.vardic['other']['var']:
			self.numjl = intup[11]
			self.nother = intup[12]
			self.oswitch = 1
		else:
			self.oswitch = 0
		self.tstates = self.nendo+self.nexo
		self.ssigma = self.mkssigma(self.tstates,self.nexo,self.sigma)

	def solve(self):
		A = self.A
		B = self.B
		tstates = self.tstates
		(F,P,retcon) = isolab(A,B,tstates,MAT.shape(A)[0])
		if MAT.sum(P.reshape(-1,1)) == 0.0:
			return
		else:
			self.P = N.matrix(P)
			self.F = N.matrix(F)

	def sim(self,tlen,sntup=None):
		# Add 1000 observations, as they will be thrown away
		# Add one more observation to start first-order vector
		exoli = [x[1] for x in self.vardic['exo']['var']]
		indx=[]
		if sntup == None:
			indx = range(len(exoli))
		else:
			for name in sntup:
				if name not in exoli:
					return 'Error: '+name+' is not a valid exo variable for this model!'
				else:
					indx.append(exoli.index(name))
		indx.sort()
		ncon = self.ncon
		nexo = self.nexo
		nendo = self.nendo
		tstates = self.tstates
		tlena = 1000+tlen
		sigma = self.sigma
		ssigma = self.ssigma
		if self.oswitch:
			numjl = self.numjl

		pp = self.P
		ff = self.F

		count = 0
		for varia in MAT.diag(sigma):
			if locals().has_key('ranvec'):
				if count in indx:
					ranvec = MAT.vstack((ranvec,N.sqrt(varia)*MAT.matrix(N.random.standard_normal(tlena))))
				else:
					ranvec = MAT.vstack((ranvec,MAT.zeros((1,tlena))))		
			else:
				if count in indx:
					ranvec = N.sqrt(varia)*MAT.matrix(N.random.standard_normal(tlena))
				else:
					ranvec = MAT.zeros((1,tlena))
			count = count + 1

		ranvec = MAT.vstack((ranvec,MAT.zeros((nendo,tlena))))

		x_one_m1 = ranvec[:,0]
		y_one_0 = ff*x_one_m1
		y_one_m1 = MAT.zeros(y_one_0.shape)
		x_one_0 = pp*x_one_m1
		if self.oswitch:
			nother = self.nother
			o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

		x_one_c = COP.deepcopy(x_one_m1)
		y_one_c = COP.deepcopy(y_one_0)
		if self.oswitch:
			o_one_c = COP.deepcopy(o_one_0)

		x_one = COP.deepcopy(x_one_c)
		y_one = COP.deepcopy(y_one_c)
		if self.oswitch:
			o_one = COP.deepcopy(o_one_c)

		for i1 in range(1,ranvec.shape[1],1):

			x_one_n = pp*x_one_c+ranvec[:,i1]
			y_one_n = ff*x_one_c
			if self.oswitch:
				o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

			x_one = MAT.hstack((x_one,x_one_c))
			y_one = MAT.hstack((y_one,y_one_n))
			if self.oswitch:
				o_one = MAT.hstack((o_one,o_one_n))

			x_one_c = x_one_n
			y_one_c = y_one_n

		# Throw away first 1000 obs
		x_one = x_one[:,1000:]
		y_one = y_one[:,1000:]
		if self.oswitch:
			o_one = o_one[:,1000:]

		self.sim_x_one = x_one
		self.sim_y_one = y_one
		self.insim = [self.sim_x_one,self.sim_y_one]
		if self.oswitch:
			self.sim_o_one = o_one
			self.insim = self.insim + [self.sim_o_one,]

	def irf(self,tlen,sntup):
		tlen = tlen + 1
		ncon = self.ncon
		nexo = self.nexo
		nendo = self.nendo
		tstates = self.tstates
		sigma = self.sigma
		ssigma = self.ssigma
		tlen = tlen
		if self.oswitch:
			numjl = self.numjl

		pp = self.P
		ff = self.F


		sposli=[]
		exoli = [x[1] for x in self.vardic['exo']['var']]
		# Check if names are valid
		for name in sntup:
			if name not in exoli:
				return 'Error: '+name+' is not a valid exoshock name for this model!'
		for name in sntup:
			sposli.append(exoli.index(name))
		# Expose spos choice to self for show_irf
		self.spos = (sntup,sposli)

		shock = MAT.zeros((nexo,1))
		sendo = MAT.zeros((nendo,1))

		for spos in sposli:
			shock[spos,0] = 1.0
		shock = MAT.vstack((shock,sendo))

		x_one_m1 = shock
		y_one_0 = ff*x_one_m1
		y_one_m1 = MAT.zeros(y_one_0.shape)
		x_one_0 = pp*x_one_m1
		if self.oswitch:
			nother = self.nother
			o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

		x_one_c = COP.deepcopy(x_one_m1)
		y_one_c = COP.deepcopy(y_one_0)
		if self.oswitch:
			o_one_c = COP.deepcopy(o_one_0)

		x_one = COP.deepcopy(x_one_c)
		y_one = COP.deepcopy(y_one_c)
		if self.oswitch:
			o_one = COP.deepcopy(o_one_c)

		for i1 in range(tlen):
			x_one_n = pp*x_one_c
			y_one_n = ff*x_one_c
			if self.oswitch:
				o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

			x_one = MAT.hstack((x_one,x_one_c))
			y_one = MAT.hstack((y_one,y_one_n))
			if self.oswitch:
				o_one = MAT.hstack((o_one,o_one_n))

			x_one_c = x_one_n
			y_one_c = y_one_n

		# Throw away first observation
		self.irf_x_one = x_one[:,1:-1]
		self.irf_y_one = y_one[:,1:-1]
		self.inirf = [self.irf_x_one,self.irf_y_one]
		if self.oswitch:
			self.irf_o_one = o_one[:,2:]
			self.inirf = self.inirf + [self.irf_o_one,]
#------------------------currently not called and really not working, removed from init----------------
class FairTaylor:

	def __init__(self,algtype='type1',intup=None):
		self.algtype = algtype
		self.data = intup[0]
		self.param = intup[1]
		self.sstates_list = intup[2]
		self.vardic = intup[3]
		self.vardic2 = intup[4]
		self.modeq = intup[5]
		self.maxlags = intup[6]
		self.maxleads = intup[7]
		self.maxshocks = intup[8]
		self.maxinfosets = intup[9]
		self.shock_vars = intup[10]
		self.iid_vars = intup[11]
		self.vars = intup[12]
		self.re_var = intup[13]
		self.re_var2 = intup[14]
		self.sstate = intup[15]

		self.prepdat(self.data)
		self.prepmodeq(self.modeq)

	def prepdat(self,data=None):
		self.cur_pres = {}
		self.cur_past = {}
		self.cur_fut = {}
		for x1 in data.items():
			self.cur_pres[x1[0]] = x1[1][self.maxlags:len(x1[1])-self.maxleads]
			self.cur_past[x1[0]] = x1[1][0:-2]
			self.cur_fut[x1[0]] = x1[1][2:]

	def prepmodeq(self,modeq=None):
		re_var2 = self.re_var2
		self.modeq1 = COP.deepcopy(modeq)
		tmp_list = []
		i1 = 0
		for x1 in modeq:
			itera = re_var2.finditer(x1)
			for x2 in itera:
				if x2.group('svar'):
					tmp_list = tmp_list + [[i1,x2.group('svar'),x2.start(),x2.end()],]
			i1 = i1 + 1
		tmp_list.reverse()
		if tmp_list:
			for x2 in tmp_list:
				i1 = x2[0]
				expr = x2[1]
				estart = x2[2]
				eend = x2[3]
				x1 = self.modeq1[i1][0:estart]+expr.capitalize()+\
				   '_bar'+self.modeq1[i1][eend:]
				self.modeq1[i1] = x1
		self.modeq2 = self.modeq1[0:len(modeq)-self.maxshocks]


		self.modeq3 = COP.deepcopy(self.modeq2)
		# Create varlist from vardic to establish fixed order
		self.varlist = self.vardic.items()
		# Create ordered vars and varnames
		varnames = []
		variables = []
		sub_list1 = []
		sub_list2 = []
		tmp_index = {}
		for x1 in self.varlist:
			varnames = varnames + [x1[1],]
			variables = variables + [x1[0],]
		self.varnames = varnames
		self.variables = variables
		re_var2 = self.re_var2
		i1 = 0
		i2 = 0
		for x1 in self.modeq2:
			itera = re_var2.finditer(x1)
			for x2 in itera:
				if x2.group('vvar') and not x2.group('exp') and not x2.group('lagt'):
					tmp_match = variables.index(x2.group('var'))
					if not tmp_index.has_key(x2.group('var')):
						tmp_index[x2.group('var')] = (i2,x2.group(0))
						i2 = i2 + 1
					sub_list1 = sub_list1 + [['invar['+str(tmp_index[x2.group('var')][0])+']',
								  x2.start(),x2.end(),x2.group(0),tmp_match,i1],]

			sub_list1.reverse()
			if sub_list1:
				for x4 in sub_list1:
					str_tmp = x4[0]
					pstart = x4[1]
					pend = x4[2]
					self.modeq3[i1] = self.modeq3[i1][0:pstart]+str_tmp+self.modeq3[i1][pend:]
				sub_list2 = sub_list2 + [sub_list1,]
				sub_list1 = []
			i1 = i1 + 1
		self.index = tmp_index
		# Create inverted index2
		index2 = {}
		for x1 in self.index.items():
			index2[x1[1][0]] = (x1[0],x1[1][1])
		self.index2 = index2
		return sub_list2

	def makecur(self,tpos,curlag,curpres,curfut):
		lag_dic={}
		pres_dic={}
		fut_dic={}
		all_dic={}
		for x1 in curlag.items():
			if float(x1[1][tpos].data) > 0:
				lag_dic[self.vardic2[x1[0]]+'_L'] = float(x1[1][tpos].data)
			else:
				lag_dic[self.vardic2[x1[0]]+'_L'] = 0.01
		for x1 in curpres.items():
			if float(x1[1][tpos].data) > 0:
				pres_dic[self.vardic2[x1[0]]+'_N'] = float(x1[1][tpos].data)
			else:
				pres_dic[self.vardic2[x1[0]]+'_N'] = 0.01
		for x1 in curfut.items():
			if float(x1[1][tpos].data) > 0:
				fut_dic['INFO_N_'+self.vardic2[x1[0]]+'_F'] = float(x1[1][tpos].data)
			else:
				fut_dic['INFO_N_'+self.vardic2[x1[0]]+'_F'] = 0.01
		all_dic.update(lag_dic)
		all_dic.update(pres_dic)
		all_dic.update(fut_dic)
		self.alldic =all_dic
		return all_dic

	def solveone(self,tpos,curlag,curpres,curfut):

		#Define the function to be handed over
		#to fsolve

		def func(invar):
			locals().update(self.param)
			locals().update(self.sstate[0])
			locals().update(self.makecur(tpos,self.cur_past,self.cur_pres,self.cur_fut))
			fdot1 = S.zeros((len(self.modeq3)),'d')
			i1=0
			for x in self.modeq3:
				fdot1[i1] = eval(x)
				i1 = i1 + 1
			return fdot1


		#def func(invar):
			#locals().update(self.param)
			#locals().update(self.sstate[0])
			#locals().update(self.makecur(tpos,self.cur_past,self.cur_pres,self.cur_fut))
			#fdot1 = S.zeros((len(self.modeq3)),'d')
			#i1=0
			#for x in self.modeq3:
				#fdot1[i1] = abs(eval(x))
				#i1 = i1 + 1
			#return sum(fdot1)

		#def fprime(invar):
			#return N.array(([1,-1]))


		# Define the initial values and
		# start the non-linear solver
		init_val = len(self.index)*[0,]
		for x1 in self.index.items():
			varia = self.vardic[x1[0]]['var']
			init_val[x1[1][0]] = float(curpres[varia][tpos].data)

		# Determine bounds
		#in_bounds = len(self.index)*[(0.01,None),]

		(output,infodict,ier,mesg) = O.fsolve(func,init_val,full_output=1)

		#(output,f_out,d_out) = O.fmin_l_bfgs_b(func,init_val,
															#fprime,approx_grad=True,bounds=in_bounds)


		# Attach the outputs of the solver as attributes
		fsout={}
		self.output = output
		i1 = 0
		for x1 in self.index2.items():
			fsout[x1[1][1]] = (output[x1[0]],x1[1][0])
			i1 = i1 + 1
		self.fsout = fsout
		self.output = output
		self.infodict = infodict
		self.ier = ier
		self.mesg = mesg

		return fsout

	def solveall(self):
		# Calculate initial convergence criterion
		criterion = 10
		loop_c = 0
		while criterion > 0.01:
			cur_past_t = {}
			cur_pres_t = {}
			cur_fut_t = {}
			for x1 in self.vardic2.items():
				cur_past_t[x1[0]] = N.zeros(len(self.cur_pres.items()[0][1]))
				cur_past_t[x1[0]][0] = float(self.cur_past[x1[0]][0].data)
				cur_pres_t[x1[0]] = N.zeros(len(self.cur_pres.items()[0][1]))
				cur_fut_t[x1[0]] = N.zeros(len(self.cur_pres.items()[0][1]))
				cur_fut_t[x1[0]][-1] = float(self.cur_fut[x1[0]][-1].data)
			for i1 in range(0,len(self.cur_pres.items()[0][1])):
				out_tmp = self.solveone(i1,self.cur_past,self.cur_pres,self.cur_fut)
				for x2 in out_tmp.items():
					# Standard
					if self.algtype == 'type1':
						cur_pres_t[self.vardic[x2[1][1]['var']]][i1] = x2[1][0]
					# Fast-Gauss Seidel
					elif self.algtype == 'type2':
						cur_pres_t[self.vardic[x2[1][1]['var']]][i1] = x2[1][0]
						try:
							self.cur_past[self.vardic[x2[1][1]['var']]][i1+1] = x2[1][0]
						except:
							pass

			# Recalculate convergence criterion
			cr_list = [0,]*len(self.vardic)
			va_list = [0,]*len(self.vardic)
			for x1,x2 in zip(self.vardic.items(),range(0,len(self.vardic),1)):
				cr_list[x2] = cr_list[x2] + sum(abs(cur_fut_t[x1[1]][0:-1]-cur_pres_t[x1[1]][1:]))
				va_list[x2] = x1[1]
			criterion = sum(cr_list)
			loop_c = loop_c + 1
			print ('Loop: ',loop_c)
			print ('Var:  ',va_list)
			print ('Crit: ',cr_list)


			for x1 in self.vardic2.items():
				cur_past_t[x1[0]][1:] = cur_pres_t[x1[0]][:-1]
				cur_fut_t[x1[0]][0:-1] = cur_pres_t[x1[0]][1:]

			for x1,x2,x3 in zip(cur_past_t.items(),cur_pres_t.items(),cur_fut_t.items()):
				self.cur_past[x1[0]] = TSS.time_series\
				    (x1[1],start_date=TSS.Date(freq=self.cur_past[x1[0]].start_date.freqstr,
							       year=self.cur_past[x1[0]].start_date.year,
							       quarter=self.cur_past[x1[0]].start_date.quarter),
							       freq=self.cur_past[x1[0]].start_date.freqstr)

				self.cur_pres[x1[0]] = TSS.time_series\
				    (x2[1],start_date=TSS.Date(freq=self.cur_pres[x1[0]].start_date.freqstr,
							       year=self.cur_pres[x1[0]].start_date.year,
							       quarter=self.cur_pres[x1[0]].start_date.quarter),
							       freq=self.cur_pres[x1[0]].start_date.freqstr)

				self.cur_fut[x1[0]] = TSS.time_series\
				    (x3[1],start_date=TSS.Date(freq=self.cur_fut[x1[0]].start_date.freqstr,
							       year=self.cur_fut[x1[0]].start_date.year,
							       quarter=self.cur_fut[x1[0]].start_date.quarter),
							       freq=self.cur_fut[x1[0]].start_date.freqstr)
"""***********************************************************"""
##########THE STEADY STATE SOLVER CLASS AND IT'S SUBCLASSES (WORKS)#########
class SSsolvers:
	def __init__(self):
		pass
#----------------------------------------------------------------------------------------------------------------------
class Manss():
	def __init__(self,intup):
		self.manss_sys = intup[0]
		self.paramdic = intup[1]

	def solve(self):
		list_tmp1 = COP.deepcopy(self.manss_sys)
		# Create manual (closed-form) steady state dictionary
		_rlog = 'LOG\('
		_rexp = 'EXP\('
		_fsolve = 'ROOT\((?P<nlexp>.*),\s*(?P<vari>.*)\s*=\s*(?P<init>.*)\s*,\s*fail\s*=\s*(?P<fval>.*)\);'
		rlog = RE.compile(_rlog)
		rexp = RE.compile(_rexp)
		fexp = RE.compile(_fsolve)
		manss={}
		locals().update(self.paramdic)
		globals().update(self.paramdic)
		for x in list_tmp1:
			str_tmp3 = x[:]
			while rlog.search(str_tmp3):
				str_tmp3 = RE.sub(rlog,'N.log(',str_tmp3)
			while rexp.search(str_tmp3):
				str_tmp3 = RE.sub(rexp,'N.exp(',str_tmp3)			

			# Do fsolve root finding, if ROOT detected
			if fexp.search(str_tmp3):
				asvari = str_tmp3.split('=')[0].strip()
				ma = fexp.search(str_tmp3)
				nlexp = ma.group('nlexp')
				vari = ma.group('vari')
				init = N.float(ma.group('init'))
				fval = N.float(ma.group('init'))
				solu,infodict,ier,mesg = O.fsolve(eval('lambda '+vari+':'+nlexp),init,full_output=1)
				if ier == 1:
					manss[asvari] = solu
					locals()[asvari] = solu
					globals().update(manss)
				else:
					manss[asvari] = fval
					locals()[asvari] = fval
					globals().update(manss)
				continue
			str_tmp = str_tmp3.split(';')[0]
			list_tmp = str_tmp.split('=')
			str_tmp1 = list_tmp[0].strip()
			str_tmp2 = list_tmp[1].strip()
			manss[str_tmp1] = eval(str_tmp2)
			locals()[str_tmp1] = eval(str_tmp2)
			globals().update(manss)
		self.sstate = manss
#----------------------------------------------------------------------------------------------------------------------
class Fsolve:
	def __init__(self,intup):
		self.ssm = intup[0]
		self.ssi = intup[1]
		self.paramdic = intup[2]

	def solve(self):

		# Turn the non-linear system into a representation
		# that is suitable for fsolve(invar) !

		ssi = self.ssi
		subdic = {}
		for y,z in zip(ssi.items(),range(len(ssi.items()))):
			subdic[y[0]] = (y[0],y[1],z)
		list_tmp1 = COP.deepcopy(self.ssm)
		for var in self.ssi.keys():
			_mreg = '(\+|\*|-|/|^|[(])'+var
			mreg = RE.compile(_mreg)
			for i1,line in enumerate(list_tmp1):
				while mreg.search(list_tmp1[i1]):
					ma = mreg.search(list_tmp1[i1])
					if len(ma.group(1))>0:
						var = ma.group()[1:]
						pos = ma.span()[0]+1
					else:
						var = ma.group()
						pos = ma.span()[0]						
					poe = ma.span()[1]
					list_tmp1[i1] = list_tmp1[i1][:pos]+'invar['+str(subdic[var][2])+']'+list_tmp1[i1][poe:]

		func_repr = list_tmp1

		# Define the function to be handed over
		# to fsolve
		def func(invar):
			locals().update(self.paramdic)
			fdot = S.zeros((len(func_repr)),'d')
			i1=0
			for x in func_repr:
				fdot[i1] = eval(x)
				i1 = i1 + 1
			return fdot

		# Define the initial values and
		# start the non-linear solver
		inlist = []
		for x in subdic.values():
			inlist.append([x[2],x[1],x[0]])
		inlist.sort()
		init_val = [x[1] for x in inlist]
		(output,infodict,ier,mesg) = O.fsolve(func,init_val,full_output=1)
		# Attach the outputs of the solver as attributes
		self.fsout={}
		for x,y in zip(output,inlist):
			self.fsout[y[2]] = x
		self.infodict = infodict
		self.ier = ier
		self.mesg = mesg
"""***********************************************************"""
###########THE MODEL EVALUATION CLASS AND IT'S SUBLCASSES (INCOMPLETE)#######
class MODeval:

	def __init__(self):
		self.test = 'tester'
#----------------------------------------------------------------------------------------------------------------------
class Minford:

	def __init__(self):
		self.test = 'tester'
"""***********************************************************"""
######################ERROR (EXCEPTION) CLASSES####################
#class MyErr(mlabraw.error):
#	pass
class MyErr(ValueError):
    pass
"""***********************************************************"""
###########################GUI CLASSES#########################
class SimpleGrid(gridlib.Grid): 
	def __init__(self, parent):
		gridlib.Grid.__init__(self, parent, -1)
		##mixins.GridAutoEditMixin.__init__(self)
		self.moveTo = None
		self.Bind(wx.EVT_IDLE, self.OnIdle)
		self.CreateGrid(25, 25)#, gridlib.Grid.SelectRows)
		##self.EnableEditing(False)
		# simple cell formatting
		self.SetColSize(3, 200)
		self.SetRowSize(4, 45)
		self.SetCellValue(0, 0, "First cell")
		self.SetCellValue(1, 1, "Another cell")
		self.SetCellValue(2, 2, "Yet another cell")
		self.SetCellValue(3, 3, "This cell is read-only")
		self.SetCellFont(0, 0, wx.Font(12, wx.ROMAN, wx.ITALIC, wx.NORMAL))
		self.SetCellTextColour(1, 1, wx.RED)
		self.SetCellBackgroundColour(2, 2, wx.CYAN)
		self.SetReadOnly(3, 3, True)
		self.SetCellEditor(5, 0, gridlib.GridCellNumberEditor(1,1000))
		self.SetCellValue(5, 0, "123")
		self.SetCellEditor(6, 0, gridlib.GridCellFloatEditor())
		self.SetCellValue(6, 0, "123.34")
		self.SetCellEditor(7, 0, gridlib.GridCellNumberEditor())
		self.SetCellValue(6, 3, "You can veto editing this cell")
		#self.SetRowLabelSize(0)
		#self.SetColLabelSize(0)
		# attribute objects let you keep a set of formatting values
		# in one spot, and reuse them if needed
		attr = gridlib.GridCellAttr()
		attr.SetTextColour(wx.BLACK)
		attr.SetBackgroundColour(wx.RED)
		attr.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
		# you can set cell attributes for the whole row (or column)
		self.SetRowAttr(5, attr)
		self.SetColLabelValue(0, "Custom")
		self.SetColLabelValue(1, "column")
		self.SetColLabelValue(2, "labels")
		self.SetColLabelAlignment(wx.ALIGN_LEFT, wx.ALIGN_BOTTOM)
		#self.SetDefaultCellOverflow(False)
		#r = gridlib.GridCellAutoWrapStringRenderer()
		#self.SetCellRenderer(9, 1, r)
		# overflow cells
		self.SetCellValue( 9, 1, "This default cell will overflow into neighboring cells, but not if you turn overflow off.");
		self.SetCellSize(11, 1, 3, 3);
		self.SetCellAlignment(11, 1, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE);
		self.SetCellValue(11, 1, "This cell is set to span 3 rows and 3 columns");
		editor = gridlib.GridCellTextEditor()
		editor.SetParameters('10')
		self.SetCellEditor(0, 4, editor)
		self.SetCellValue(0, 4, "Limited text")
		renderer = gridlib.GridCellAutoWrapStringRenderer()
		self.SetCellRenderer(15,0, renderer)
		self.SetCellValue(15,0, "The text in this cell will be rendered with word-wrapping")

		# test all the events
		self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.OnCellLeftClick)
		self.Bind(gridlib.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClick)
		self.Bind(gridlib.EVT_GRID_CELL_LEFT_DCLICK, self.OnCellLeftDClick)
		self.Bind(gridlib.EVT_GRID_CELL_RIGHT_DCLICK, self.OnCellRightDClick)
		self.Bind(gridlib.EVT_GRID_LABEL_LEFT_CLICK, self.OnLabelLeftClick)
		self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_CLICK, self.OnLabelRightClick)
		self.Bind(gridlib.EVT_GRID_LABEL_LEFT_DCLICK, self.OnLabelLeftDClick)
		self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_DCLICK, self.OnLabelRightDClick)
		self.Bind(gridlib.EVT_GRID_ROW_SIZE, self.OnRowSize)
		self.Bind(gridlib.EVT_GRID_COL_SIZE, self.OnColSize)
		self.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.OnRangeSelect)
		self.Bind(gridlib.EVT_GRID_CELL_CHANGE, self.OnCellChange)
		self.Bind(gridlib.EVT_GRID_SELECT_CELL, self.OnSelectCell)
		self.Bind(gridlib.EVT_GRID_EDITOR_SHOWN, self.OnEditorShown)
		self.Bind(gridlib.EVT_GRID_EDITOR_HIDDEN, self.OnEditorHidden)
		self.Bind(gridlib.EVT_GRID_EDITOR_CREATED, self.OnEditorCreated)
	def OnCellLeftClick(self, evt):
		evt.Skip()
	def OnCellRightClick(self, evt):
		evt.Skip()
	def OnCellLeftDClick(self, evt):
		evt.Skip()
	def OnCellRightDClick(self, evt):
		evt.Skip()
	def OnLabelLeftClick(self, evt):
		evt.Skip()
	def OnLabelRightClick(self, evt):
		evt.Skip()
	def OnLabelLeftDClick(self, evt):
		evt.Skip()
	def OnLabelRightDClick(self, evt):
		evt.Skip()
	def OnRowSize(self, evt):
		evt.Skip()
	def OnColSize(self, evt):
		evt.Skip()
	def OnRangeSelect(self, evt):
		if evt.Selecting():
			pass
		evt.Skip()
	def OnCellChange(self, evt):
		# Show how to stay in a cell that has bad data.  We can't just
		# call SetGridCursor here since we are nested inside one so it
		# won't have any effect.  Instead, set coordinates to move to in
		# idle time.
		value = self.GetCellValue(evt.GetRow(), evt.GetCol())
		if value == 'no good':
			self.moveTo = evt.GetRow(), evt.GetCol()
	def OnIdle(self, evt):
		if self.moveTo != None:
			self.SetGridCursor(self.moveTo[0], self.moveTo[1])
			self.moveTo = None
		evt.Skip()
	def OnSelectCell(self, evt):
		# Another way to stay in a cell that has a bad value...
		row = self.GetGridCursorRow()
		col = self.GetGridCursorCol()
		if self.IsCellEditControlEnabled():
			self.HideCellEditControl()
			self.DisableCellEditControl()
		value = self.GetCellValue(row, col)
		if value == 'no good 2':
			return  # cancels the cell selection
		evt.Skip()
	def OnEditorShown(self, evt):
		if evt.GetRow() == 6 and evt.GetCol() == 3 and \
		   wx.MessageBox("Are you sure you wish to edit this cell?",
				 "Checking", wx.YES_NO) == wx.NO:
			evt.Veto()
			return
		evt.Skip()
	def OnEditorHidden(self, evt):
		if evt.GetRow() == 6 and evt.GetCol() == 3 and \
		   wx.MessageBox("Are you sure you wish to  finish editing this cell?",
				 "Checking", wx.YES_NO) == wx.NO:
			evt.Veto()
			return
		evt.Skip()
	def OnEditorCreated(self, evt):
		pass
#----------------------------------------------------------------------------------------------------------------------
class TestFrame(wx.Frame):
	def __init__(self, parent):
		wx.Frame.__init__(self, parent, -1, "Simple Grid Demo", size=(640,480))
		self.grid = SimpleGrid(self)
##############COLLECTION OF SUPPORTING FUNCTIONS AND CLASSES#############
def locate(stringlines,varlist):
	'''
	A function, takes stringlines (split file) and a varlist
	and then creates a location dictionary for location of strings in file
	'''
	locdic = {}
	for x in varlist:
		row_iter=0
		while row_iter < len(stringlines):
			if x[0] in stringlines[row_iter]:
				locdic[x[1]] = row_iter
				break
			else:
				row_iter=row_iter+1
	return locdic

class dicwrap:
	def __init__(self,other,initlev,nreg):
		self.other = other
		self.initlev = initlev
		self.nreg = nreg
		self.wrapdic = other.paramdic
	def __getattr__(self,attrname):
		return getattr(self.wrapdic,attrname)

	def __setitem__(self,key,value):
		other = self.other
		initlev = self.initlev
		nreg = self.nreg
		if self.wrapdic[key] != value:
			self.wrapdic[key] = value
			################## STEADY STATE CALCULATIONS !!! ##############
			# Delete the models previous results
			del other.modsolvers
			del other.sssolvers
			other.sssolvers = SSsolvers()
			# Reset switch back to zero
			other.switches['ss_suc'] = ['0','0']
			# Solve for steady-state using fsolve
			if sum([nreg.search(x)!=None for x in other.txtpars.secs['ssm'][0]]) == 0:
				intup = (other.ssys_list,other.ssidic,other.paramdic)
				other.sssolvers.fsolve = Fsolve(intup)
				other.sssolvers.fsolve.solve()
				if other.sssolvers.fsolve.ier == 1:
					nsstate = other.sssolvers.fsolve.fsout
					mreg = RE.compile('(.*?)_bar')
					mreg2 = RE.compile('(.*?)_F|B.*?_bar')
					for key in nsstate.keys():
						ma = mreg.search(key)
						if ma:
							vari = ma.group(1)
							for okey in other.sstate:
								ma2 = mreg2.search(okey)
								if ma2:
									vari2 = ma2.group(1)
									if vari == vari2:
										other.sstate[okey] = nsstate[key]
					other.sstate.update(nsstate)
					other.numssdic.update(nsstate)
					other.switches['ss_suc'] = ['1','1']
				else:
					other.switches['ss_suc'] = ['1','0']
			# Solve for steady-state using manss
			if sum([nreg.search(x)!=None for x in other.txtpars.secs['sss'][0]]) == 0:
				if other.switches['ss_suc'] == ['1','1']:
					alldic = {}
					alldic.update(other.sstate)
					alldic.update(other.paramdic)
					intup = (other.manss_sys,alldic)
					other.sssolvers.manss = Manss(intup)
					other.sssolvers.manss.solve()
					nsstate = other.sssolvers.manss.sstate
					mreg = RE.compile('(.*?)_bar')
					mreg2 = RE.compile('(.*?)_F|B.*?_bar')
					for key in nsstate.keys():
						ma = mreg.search(key)
						if ma:
							vari = ma.group(1)
							for okey in other.sstate:
								ma2 = mreg2.search(okey)
								if ma2:
									vari2 = ma2.group(1)
									if vari == vari2:
										other.sstate[okey] = nsstate[key]
					other.sstate.update(nsstate)
				else:
					intup = (other.manss_sys,other.paramdic)
					other.sssolvers.manss = Manss(intup)
					other.sssolvers.manss.solve()
					nsstate = other.sssolvers.manss.sstate
					mreg = RE.compile('(.*?)_bar')
					mreg2 = RE.compile('(.*?)_F|B.*?_bar')
					for key in nsstate.keys():
						ma = mreg.search(key)
						if ma:
							vari = ma.group(1)
							for okey in other.sstate:
								ma2 = mreg2.search(okey)
								if ma2:
									vari2 = ma2.group(1)
									if vari == vari2:
										other.sstate[okey] = nsstate[key]
					other.sstate.update(nsstate)
			if initlev == 1: return
	
			# Open the model solution tree branch
			other.modsolvers = MODsolvers()
######################## LINEAR METHODS !!! ############################
			if sum([nreg.search(x)!=None for x in other.txtpars.secs['modeq'][0]]) == 0:
				# Open the matlab Uhlig object
				intup = ((other.nendo,other.ncon,other.nexo),
					 other.eqindx,
					 other.vreg,
					 other.llsys_list,
					 other.diffli1,
					 other.diffli2,
					 sess1,
					 other.vardic)
				other.modsolvers.matuhlig = MatUhlig(intup)
				# Open the native Uhlig object
				intup = ((other.nendo,other.ncon,other.nexo),
					 other.eqindx,
					 other.vreg,
					 other.llsys_list,
					 other.diffli1,
					 other.diffli2,
					 sess1)
				other.modsolvers.pyuhlig = PyUhlig(intup)
				# Open the matlab Klein object
				intup = ((other.nendo,other.ncon,other.nexo),
					 other.eqindx,
					 other.vreg,
					 other.llsys_list,
					 other.diffli1,
					 other.diffli2,
					 sess1)
				other.modsolvers.matklein = MatKlein(intup)
				# Open the Fortran Klein object
				intup = ((other.nendo,other.ncon,other.nexo),
					 other.eqindx,
					 other.vreg,
					 other.llsys_list,
					 other.diffli1,
					 other.diffli2,
					 sess1)
				other.modsolvers.forklein = ForKlein(intup)
################## 1ST-ORDER NON-LINEAR METHODS !!! ##################
			if sum([nreg.search(x)!=None for x in other.txtpars.secs['focs'][0]]) == 0:
	
				# First, create the Jacobian and (possibly-->mk_hessian==True?) Hessian
				if use_anaderiv:
					if ncpus > 1 and mk_hessian:
						other.mkjaheppn()
					elif ncpus > 1 and not mk_hessian:
						other.mkjaheppn()
					else:
						other.mkjahen()
				else:
					other.mkjahenmat()
	
				# Open the MatWood object
				intup = (other.jAA,other.jBB,
					 other.nexo,other.ncon,
					 other.nendo,sess1)
				other.modsolvers.matwood = MatWood(intup)
				# Open the Fortran KleinD object
				if 'nlsubsys' in dir(other):
					intup = (other.numj,
						 other.nendo,other.nexo,
						 other.ncon,other.sigma,
						 other.jAA,other.jBB,
						 other.vardic,other.vdic,
						 other.modname,other.audic,
						 other.numjl,
						 other.nother)
				else:
					intup = (other.numj,
						 other.nendo,other.nexo,
						 other.ncon,other.sigma,
						 other.jAA,other.jBB,
						 other.vardic,other.vdic,
						 other.modname,other.audic)
				other.modsolvers.forkleind = ForKleinD(intup)
################## 2ND-ORDER NON-LINEAR METHODS !!! ##################
				if sum([nreg.search(x)!=None for x in other.txtpars.secs['vcvm'][0]]) == 0 and\
				   'numh' in dir(other):
					# Open the MatKlein2D object
					if 'nlsubsys' in dir(other):
						intup = (other.numj,other.numh,
							 other.nendo,other.nexo,
							 other.ncon,other.sigma,
							 other.jAA,other.jBB,
							 other.vardic,other.vdic,
							 other.modname,other.audic,
							 other.numjl,other.numhl,
							 other.nother,sess1)
					else:
						intup = (other.numj,other.numh,
							 other.nendo,other.nexo,
							 other.ncon,other.sigma,
							 other.jAA,other.jBB,
							 other.vardic,other.vdic,
							 other.modname,other.audic,
							 sess1)
					other.modsolvers.matklein2d = MatKlein2D(intup)
					# Open the PyKlein2D object
					intup = intup[:-1]
					other.modsolvers.pyklein2d = PyKlein2D(intup)
	
			if 'jAA' in dir(other):
				other.mkeigv()

	def __getitem__(self,key):
		return self.wrapdic[key]
	def __repr__(self):
		return self.wrapdic.__repr__()
	def __str__(self):
		return self.wrapdic.__str__()
