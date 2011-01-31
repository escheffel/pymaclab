#from __future__ import division
import wx
import wx.grid as gridlib
import trace
import time
import datetime
from pymaclab.filters._hpfilter import hpfilt
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

from ..stats.var import VAR #TODO: remove for statsmodels version
from solvers.steadystate import SSsolvers, Manss, Fsolve
from solvers.modsolvers import (MODsolvers, PyUhlig, MatUhlig, MatKlein,
        MatKleinD, MatWood, ForKlein, PyKlein2D, MatKlein2D, ForKleinD,
        FairTaylor)

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
try:    # no mp in 2.5
    from multiprocessing import cpu_count
    ncpus = cpu_count()
except:
    ncpus = 1
ncpus = 1 # debug without parallel stuff
#NOTE: also see when this would actually speed things up.  doesn't seem to 
# now

# Update level
updlev = 0

# Empty locdic for the locate helper function
locdic = {}


# Open mlab session
if use_matlab:
    sess1 = mlabraw.open('matlab - nojvm -nosplash')

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
    def __init__(self,ffile=None,dbase=None,initlev=2):
        """
        DSGE Model Class

        Parameters
        ----------
        ffile : str
            The modfile
        dbase : DataBase object
            Not sure what this is for yet
        initlev : int
            0 just parse the modfile
            1 parse modfile and solve steady state
            2 Jacobian (and Hessian)
        """
        self._initlev = initlev #TODO: remove this as an option
        # Set no author
        self.setauthor()
        # Open the switches dic and initiate
        #NOTE this is _not_ pythonic
        self.switches = {}
        self.switches['errocc'] = '0'
        self.switches['ss_suc'] = ['0','0']
        # Attach some model attributes
        self.modfile = ffile
        self.dbase = dbase

    # Initializes all of the rest, errors can occur here ! (steady state, jacobian, hessian)
    def init2(self):
        initlev = self._initlev

        # Create None tester regular expression
        _nreg = '^\s*None\s*$'
        nreg = RE.compile(_nreg)

        self.txtpars = TXTparser(self.modfile)
        # Start to populate the model from file
        self.pop(input=self.txtpars)
        # Attach the data from database
        if self.dbase != None:
            self.getdata(dbase=self.dbase)
        if initlev == 0: 
            return

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
#NOTE: don't do this?
# or don't modify _existing_ contents?
# see: http://docs.python.org/library/functions.html#locals
#also locals() is no longer a dictionary in Python 3
#NOTE: just use the dictionary itself
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
        from sympycore import FunctionRing
        imports = ('numpy','copy','numpy.matlib',)
#job_server.submit(self, func, args=(), depfuncs=(), modules=(), callback=None, callbackargs=(), group='default', globals=None) 
#    Submits function to the execution queue                            
   
#    func - function to be executed
#    args - tuple with arguments of the 'func'
#    depfuncs - tuple with functions which might be called from 'func'
#    modules - tuple with module names to import
#    callback - callback function which will be called with argument
#            list equal to callbackargs+(result,)
#            as soon as calculation is done
#    callbackargs - additional arguments for callback function
#    group - job group, is used when wait(group) is called to wait for
#    jobs in a given group to finish
#    globals - dictionary from which all modules, functions and classes
#    will be imported, for instance: globals=globals()
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
#TODO: make sure that everything gets evaluated in float terms
#for instance, in old cee.txt 1.03**(1/4) does int division

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
#   pass
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
