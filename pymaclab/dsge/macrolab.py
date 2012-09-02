#from __future__ import division
from scikits import timeseries as ts
from copy import deepcopy
import os
import re
from helpers import now_is
import numpy as np
from numpy import matlib as mat
from scipy.linalg import eig as scipyeig
#from scikits.timeseries.lib import plotlib as tsplt
from tempfile import mkstemp
import trace
import time
import datetime
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
    # Import Sympy
    import sympy
except:
    print "You need to have Sympy installed for PyMacLab to work properly"
    print "Get is via: sudo pip install sympy"
try:
    import pp
except:
    print "You need to install pp (easy_install)"
    print "this is optional"

#NOTE: Imports from the refactor
from pymaclab.filters._hpfilter import hpfilter
from pymaclab.filters._bkfilter import bkfilter
from ..stats.var import VAR #TODO: remove for statsmodels version
from solvers.steadystate import SSsolvers, ManualSteadyState, Fsolve
#TODO: delay above and only import if needed
from solvers.modsolvers import (MODsolvers, PyUhlig, MatUhlig, MatKlein,
        MatKleinD, MatWood, ForKlein, PyKlein2D, MatKlein2D, ForKleinD,
        FairTaylor)
from parsers._modparser import parse_mod
from parsers._dsgeparser import populate_model_stage_one,populate_model_stage_one_b,populate_model_stage_two
from tools import dicwrap

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
    import pp
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
        str_tmp1 = open(os.path.join(datapath,filename),'r')
        flist = str_tmp1.read().splitlines()[:]
        dat_Start = flist[0].split(',')[1].strip()[:]
        dat_End = flist[1].split(',')[1].strip()[:]
        dat_Freq = re.sub('"','',flist[2].split(',')[1]).strip()[:]
        dat_Name = re.sub('"','',flist[3].split(',')[1]).strip()[:]
        dat_Code = re.sub('"','',flist[4].split(',')[1]).strip()[:]     
        dat_Curr = re.sub('"','',flist[5].split(',')[1]).strip()[:]
        datadate = flist[6:]
        datacsv = ''
        i1=0
        for x in datadate:
            str_tmp1 = x.split(',')[1]
            datacsv = datacsv[:] + str_tmp1[:]+'\n'
            i1 = i1 + 1
        output = open(os.getcwd()+'/'+'tmp_001.csv', 'w')
        output.write(datacsv)
#        pdata = P.load('tmp_001.csv',delimiter=',')
        pdata = np.loadtxt('tmp_001.csv', delimiter=",")
        os.remove(os.getcwd()+'/'+'tmp_001.csv')
        if dat_Freq == 'D':
            tsdata = data = ts.time_series(pdata,start_date=ts.Date(freq=dat_Freq,
                                          year=int(dat_Start[6:10]),
                                          month=int(dat_Start[3:5]),
                                          day=int(dat_Start[0:2])))
        elif dat_Freq == 'M':
            tsdata = data = ts.time_series(pdata,start_date=ts.Date(freq=dat_Freq,
                                          year=int(dat_Start[6:10]),
                                          month=int(dat_Start[3:5])))
        elif dat_Freq == 'Q':
            tsdata = data = ts.time_series(pdata,start_date=ts.Date(freq=dat_Freq,
                                          year=int(dat_Start[6:10]),
                                          quarter=int(dat_Start[3:4])))

        self.datdic[dat_Code]={}
        self.datdic[dat_Code]['infile'] = tsdata
        self.datdic[dat_Code]['is_DatStr'] = True
        self.datdic[dat_Code]['Dat_Desc'] = dat_Name
        self.datdic[dat_Code]['DatStr_Code'] = dat_Code
        self.datdic[dat_Code]['imTimeStamp'] = now_is()
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
        str_tmp1 = deepcopy(csvfile)
        flist = str_tmp1.splitlines()[:]
        dat_Start = flist[0].split(',')[1].strip()[:]
        dat_End = flist[1].split(',')[1].strip()[:]
        dat_Freq = re.sub('"','',flist[2].split(',')[1]).strip()[:]
        dat_Name = re.sub('"','',flist[3].split(',')[1]).strip()[:]
        dat_Code = re.sub('"','',flist[4].split(',')[1]).strip()[:]     
        dat_Curr = re.sub('"','',flist[5].split(',')[1]).strip()[:]
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
        output = open(os.getcwd()+'/'+'tmp_001.csv', 'w')
        output.write(datacsv)
#        pdata = P.load('tmp_001.csv',delimiter=',')
        pdate = np.loadtxt('tmp_001.csv', delimiter=",")
        os.remove(os.getcwd()+'/'+'tmp_001.csv')
        if dat_Freq == 'D':
            tsdata = data = time_series(pdata,start_date=ts.Date(freq=dat_Freq,
                                          year=dat_Start[0:1],
                                          month=dat_Start[3:4],
                                          day=dat_Start[6:9]))
        elif dat_Freq == 'M':
            tsdata = data = time_series(pdata,start_date=ts.Date(freq=dat_Freq,
                                          year=dat_Start[0:1],
                                          month=dat_Start[3:4]))
        elif dat_Freq == 'Q':
            tsdata = data = time_series(pdata,start_date=ts.Date(freq=dat_Freq,
                                          year=dat_Start[0:1],
                                          quarter=dat_Start[3:4]))

        self.datdic[dat_Name]={}
        self.datdic[dat_Name]['infile'] = tsdata
        self.datdic[dat_Name]['is_DatStr'] = True
        self.datdic[dat_Name]['DatStr_Code'] = dat_Code
        self.datdic[dat_Name]['imTimeStamp'] = now_is()

    def tsim(self,datname,tsdesc,tsin):
        self.datdic[datname]={}
        self.datdic[datname]['infile'] = tsin
        self.datdic[datname]['is_DatStr'] = False
        self.datdic[datname]['Dat_Desc'] = tsdesc
        self.datdic[datname]['DatStr_Code'] = None
        self.datdic[datname]['imTimeStamp'] = now_is()
        self.datdic[datname]['start_date'] = tsin.start_date
        self.datdic[datname]['end_date'] = tsin.end_date
        self.datdic[datname]['freq'] = tsin.freqstr[0:1]
        self.datdic[datname]['infile'].TSdesc = tsdesc
        self.datdic[datname]['infile'].TSname = datname
        self.datdic[datname]['svalue'] = tsin.dates[0].value
        self.datdic[datname]['evalue'] = tsin.dates[-1].value
        self.datdic[datname]['alias'] = False

    def isDatStreamCSV(self,csvfile):
        str_tmp1 = deepcopy(csvfile)
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
        tsoutf = hpfilter(tsinf,N.zeros((N.shape(tsinf)[0],3)),
            N.shape(tsinf)[0],lam,0)
        tsoutf = TSS.time_series\
               (tsoutf,start_date=tsinf.start_date,freq=tsin['freq'])
        self.tsim(tsout,tsin['Dat_Desc']+',hp-filtered',tsoutf)

    def mklog(self,tsname,tsout):
        tsinf = self.datdic[tsname]['infile']
        tsin = self.datdic[tsname]
        tsoutf = np.log(tsinf)
        self.tsim(tsout,tsin['Dat_Desc']+',logarithmic',tsoutf)

    def mkexp(self,tsname,tsout):
        tsinf = self.datdic[tsname]['infile']
        tsin = self.datdic[tsname]
        tsoutf = np.exp(tsinf)
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
        s_date = ts.Date(freq,value=max_sval)
        e_date = ts.Date(freq,value=min_eval)
        list_tmp = []
        for x1 in self.datdic.items():
            if x1[1]['alias'] != 'None' and x1[1]['freq'] == freq:
                self.modif[x1[1]['alias']] = x1[1]['infile'][s_date:e_date]

    def mkvar(self,varlag=1,varord='None',spos='None',useconst='const'):
        tmp_dic={}
        for x1 in varord:
            tmp_dic[x1[0]] = self.modif[x1[0]]
        dat_mat = mat.zeros((np.shape(np.mat(tmp_dic.items()[0][1].data).T)[0],len(varord)))
        self.datmat = dat_mat
        for x1 in varord:
            self.datmat[:,x1[1]-1] = np.mat(tmp_dic[x1[0]].data).T
        self.VAR = VAR(varlag,self.datmat,useconst)
        self.VAR.ols()
        self.VAR.ols_comp()
        if spos != 'None':
            self.VAR.do_irf(spos)
            self.VAR.IRF_attr['shock_var'] = varord[[x1[1] for x1 in varord].index(spos-1)][0]
"""***********************************************************"""
##################THE DSGE MODEL CLASS (WORKS)#####################
class DSGEmodel(object):
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
            
        Returns
        -------
        A DSGE model instance
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
        '''
        Model's second intialisation method called by newMOD() function call. The model's
        __init__() method only creates the instance and adds information, but does no
        parsing and computations whatsoever. init2() parses and ready's the model for
        calculating the steady state and the dynamic solution to the model.
        
        init_lev = 0: only parsing, no calculations
        init_lev = 1: parsing and steady state calculations
        init_lev = 2: parsing, steady state calculations and dynamic solution computation
        '''
        initlev = self._initlev

        # Create None tester regular expression
#        _nreg = '^\s*None\s*$'
#        nreg = re.compile(_nreg)

        txtpars = parse_mod(self.modfile)
        self.txtpars = txtpars  # do we need txtpars attached for anything else?
        secs = txtpars.secs # do we need txtpars attached for anything else?

        # Initial population method of model, does NOT need steady states
        self = populate_model_stage_one(self, secs)
        # This is an additional populator which creates subsitution dictionary
        # and uses this to already get rid of @s in the manuall sstate list
        self = populate_model_stage_one_b(self,secs)

        # Attach the data from database
        if self.dbase != None:
            self.getdata(dbase=self.dbase)
        if initlev == 0: 
            return

################## STEADY STATE CALCULATIONS !!! ####################
        self.sssolvers = SSsolvers()
        # Solve for steady-state using fsolve
        # check if the Steady-State Non-Linear system .mod section has an entry
        if any([False if 'None' in x else True for x in secs['manualss'][0]]):
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
        # Check if the Steady States [Closed Form] has an entry
        if any([False if 'None' in x else True for x in secs['closedformss'][0]]):
            if self.switches['ss_suc'] == ['1','1']:
                alldic = {}
                alldic.update(self.sstate)
                alldic.update(self.paramdic)
                intup = (self.manss_sys,alldic)
                self.sssolvers.manss = ManualSteadyState(intup)
                self.sssolvers.manss.solve()
                self.sstate.update(self.sssolvers.manss.sstate)
            else:
                intup = (self.manss_sys,self.paramdic)
                self.sssolvers.manss = ManualSteadyState(intup)
                self.sssolvers.manss.solve()
                self.sstate = self.sssolvers.manss.sstate
        if initlev == 1: return

        # No populate more with stuff that needs steady state
        self = populate_model_stage_two(self, secs)

        # Open the model solution tree branch
        self.modsolvers = MODsolvers()
        ######################## LINEAR METHODS !!! ############################
        # see if there are any log-linearized equations
        if any([False if 'None' in x else True for x in secs['modeq'][0]]):
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
        if any([False if 'None' in x else True for x in secs['focs'][0]]):

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
                     self.mod_name,self.audic,
                     self.numjl,
                     self.nother)
            else:
                intup = (self.numj,
                     self.nendo,self.nexo,
                     self.ncon,self.sigma,
                     self.jAA,self.jBB,
                     self.vardic,self.vdic,
                     self.mod_name,self.audic)
            self.modsolvers.forkleind = ForKleinD(intup)

    ################## 2ND-ORDER NON-LINEAR METHODS !!! ##################
        if any([False if 'None' in x else True for x in secs['vcvm'][0]]) and 'numh' in dir(self):
            # Open the MatKlein2D object
            if 'nlsubsys' in dir(self):
                intup = (self.numj,self.numh,
                     self.nendo,self.nexo,
                     self.ncon,self.sigma,
                     self.jAA,self.jBB,
                     self.vardic,self.vdic,
                     self.mod_name,self.audic,
                     self.numjl,self.numhl,
                     self.nother,sess1)
            else:
                intup = (self.numj,self.numh,
                     self.nendo,self.nexo,
                     self.ncon,self.sigma,
                     self.jAA,self.jBB,
                     self.vardic,self.vdic,
                     self.mod_name,self.audic,
                     sess1)
            self.modsolvers.matklein2d = MatKlein2D(intup)
            # Open the PyKlein2D object
            intup = intup[:-1]
            self.modsolvers.pyklein2d = PyKlein2D(intup)

        if 'jAA' in dir(self):
            self.mkeigv()
        
        # Wrap the paramdic at the end of model initialization
        self.params = dicwrap(self,initlev)

    # html info of model opened with webbrowser
    def info(self):
        '''
        A convience method for collecting and displaying the model's properties in a browser
        using html language.
        '''
        tmplist = glob.glob('tempz*.html')
        for x in tmplist:
            os.remove(x)
        modname = self.mod_name
        secs = self.txtpars.secs
        direc = os.getcwd()
        fd,fpath = mkstemp(prefix='tempz',suffix='.html',dir=direc)
        htmlfile = os.fdopen(fd,'w+b')
        htmlfile.write('<HTML><BODY BGCOLOR="white">\n')
        htmlfile.write('<H2>%s</H2>\n'%modname)
        htmlfile.write('\n\n')
        for x1 in secs['info'][1]:
            htmlfile.write('<P>'+x1+'\n')
        htmlfile.write('<P>'+'<H4>Model Parameters</H4>\n')
        for x1 in secs['params'][1]:
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
        os.system(cmd)
        tmplist = glob.glob('tempz*.html')
        for x in tmplist:
            os.remove(x)
        return 'Model Website opened!'
    # Set the author of the current model
    def setauthor(self,author=None):
        '''
        Convience method for setting the model's author's name.
        '''
        if author == None:
            self.author = 'No author'
        else:
            self.author = author
    # pdflatex the model tex file and open with pdf viewer
    def pdf(self):
        '''
        This will pdflatex the model's tex file and then view it in the system's pdf viewer.
        '''
        modfile = self.modfile
        modfname = modfile.split('.')[0]
        if modfname+'.pdf' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.pdf'))
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.aux' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.aux'))
        if modfname+'.dvi' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.dvi'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))

        # Does the model tex file even exist?
        if modfname+'.tex' not in os.listdir(modfpath):
            print 'Error: The model tex file does not exist!'
            print 'Use model.texed() to create a new one!'
            return

        # First check for Chktex syntax errors
        if 'texer.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),'texer.log'))
        args = '--quiet -v2 -n3 -n25 -n12 -n35'
        cmd = 'chktex '+args+' '+os.path.join(modfpath,modfname+'.tex'+' > texer.log')
        os.system(cmd)
        file = open(os.path.join(os.getcwd(),'texer.log'),'rU')
        errlist = file.readlines()
        if len(errlist)== 0:
            file.close()
            os.remove(os.path.join(os.getcwd(),'texer.log'))
        else:
            print 'There were errors in the model tex file! Abort!\n'

            for x in errlist:
                print x
            file.close()
            os.remove(os.path.join(os.getcwd(),'texer.log'))
            return

        # PDFLatex the tex file and open
        cmd = 'pdflatex '+'-output-directory='+modfpath+' '+os.path.join(modfpath,modfname+'.tex'+' > out.log')
        os.system(cmd)
        cmd2 = pdfpath + ' ' + os.path.join(modfpath,modfname+'.pdf')
        os.system(cmd2)
        if modfname+'.pdf' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.pdf'))
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.aux' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.aux'))
        if modfname+'.dvi' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.dvi'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))
    # Opens the model text file for editing
    def txted(self):
        '''
        A convience method for launching an editor editing the model's associated
        model txt file. If the file gets saved (even if unaltered) the model is
        re-initialized.
        '''
        modfile = self.modfile
        modfname = modfile.split('.')[0]
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))
        cmd = txtedpath+' '+os.path.join(modfpath,self.modfile+' > out.log')
        timestamp0 = os.stat(os.path.join(modfpath,self.modfile))[8]
        os.system(cmd)
        timestamp1 = os.stat(os.path.join(modfpath,self.modfile))[8]
        if timestamp0 != timestamp1:
            self.__init__(self.modfile,self.dbase)
            self.init2()
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))
    # Opens the model latex file for editing
    def texed(self):
        '''
        A convenience method which allows users to launch the model's tex file in an
        editor specified in the configuration settings of pymaclab.
        '''
        modfile = self.modfile
        modfname = modfile.split('.')[0]
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))
        if modfname+'.tex' not in os.listdir(modfpath):
            file = open(os.path.join(modfpath,modfname+'.tex'),'wb')
            file.write('\\documentclass[a4paper,11pt]{article}\n')
            file.write('\\title{%s}\n'% self.mod_name)
            file.write('\\author{%s}\n'% self.author)
            file.write('\\begin{document}\n')
            file.write('\\maketitle\n\n')
            file.write('Start writing your model in Latex here!')
            file.write(8*'\n')
            file.write('\\end{document}\n')
            file.close()
        cmd = texedpath+' '+os.path.join(modfpath,self.modfile[:-3]+'tex'+' > out.log')
        os.system(cmd)
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))
    # Choose to delete the current model's latex file. Your are asked again whether yes/no
    def deltex(self):
        '''
        A simple method which allows you to delete the current tex file associated with this model.
        Will ask for confirmation though, so cannot be used in a fire-and-forget batch file.
        '''
        answ = 0
        while answ not in ['Y','y','N','n']:
            answ = raw_input("REALLY delete this model's tex file ? (Y/N)")
        if answ in ['Y','y']:
            modfile = self.modfile
            modfname = modfile.split('.')[0]
            if modfname+'.tex' in os.listdir(modfpath):
                os.remove(os.path.join(modfpath,modfname+'.tex'))
        elif answ in ['N','n']:
            return
    # Change the 'current view' of model solution algorithms
    def ccv(self,instring=None):
        '''
        This is just a convenience method which allows the current preferred dynamic
        solution method to be stored (attached) here and made accessible without having
        to traverse too deeply into the model instance tree structure.
        '''
        if instring == None:
            return 'Error: Your have not specified a current view name!'
        elif instring not in dir(self.modsolvers):
            return 'Error: Your chosen view does not exist for this model!'
        else:
            self.cv = eval('self.modsolvers.'+instring)
    # Get data method, in case model has not been loaded with database object
    def getdata(self,datafile=None,dbase=None):
        '''
        A separate method useful for importing and attaching to the model data which may
        be used for estimation purposes, which has however not been implemented yet.
        '''
        self.data = {}
        if datafile != None:
            input = open(os.path.join(os.getcwd(),datafile), 'r')
            output = open(os.path.join(os.getcwd(),'tmp_001.csv'), 'w')
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
#            rawdata = P.load('tmp_001.csv',delimiter=',')
            rawdata = np.loadtxt('tmp_001.csv', delimiter=",")
            os.remove(os.getcwd()+'/'+'tmp_001.csv')
            self.rawdata = rawdata
            for x1 in self.dataprop.items():
                if type(x1[1]) != type('x'):
                    exec(x1[0]+'='+str(x1[1]))
                elif type(x1[1]) == type('x'):
                    exec(x1[0]+'='+"'"+x1[1]+"'")

            if np.shape(rawdata)[1] != len(self.varnames):
                print 'DAT IMPORT ERR: Len(varnames) != Len(data)'
                data = TimeSeries(rawdata,start_date=ts.Date(freq=freq,year=year,month=month))
            else:
                data = TimeSeries(rawdata,start_date=ts.Date(freq=freq,year=year,month=month))
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
        instring = deepcopy(cinstring)

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

        reva = re.compile('(?P<pre>^[^a-zA-Z^_]{0,1}|[^a-zA-Z^_]{1,1})('+iti+varp+ti+')')

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
        An unparallelized method using native Python and Sympycore in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        
        Parameters
        ----------
        self: object instance
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.jdic: attaches the analytical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.hdic: attaches the analytical 3D Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
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
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)

        # Create substitution var list and dictionary
        tmpli = []
        for i,x in enumerate(intup):
            nolen = len(str(i))
            inds = (5-nolen)*'0'+str(i)
            tmpli.append([x,'SUB'+inds])
            dicli = dict(tmpli)
            dicli2 = dict([[x[1],x[0]] for x in tmpli])
        self.subs_li = deepcopy(tmpli)
        self.var_li = [x[0] for x in tmpli]
        self.subs_dic = deepcopy(dicli)
        self.subs_dic2 = deepcopy(dicli2)

        func = []
        subli = []
        symdic = {}
        for x in dicli.values():
            symdic[x] = sympy.Symbol(x)
        locals().update(symdic)
#NOTE: don't do this?
# or don't modify _existing_ contents?
# see: http://docs.python.org/library/functions.html#locals
#also locals() is no longer a dictionary in Python 3
#NOTE: just use the dictionary itself
        for x in self.paramdic.keys():
            locals()[x] = sympy.Symbol(x)
        for x in self.sstate.keys():
            locals()[x] = sympy.Symbol(x)
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
            elog = re.compile('LOG\(')
            while elog.search(str_tmp):
                ma = elog.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympy.log('+str_tmp[poe:]
            eexp = re.compile('EXP\(')
            while eexp.search(str_tmp):
                ma = eexp.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympy.exp('+str_tmp[poe:]
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
        mreg = re.compile(_mreg)
        func2 = []
        for i,x in enumerate(func):
            func2.append(func[i])
            liobj = list(mreg.finditer(str(func2[i])))
            doneli = []
            for ma in reversed(liobj):
                suba = ma.group()
                if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
                   and suba not in doneli:
                    func2[i] = func2[i].subs(sympy.Symbol(suba),sympy.exp(sympy.Symbol(suba)))
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
                    self.func2[i1] = sympy.Symbol(vmavar)*(self.func2[i1])

        # Create list with log of steady states
        evalli = []
        alldic = deepcopy(self.paramdic)
        alldic.update(self.sstate)
        for x in intup:
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            if 'log' in moddic[lookdic[vari]]:
                evalli.append(np.log(alldic[vari+'_bar']))
            else:
                evalli.append(alldic[vari+'_bar'])
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Make 2D symbolic and numeric Jacobian
        def mkjac(jrows=jrows,jcols=jcols):
            rdic = dict([[x,'0'] for x in range(jcols)])
            jdic = dict([[x,deepcopy(rdic)] for x in range(jrows)])
            jdicc = deepcopy(jdic)
            jcols = len(jdic[0])
            jrows = len(jdic)
            carry_over_dic = {}
            numj = mat.zeros((jrows,jcols))
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            for x in range(jrows):
                jdicc[x] = {}
                for y in range(jcols):
                    jdic[x][y] = func2[x].diff(symdic[tmpli[y][1]])
                    suba_dic = self.subs_dic2
                    differo_var = str(symdic[tmpli[y][1]])
                    differo_var = suba_dic[differo_var]
                    carry_over_dic[y] = differo_var
                    jdicc[x][differo_var] = str(jdic[x][y])
                    if mreg.search(jdicc[x][differo_var]):
                        for keyo in suba_dic.keys():
                            jdicc[x][differo_var] = jdicc[x][differo_var].replace(keyo,suba_dic[keyo])
                    else:
                        jdicc[x][differo_var] = jdicc[x][differo_var]
                    evalfo = jdic[x][y].evalf()
                    if 'log(' not in str(evalfo) and 'exp(' not in str(evalfo):
                        numj[x,y] = eval(str(evalfo))
                    elif 'exp(' in str(evalfo):
                        numj[x,y] = eval(str(evalfo).replace('exp(','np.exp('))
                    elif 'log(' in str(evalfo):
                        numj[x,y] = eval(str(evalfo).replace('log(','np.log('))
            # Take out the elements from the variable substitution equations
            lengor = len(self.nlsys_list)
            for keyo in jdicc.keys():
                if keyo > (lengor-1):
                    jdicc.pop(keyo)
            return numj,jdic,jdicc,carry_over_dic

        # Now make 3D symbolic and numeric Hessian
        def mkhes(jrows=jrows,jcols=jcols,trans_dic=None):
            rdic = dict([[x,'0'] for x in range(jcols)])
            rdic1 = dict([[x,deepcopy(rdic)] for x in range(jcols)])
            hdic = dict([[x,deepcopy(rdic1)] for x in range(jrows)])
            hdicc = deepcopy(hdic)
            hcols = len(hdic[0])
            hrows = len(hdic[0])*len(hdic)
            numh = mat.zeros((hrows,hcols))
            jdic = self.jdic
            jdicc = self.jdicc
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            count = 0
            for x in range(jrows):
                hdicc[x] = {}
                for y in range(jcols):
                    hdicc[x][trans_dic[y]] = {}
                    for z in range(jcols):
                        hdic[x][y][z] = jdic[x][y].diff(symdic[tmpli[z][1]])
                        suba_dic = self.subs_dic2
                        differo_var = str(symdic[tmpli[z][1]])
                        differo_var = suba_dic[differo_var]
                        hdicc[x][trans_dic[y]][differo_var] = str(hdic[x][y][z])
                        if mreg.search(hdicc[x][trans_dic[y]][differo_var]):
                            for keyo in suba_dic.keys():
                                hdicc[x][trans_dic[y]][differo_var] = hdicc[x][trans_dic[y]][differo_var].replace(keyo,suba_dic[keyo])
                        else:
                            hdicc[x][trans_dic[y]][differo_var] = hdicc[x][trans_dic[y]][differo_var]
                        evalfo = hdic[x][y][z].evalf()
                        if 'log(' not in str(evalfo) and 'exp(' not in str(evalfo):
                            numh[count,z] = eval(str(evalfo))
                        elif 'exp(' in str(evalfo):
                            numh[count,z] = eval(str(evalfo).replace('exp(','np.exp('))
                        elif 'log(' in str(evalfo):
                            numh[count,z] = eval(str(evalfo).replace('log(','np.log('))
                    count = count + 1
            # Take out the elements from the variable substitution equations
            lengor = len(self.nlsys_list)
            for keyo in hdicc.keys():
                if keyo > (lengor-1):
                    hdicc.pop(keyo)
            return numh,hdic,hdicc

        self.numj,self.jdic,self.jdicc,carry_over_dic = mkjac()
        numj = self.numj
        if mk_hessian:
            self.numh,self.hdic,self.hdicc = mkhes(trans_dic=carry_over_dic)
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
        An parallelized method using native Python and Sympycore in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        This is the parallelized version of method self.mkjahe using the Python pp library.
        
        Parameters
        ----------
        self: object instance
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.jdic: attaches the analytical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.hdic: attaches the analytical 3D Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
        '''
        # import local sympy
        import sympy

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
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)

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
            symdic[x] = sympy.Symbol(x)
        locals().update(symdic)
        for x in self.paramdic.keys():
            locals()[x] = sympy.Symbol(x)
        for x in self.sstate.keys():
            locals()[x] = sympy.Symbol(x)
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
            elog = re.compile('LOG\(')
            while elog.search(str_tmp):
                ma = elog.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympy.log('+str_tmp[poe:]
            eexp = re.compile('EXP\(')
            while eexp.search(str_tmp):
                ma = eexp.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympy.exp('+str_tmp[poe:]
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
        mreg = re.compile(_mreg)
        func2 = []
        for i,x in enumerate(func):
            func2.append(func[i])
            liobj = list(mreg.finditer(str(func2[i])))
            doneli = []
            for ma in reversed(liobj):
                suba = ma.group()
                if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
                   and suba not in doneli:
                    func2[i] = func2[i].subs(sympy.Symbol(suba),sympy.exp(sympy.Symbol(suba)))
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
                    self.func2[i1] = sympy.Symbol(vmavar)*(self.func2[i1])

        # Create list with log of steady states
        evalli = []
        alldic = {}
        alldic.update(self.paramdic)
        alldic.update(self.sstate)
        for x in intup:
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            if 'log' in moddic[lookdic[vari]]:
                evalli.append(np.log(alldic[vari+'_bar']))
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
                hdic = dict([[x,deepcopy(rdic)] for x in range(jcols)])
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
                evalfo = jdic[y].evalf()
                if 'log(' not in str(evalfo) and 'exp(' not in str(evalfo):
                    numj[0,y] = eval(str(evalfo))
                elif 'exp(' in str(evalfo):
                    numj[0,y] = eval(str(evalfo).replace('exp(','np.exp('))
                elif 'log(' in str(evalfo):
                    numj[0,y] = eval(str(evalfo).replace('log(','np.log('))
                if mk_hessian:
                    for z in range(jcols):
                        hdic[y][z] = jdic[y].diff(symdic[tmpli[z][1]])
                        evalfo2 = hdic[y][z].evalf()
                        if 'log(' not in str(evalfo2) and 'exp(' not in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2))
                        elif 'exp(' in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2).replace('exp(','np.exp('))
                        elif 'log(' in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2).replace('log(','np.log('))
                    count = count + 1
            if mk_hessian:
                return (numj,jdic,numh,hdic)
            else:
                return (numj,jdic)

        # Start parallel Python job server
        ppservers = ()
        inputs = [x for x in xrange(len(self.func2))]
        job_server = pp.Server(ncpus,ppservers=ppservers)
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
                numj = mat.vstack((numj,job()[0]))
                jdic[i1+1] = job()[1]
                numh = mat.vstack((numh,job()[2]))
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
                numj = mat.vstack((numj,job()[0]))
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
        '''
        A method using mlabwrap to call an external Matlab code in oder to calculate
        the numerical Jacobian and Hessian of the model.
        
        Parameters
        ----------
        self: object instance
        msess: active mlabwrap session to be used
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
        '''

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
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)

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
            elog = re.compile('LOG\(')
            while elog.search(str_tmp):
                ma = elog.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'log('+str_tmp[poe:]
            eexp = re.compile('EXP\(')
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
        mreg = re.compile(_mreg)
        func2 = []
        for i,x in enumerate(func):
            func2.append(func[i])
            liobj = list(mreg.finditer(str(func2[i])))
            doneli = []
            for ma in reversed(liobj):
                suba = ma.group()
                if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
                   and suba not in doneli:
                    func2[i] = func2[i].subs(sympy.Symbol(suba),sympy.exp(sympy.Symbol(suba)))
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
                    self.func2[i1] = sympy.Symbol(vmavar)*(self.func2[i1])

                    vma = self.vreg(patup,vari0,False,'min')
                    vmavar = vma[2][1]+'_bar'
                    self.func2[i1] = sympy.Symbol(vmavar)*(self.func2[i1])


        # Transform python exponentiation into matlab exponentiation
        _mreg='\*{2,2}'
        mreg=re.compile(_mreg)
        for i,x in enumerate(self.func2):
            self.func2[i] = mreg.sub('^',self.func2[i])

        # Create matlab function
        if 'mfunc.m' in os.listdir(os.path.join(mlabpath,'Klein')):
            os.remove(os.path.join(mlabpath,'Klein/mfunc.m'))

        mfunc = open(os.path.join(mlabpath,'Klein/mfunc.m'),'w')
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
        inmat = mat.zeros((len(tmpli),1))
        alldic={}
        alldic.update(self.paramdic)
        alldic.update(self.sstate)
        for i,x in enumerate(tmpli):
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            inmat[i,0] = alldic[vari+'_bar']
            if 'log' in moddic[lookdic[vari]]:
                inmat[i,0] = np.log(inmat[i,0])


        #Make Jacobian and Hessian
        sess1 = msess
        directory = os.path.join(mlabpath,'Klein')
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

        self.jAA = np.matrix(self.numj[:,:int(len(intup)/2)])
        self.jBB = np.matrix(-self.numj[:,int(len(intup)/2):])

        del self.func1
        del self.func2

        os.remove(os.path.join(mlabpath,'Klein/mfunc.m'))

###NUMERICAL ONLY JACOBIAN AND HESSIAN METHODS - FOR UPDATING###############
    # The unparallelized Jacobian, Hessian computation method
    def mkjahen(self):
        '''
        An unparallelized method using native Python and Sympycore in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        This is the corresponding method to "self.mkjahe()" but is only called when
        some parameters were changed and the model needs dynamic updating. In particular
        in this case the expensive computation of the analytical expressions can be
        avoided, as only parameters affecting the steady state may have been altered.
        Uses existing self.jdic for Jacobian and self.hdic for Hessian of model.
        
        Parameters
        ----------
        self: object instance
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
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
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)

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
                evalli.append(np.log(alldic[vari+'_bar']))
            else:
                evalli.append(alldic[vari+'_bar'])
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Make 2D symbolic and numeric Jacobian
        def mkjac(jrows=jrows,jcols=jcols):
            jdic = self.jdic
            numj = mat.zeros((jrows,jcols))
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            for x in range(jrows):
                for y in range(jcols):
                    evalfo = jdic[x][y].evalf()
                    if 'exp(' not in str(evalfo) and 'log(' not in str(evalfo):
                        numj[x,y] = eval(str(evalfo))
                    elif 'exp(' in str(evalfo):
                        numj[x,y] = eval(str(evalfo).replace('exp(','np.exp('))
                    elif 'log(' in str(evalfo):
                        numj[x,y] = eval(str(evalfo).replace('log(','np.log('))
            return numj

        # Now make 3D symbolic and numeric Hessian
        def mkhes(jrows=jrows,jcols=jcols):
            hdic = self.hdic
            hrows = jrows*jcols
            numh = mat.zeros((hrows,jcols))
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            count = 0
            for x in range(jrows):
                for y in range(jcols):
                    for z in range(jcols):
                        evalfo = hdic[x][y][z].evalf()
                        if 'exp(' not in str(evalfo) and 'log(' not in str(evalfo):
                            numh[count,z] = eval(str(evalfo))
                        elif 'exp(' in str(evalfo):
                            numh[count,z] = eval(str(evalfo).replace('exp(','np.exp('))
                        elif 'log(' in str(evalfo):
                            numh[count,z] = eval(str(evalfo).replace('log(','np.log('))
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
        A parallelized method using native Python and Sympycore in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        This is the corresponding method to "self.mkjahe()" but is only called when
        some parameters were changed and the model needs dynamic updating. In particular
        in this case the expensive computation of the analytical expressions can be
        avoided, as only parameters affecting the steady state may have been altered.
        Uses existing self.jdic for Jacobian and self.hdic for Hessian of model.
        
        Parameters
        ----------
        self: object instance
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
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
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)

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
                evalli.append(np.log(alldic[vari+'_bar']))
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
                evalfo = jdic[lcount][y].evalf()
                if 'exp(' not in str(evalfo) and 'log(' not in str(evalfo):
                    numj[0,y] = eval(str(evalfo))
                elif 'exp(' in str(evalfo):
                    numj[0,y] = eval(str(evalfo).replace('exp(','np.exp('))
                elif 'log(' in str(evalfo):
                    numj[0,y] = eval(str(evalfo).replace('log(','np.log('))
                if mk_hessian:
                    for z in range(jcols):
                        evalfo2 = hdic[lcount][y][z].evalf()
                        if 'exp(' not in str(evalfo2) and 'log(' not in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2))
                        elif 'exp(' in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2).replace('exp(','np.exp('))
                        elif 'log(' in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2).replace('log(','np.log('))
                    count = count + 1
            if mk_hessian:
                return (numj,numh)
            else:
                return numj

        # Start parallel Python job server
        ppservers = ()
        inputs = [x for x in xrange(len(self.func2))]
        job_server = pp.Server(ncpus,ppservers=ppservers)
        from sympycore import exp
        from sympycore import log
        imports = ('numpy','copy','numpy.matlib',)
        jobs = [job_server.submit(mkjaheseq,(input,self.func2,jcols,tmpli,self.paramdic,self.sstate,evaldic,mk_hessian,exp,log,self.jdic,self.hdic),(),imports) for input in inputs]
        if mk_hessian:
            job_0 = jobs[0]
            numj = job_0()[0]
            numh = job_0()[1]
            for i1,job in enumerate(jobs[1:len(jobs)]):
                numj = mat.vstack((numj,job()[0]))
                numh = mat.vstack((numh,job()[1]))
            self.numj = numj
            self.numh = numh
        else:
            job_0 = jobs[0]
            numj = job_0()
            for i1,job in enumerate(jobs[1:len(jobs)]):
                numj = mat.vstack((numj,job()))
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
        '''
        Method useful for completely re-initializing model based on a different or else
        altered external modfile. It is also perceivable to load the existing modfile at
        runtime, change it in a program, save it as a temporary file and then re-load it
        using this method.
        '''
        self.txtpars.__init__(self.modfile)

    # Method updating model after anything has been changed manually, like paramvalue, etc !
    def updm(self):
        '''
        Method useful for upating model manually if any parameter value has been changed
        manually as well.
        '''
        # Create None tester regular expression
        _nreg = '^\s*None\s*$'
        nreg = re.compile(_nreg)
################## STEADY STATE CALCULATIONS !!! ####################
        self.sssolvers = SSsolvers()
        # Solve for steady-state using fsolve
        if sum([nreg.search(x)!=None for x in self.txtpars.secs['manualss'][0]]) == 0:
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
        if sum([nreg.search(x)!=None for x in self.txtpars.secs['closedformss'][0]]) == 0:
            intup = (self.manss_sys,self.paramdic)
            self.sssolvers.manss = ManualSteadyState(intup)
            self.sssolvers.manss.solve()
            self.sstate = self.sssolvers.manss.sstate
        if self._initlev == 1: return

        # No populate more with stuff that needs steady state
        secs = self.txtpars.secs
        self = populate_model_stage_two(self, secs)

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
                     self.mod_name,
                     self.numjl,
                     self.nother)
            else:
                intup = (self.numj,
                     self.nendo,self.nexo,
                     self.ncon,self.sigma,
                     self.jAA,self.jBB,
                     self.vardic,self.vdic,
                     self.mod_name)
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
                         self.mod_name,
                         self.numjl,self.numhl,
                         self.nother,sess1)
                else:
                    intup = (self.numj,self.numh,
                         self.nendo,self.nexo,
                         self.ncon,self.sigma,
                         self.jAA,self.jBB,
                         self.vardic,self.vdic,
                         self.mod_name,
                         sess1)
                self.modsolvers.matklein2d = MatKlein2D(intup)
                # Open the PyKlein2D object
                intup = intup[:-1]
                self.modsolvers.pyklein2d = PyKlein2D(intup)

        if 'jAA' in dir(self):
            self.mkeigv()

    def mkeigv(self):
        '''
        A method useful for calculating the Eigenvalues as in Blanchard and Kahn
        for inspection. Recall the Eigenvalue Rule of this paper which needs to be
        meet for the model to be solvable.
        '''
        AA = self.jAA
        BB = self.jBB
        try:
            CC = BB.I*AA
            self.eigvals = scipyeig(CC)[0]
        except:
            pass

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
