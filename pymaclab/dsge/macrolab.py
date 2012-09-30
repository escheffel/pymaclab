'''
.. module:: macrolab
   :platform: Linux
   :synopsis: The macrolab module contains the the most important classes for doing work with DSGE models. In particular it supplies
              the DSGEmodel class containing most of the functionality of DSGE model instances. The module also contains the TSDataBase
              class which was supposed to be an advanced data carrier class to be passed to the DSGE model instance for estimation
              purposes, but this is deprecated and will probably be replaced with a pandas data frame in the near future.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''

#from __future__ import division
# Have taken out the dependency from this module by commenting out
# The TimeSeries Database as well as the get_data method on the DSGE class
# In future versions this has to be replaced with Pandas instead
#from scikits import timeseries as ts
import copy
from copy import deepcopy
import os
import sys
import re
from helpers import now_is
import numpy as np
from numpy import matlib as mat
from scipy.linalg import eig as scipyeig
from tempfile import mkstemp
import time
import datetime
import glob
# Import Sympycore, now comes supplied and installed with pymaclab
import sympycore
# Import PP, now comes supplied and installed with pymaclab
import pp

#NOTE: Imports from the refactor
from pymaclab.filters._hpfilter import hpfilter
from pymaclab.filters._bkfilter import bkfilter
from ..stats.var import VAR #TODO: remove for statsmodels version
from solvers.steadystate import SSsolvers, ManualSteadyState, Fsolve
from parsers._modparser import parse_mod
from parsers._dsgeparser import populate_model_stage_one,populate_model_stage_one_a,\
     populate_model_stage_one_b,populate_model_stage_one_bb,populate_model_stage_two
from updaters.one_off import Updaters, dicwrap, dicwrap_deep, listwrap, matwrap
from updaters.queued import Updaters_Queued, dicwrap_queued, dicwrap_deep_queued, listwrap_queued,\
     matwrap_queued, Process_Queue, queue

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
# Open mlab session
if use_matlab:
    sess1 = mlabraw.open('matlab - nojvm -nosplash')

# Empty locdic for the locate helper function
locdic = {}

'''
################THE TIMESERIES DATABASE CLASS (WORKS)#################
class TSDataBase:
    """
    This is the Time Series Database class.
    """
    def __init__(self):
        self.datdic = {}
        self.modif = {}

    def nameimDatStr(self,filename='none'):
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
        tsoutf = hpfilter(tsinf,np.zeros((np.shape(tsinf)[0],3)),
            np.shape(tsinf)[0],lam,0)
        tsoutf = ts.time_series\
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
'''
##################THE DSGE MODEL CLASS (WORKS)#####################
class DSGEmodel(object):
    '''
    This is the macrolab DSGEmodel class. It is the main class of the packages and instantiates
    DSGE model instances which possess many features such as model file parsers, solvers, etc. The __init__
    function first called mostly attaches the passed arguments to the DSGE model instance in form of private _X
    data fields.
    
    .. note::
    
       Notice that the various init method, i.e. init1, init1a, etc. are not called here, but they are called externally in the
       pymaclab package when instantiating a new DSGE model using pymaclab.newMOD().
    
    :param ffile:             The absolute path to the PyMacLab DSGE model file to be parsed on instantiation
    :type ffile:              str
    :param dbase:             A database with time series data for estimation purposes. Not implemented at the moment
    :type standard:           unknown
    :param initlev:           Takes values 0,1,2. Determines how deep the instantiation cascades through all methods.
                              If 0 then the file is parsed and the instance is only *prepared* for steady state solving.
                              If 1 then the file is parsed, the SS is automatically computed and instance is *prepared*
                              for dynamic solving. If 2 then the instance is solved all the way.
    :type initlev:            int
    :param mesg:              If True then lots of diagnostics are printed to the screen during instantiation.
    :type mesg:               bool
    :param ncpus:             The number of CPU core to be employed, defaults to 1. But 'auto' can also be used for detection
    :type ncpus:              int|str
    :param mk_hessian:        Whether the Hessian should be computed as this is expensive.
    :type mk_hessian:         bool
    :param use_focs:          Should the FOCs be used directly to look for the steady state? Must use list|tuple to pick equations.
    :type use_focs:           tuple
    :param ssidic:            A Python ssidic with the initial starting values for solving for the SS numerically
    :type ssidic:             dic
        
    :return self:             *(dsge_inst)* - A populated DSGE model instance with fields and methods

    '''    
    # Initializes the absolute basics, errors difficult to occur
    def __init__(self,ffile=None,dbase=None,initlev=2,mesg=False,ncpus=1,mk_hessian=True,use_focs=False,ssidic=None):
        self._initlev = initlev #TODO: remove this as an option
        self._mesg = mesg
        self._ncpus = ncpus
        self._mk_hessian = mk_hessian
        self._use_focs = use_focs
        self._ssidic = ssidic
        if self._use_focs and self._ssidic != None:
            self.ssidic = deepcopy(self._ssidic)
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
    def init1(self):
        '''
        The init1 method. Model population proceeds from the __init__function here. In particular the data gets read in
        (not implemented at the moment) and the model parsing begins.
        
        .. note::
           The self.vardic gets created and the manual as well as
           the numerical steady state sections gets parsed and attached to the DSGE model instance. So the most import fields created
           here using function populate_model_stage_one() are:
           
             * self.vardic - variable names dictionary with transform and filtering info
             * self.mod_name - short name of the DSGE model from mod file
             * self.mod_desc - longer model description from mod file
             * self.paramdic - dic of defined parameters with their numerical values
             * self.manss_sys - list of equations from the closed form steady state section
             * self.ssys_list - list of equations from the numerical steady state section
           
           Also the updaters and updaters_queued branches are opened here and the self.vardic gets wrapped for dynamic updating
           behaviour.
        
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst
        
        '''
        initlev = self._initlev
        ncpus = self._ncpus
        mk_hessian = self._mk_hessian
        mesg = self._mesg
        # Attach the data from database
        if self.dbase != None:
            self.getdata(dbase=self.dbase)        
        if mesg: print "INIT: Instantiating DSGE model with INITLEV="+str(initlev)+" and NCPUS="+str(ncpus)+"..."
        if mesg: print "INIT: Attaching model properties to DSGE model instance..."

        # Create None tester regular expression
        #        _nreg = '^\s*None\s*$'
        #        nreg = re.compile(_nreg)
        
        if mesg: print "INIT: Parsing model file into DSGE model instance..."
        txtpars = parse_mod(self.modfile)
        self.txtpars = txtpars  # do we need txtpars attached for anything else?
        secs = txtpars.secs # do we need txtpars attached for anything else?
        if mesg: print "INIT: Extraction of info into DSGE model instance Stage [1]..."
        # Initial population method of model, does NOT need steady states
        self = populate_model_stage_one(self, secs)

        # Open updaters path
        self.updaters = Updaters()        
        # Open the updaters_queued path
        self.updaters_queued = Updaters_Queued()
        # Add the empty queue
        self.updaters_queued.queue = queue

        # Wrap the vardic
        self.updaters.vardic = dicwrap_deep(self,'self.vardic',initlev)
        self.updaters_queued.vardic = dicwrap_deep_queued(self,'self.vardic',initlev)
        
    def init1a(self):
        '''
        The init1a method. Model population proceeds from the init1 method here. The only field which get created here is the raw
        (i.e. unsubstituted) substitution dictionary.
        
        .. note::
            Field which are created here using the function populate_model_stage_one_a() which in turn calls mk_subs_dic():
           
             * self.nlsubs_raw1 - a list of the @items and their replacements
             * self.nlsubsdic - the above just expressed as a keyed list
        
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''
        mesg = self._mesg
        secs = self.txtpars.secs
        # This is an additional populator which creates subsitution dictionary
        # and uses this to already get rid of @s in the manuall sstate list
        if mesg: print "INIT: Extraction of info into DSGE model instance Stage [2]..."
        # This function only creates the raw substitution dictionary and list from modfile
        self = populate_model_stage_one_a(self,secs)
        
    def init1b(self):
        '''
        The init1b method. Model population proceeds from the init1a method here.

        .. note::
        
           The only thing which gets done here purposefully *after* calling init1a() is to wrap these fields:
           
             * self.nlsubsdic - the above just expressed as a keyed list
             * self.paramdic - dic of defined parameters with their numerical values
             
           in order to give them dynamic updating behaviour. No more is done in this init method call.
        
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''        
        initlev = self._initlev
        secs = self.txtpars.secs
        self = populate_model_stage_one_b(self,secs)
        
        # Wrap the nlsubsdic
        self.updaters_queued.nlsubsdic = dicwrap_queued(self,'self.nlsubsdic',initlev)
        # Wrap the paramdic
        self.updaters_queued.paramdic = dicwrap_queued(self,'self.paramdic',initlev)
        

        # Wrap the nlsubsdic
        self.updaters.nlsubsdic = dicwrap(self,'self.nlsubsdic',initlev)
        # Wrap the paramdic
        self.updaters.paramdic = dicwrap(self,'self.paramdic',initlev)        
        
    def init1c(self):
        '''
        The init1c method. Model population proceeds from the init1b method here. We call populate_model_stage_one_bb() which does quite
        a bit of substitution/replacement of the @-prefixed variables. It does this in the numerical and closed form steady state
        calculation sections, as well as in the actual FOCs themselves.

        .. note::
        
           *After* that we can wrap various fields for updating which are:
           
             * self.foceqs - list of firs-order conditions with @ replacements done
             * self.manss_sys - the list of equations from closed form SS section
             * self.syss_list - the list of equations from numerical SS section
        
           Notice also that here we replace or generate fields in case the FOCs are supposed to be used
           directly in SS calculation because the "use_focs" parameter was not passed empty.
        
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''
        initlev = self._initlev
        secs = self.txtpars.secs
        mesg = self._mesg
        if mesg: print "INIT: Substituting out @ variables in steady state sections..."
        self = populate_model_stage_one_bb(self,secs)
        
        # Wrap foceqs
        self.updaters.foceqs = listwrap(self,'self.foceqs',initlev)
        # Wrap foceqs
        self.updaters_queued.foceqs = listwrap(self,'self.foceqs',initlev)
        
        # Allow for numerical SS to be calculated using only FOCs
        if self._use_focs and self._ssidic != None:
            list_tmp = []
            for elem in self._use_focs:
                list_tmp.append(self.foceqss[elem])
            self.ssys_list = deepcopy(list_tmp)
            self.ssidic = copy.deepcopy(self._ssidic)
            
        # Wrap manss_sys
        if 'manss_sys' in dir(self):
            self.updaters.manss_sys = listwrap(self,'self.manss_sys',initlev)
            self.updaters_queued.manss_sys = listwrap_queued(self,'self.manss_sys',initlev)
        # Wrap ssys_list
        if 'ssys_list' in dir(self):
            self.updaters.ssys_list = listwrap(self,'self.ssys_list',initlev)
            self.updaters_queued.ssys_list = listwrap_queued(self,'self.ssys_list',initlev)


    def init2(self):
        '''
        The init2 method. Model population proceeds from the init1c method here. In this initialisation method
        the only thing which is being done is to open up the sssolvers branch and pass down required objects to
        the manuals closed from solver and the numerical root-finding solver depending on whether information for
        this has been included in the DSGE model file. *No* attempt is made at solving for the steady state, the
        respective solvers are only being *prepared*.
        
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst
        
        '''
        mesg = self._mesg
        secs = self.txtpars.secs
        initlev = self._initlev
        if mesg: print "INIT: Preparing DSGE model instance for steady state solution..."

        # Attach the steady state class branch, and add Fsolve if required but do no more !
        self.sssolvers = SSsolvers()

        # check if the Steady-State Non-Linear system .mod section has an entry
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) or self._use_focs:
            intup = (self.ssys_list,self.ssidic,self.paramdic)
            self.sssolvers.fsolve = Fsolve(intup)
        # check if the Steady-State Non-Linear system closed form has an entry
        #  we do this here with ELIF because we only want to set this up for solving if manualss does not exist
        elif all([False if 'None' in x else True for x in secs['closedformss'][0]]):
            alldic = {}
            alldic.update(self.paramdic)
            intup = (self.manss_sys,alldic)
            self.sssolvers.manss = ManualSteadyState(intup)

    def init3(self):
        '''
        Initialisation method init3 is quite complex and goes through a number of logical tests in order to determine
        how to solve for the steady state of the model.
        
        .. note::
        
           There are 5 different ways a DSGE model can obtain its steady state solution depending on what information has
           been provided:
           
             1) Information has been provided using the "use_focs" parameter to use FOCs directly
             2) Information has only been provided in the numerical SS section
             3) Information has only been provided in the closed form SS section
             4) Both CF-SS and NUM-SS info are present and NUM-SS is subset if CF-SS
             5) Both CF-SS and NUM-SS info are present and CF is residual
             
           These options are better explained in the documentation to PyMacLab in the steady state solver section.
           
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst
        
        '''
        txtpars = self.txtpars
        secs = txtpars.secs
        initlev = self._initlev
################################## STEADY STATE CALCULATIONS !!! #######################################
        if self._mesg: print "INIT: Attempting to find DSGE model's steady state automatically..."
        # ONLY NOW try to solve !
        ##### OPTION 1: There is only information externally provided and we are using FOCs
        if not all([False if 'None' in x else True for x in secs['closedformss'][0]]) and\
           not all([False if 'None' in x else True for x in secs['manualss'][0]]) and self._use_focs and self._ssidic != None:
            self.sssolvers.fsolve.solve()
            
        ##### OPTION 2: There is only information provided in the numerical section NOT in closed-form
        if not all([False if 'None' in x else True for x in secs['closedformss'][0]]) and\
           all([False if 'None' in x else True for x in secs['manualss'][0]]) and not self._use_focs and self._ssidic == None:
            if self._mesg: print "SS: ONLY numerical steady state information supplied...attempting to solve SS..."
            self.sssolvers.fsolve.solve()

        ##### OPTION 3: There is only information on closed-form steady state, BUT NO info on numerical steady state
        # Solve using purely closed form solution if no other info on model is available
        if all([False if 'None' in x else True for x in secs['closedformss'][0]]) and\
           not all([False if 'None' in x else True for x in secs['manualss'][0]]):
            if self._mesg: print "SS: ONLY CF-SS information supplied...attempting to solve SS..."
            alldic = {}
            alldic.update(self.paramdic)
            intup = (self.manss_sys,alldic)
            self.sssolvers.manss = ManualSteadyState(intup)
            self.sssolvers.manss.solve()
            self.sstate = {}
            self.sstate.update(self.sssolvers.manss.sstate)
        ##### OPTION 4: There is information on closed-form AND on numerical steady state
        # Check if the numerical and closed form sections have entries
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           all([False if 'None' in x else True for x in secs['closedformss'][0]]):
            # Create unordered Set of closed from solution variables
            manss_set = set()
            for elem in self.manss_sys:
                manss_set.add(elem.split('=')[0].lstrip().rstrip())
            # If ssidic is not empty we need to make sure it perfectly overlaps with manss_set in order to replace ssidic
            if self.ssidic != {}:
                numss_set = set()
                for keyo in self.ssidic.keys():
                    numss_set.add(keyo)
                ##### OPTION 4a: If there is an ssidic and its keys are subset of manss_set, the use as suggestion for new ssi_dic
                if numss_set.issubset(manss_set):
                    if self._mesg: print "SS: CF-SS and NUM-SS (overlapping) information information supplied...attempting to solve SS..."
                    alldic = {}
                    alldic.update(self.paramdic)
                    intup = (self.manss_sys,alldic)
                    self.sssolvers.manss = ManualSteadyState(intup)
                    self.sssolvers.manss.solve()
                    for keyo in self.ssidic.keys():
                        self.ssidic[keyo] = self.sssolvers.manss.sstate[keyo]
            ##### OPTION 4b: ssidic is empty, so we have to assumed that the variables in closed form are suggestions for ssidic
            # If it is empty, then just compute the closed form SS and pass to ssidic as starting value
            elif self.ssidic == {}:
                if self._mesg: print "SS: CF-SS and NUM-SS (empty ssdic) information information supplied...attempting to solve SS..."
                alldic = {}
                alldic.update(self.paramdic)
                intup = (self.manss_sys,alldic)
                self.sssolvers.manss = ManualSteadyState(intup)
                self.sssolvers.manss.solve()
                self.ssidic.update(self.sssolvers.manss.sstate)
                # Test at least if the number of ssidic vars equals number of equations
                if len(self.ssidic.keys()) != len(self.ssys_list):
                    print "Error: Number of variables in initial starting values dictionary != number of equations"
                    sys.exit()
        ######## Finally start the numerical root finder with old or new ssidic from above
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           not all([False if 'None' in x else True for x in secs['closedformss'][0]]) and not self._use_focs:            
            self.sssolvers.fsolve.solve()
        elif all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           all([False if 'None' in x else True for x in secs['closedformss'][0]]) and not self._use_focs:           
            self.sssolvers.fsolve.solve()

        if all([False if 'None' in x else True for x in secs['manualss'][0]]):
            if self.sssolvers.fsolve.ier == 1:
                self.sstate = self.sssolvers.fsolve.fsout
                self.numssdic = self.sssolvers.fsolve.fsout
                # Attach solutions to intial variable dictionaries, for further analysis
                self.ssidic_modfile = deepcopy(self.ssidic)
                # Update old ssidic with found solutions
                self.ssidic = deepcopy(self.sssolvers.fsolve.fsout)
                self.sssolvers.fsolve.ssi = self.ssidic
                self.switches['ss_suc'] = ['1','1']
                if self._mesg: print "INIT: Steady State of DSGE model found (SUCCESS)..."
            else:
                self.switches['ss_suc'] = ['1','0']
                if self._mesg: print "INIT: Steady State of DSGE model not found (FAILURE)..."

        ########## Here we are trying to merge numerical SS solver's results with result closed-form calculations, if required
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           all([False if 'None' in x else True for x in secs['closedformss'][0]]):
            if self._mesg: print "INIT: Merging numerical with closed form steady state if needed..."
        # Check for existence of closedform AND numerical steady state
        # We need to stop the model instantiation IFF numerical solver was attempted but failed AND closed form solver depends on it.
        if all([False if 'None' in x else True for x in secs['closedformss'][0]]) and\
           all([False if 'None' in x else True for x in secs['manualss'][0]]) and self.ssidic != {}:
            if self.switches['ss_suc'] == ['1','0']:
                print "ERROR: You probably want to use numerical steady state solution to solve for RESIDUAL closed form steady states."
                print "However, the numerical steady state solver FAILED to find a root, so I am stopping model instantiation here."
                sys.exit()
            ##### OPTION 5: We have both numerical and (residual) closed form information
            # Check if a numerical SS solution has been attempted and succeeded, then take solutions in here for closed form.
            elif self.switches['ss_suc'] == ['1','1'] and not numss_set.issubset(manss_set):
                if self._mesg: print "SS: CF-SS (residual) and NUM-SS information information supplied...attempting to solve SS..."
                alldic = {}
                alldic.update(self.sstate)
                alldic.update(self.paramdic)
                intup = (self.manss_sys,alldic)
                self.sssolvers.manss = ManualSteadyState(intup)
                self.sssolvers.manss.solve()
                # Here merging takes place
                self.sstate.update(self.sssolvers.manss.sstate)

        # Double check if no steady state values are negative, as LOGS may have to be taken.
        if 'sstate' in dir(self):
            for keyo in self.sstate.keys():
                if '_bar' in keyo and float(self.sstate[keyo]) < 0.0:
                    print "WARNING: Steady state value "+keyo+ " is NEGATIVE!"
                    print "This is very likely going to either error out or produce strange results"
                    print "Re-check your model declarations carefully!"
######################################### STEADY STATE CALCULATION SECTION DONE ##################################
    def init4(self):
        '''
        This model instance sub-initializor only calls the section which use the computed steady state
        in order to compute derivatives and open dynamic solver branches on the instance. But Jacobian and Hessian
        are *not* computed here, this is postponed to the next init level.
        
        .. note::
        
           The following last field is wrapped for dynamic execution:
           
             * self.sigma - the variance-covariance matrix of the iid shocks
             
           Notice that after wrapping this last field the process_queue class is instantiated at last,
           because it needs to have access to *all* of the wrapped fields. Also in this method, the function
           populate_model_stage_two() is called which prepares the nonlinear FOCs for derivative-taking.
           
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''
        txtpars = self.txtpars
        secs = txtpars.secs
        initlev = self._initlev
        ncpus = self._ncpus
        mk_hessian = self._mk_hessian
        mesg = self._mesg
        if mesg: print "INIT: Preparing DSGE model instance for computation of Jacobian and Hessian..."
        # Now populate more with stuff that needs steady state
        self = populate_model_stage_two(self, secs)

        # Need to wrap variance covariance matrix here
        self.updaters.sigma = matwrap(self,'self.sigma',initlev)
        # Need to wrap variance covariance matrix here
        self.updaters_queued.sigma = matwrap_queued(self,'self.sigma',initlev)

        ####### All queued updaters initialized, no add processing instance
        # Add the queue process instance
        self.updaters_queued.process_queue = Process_Queue(other=self)

    def init5(self):
        '''
        This model instance initialisation step is the last substantial one in which the dynamic solution of the DSGE
        model instance is finally computed using a choice of methods which can be called at runtime.
        
         
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''
        txtpars = self.txtpars
        secs = txtpars.secs
        initlev = self._initlev
        ncpus = self._ncpus
        mk_hessian = self._mk_hessian
        mesg = self._mesg
        #TODO: delay above and only import if needed
        from solvers.modsolvers import MODsolvers
        # Open the model solution tree branch
        self.modsolvers = MODsolvers()
        ######################## LINEAR METHODS !!! ############################
        # see if there are any log-linearized equations
        if any([False if 'None' in x else True for x in secs['modeq'][0]]):
            from solvers.modsolvers import PyUhlig, MatUhlig, MatKlein, MatKleinD, ForKlein
            if mesg: print "INIT: Computing DSGE model's log-linearized solution using Uhlig's Toolbox..."

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
            from solvers.modsolvers import (MatWood, ForKleinD)
            if ncpus > 1 and mk_hessian:
                if mesg: print "INIT: Computing DSGE model's Jacobian and Hessian using parallel approach..."
                self.mkjahepp()
            elif ncpus > 1 and not mk_hessian:
                if mesg: print "INIT: Computing DSGE model's Jacobian using parallel approach..."
                self.mkjahepp()
            else:
                if mesg: print "INIT: Computing DSGE model's Jacobian and Hessian using serial approach..."
                self.mkjahe()

            # Check if the obtained matrices A and B have correct dimensions
            if self.jAA.shape[0] != self.jAA.shape[1]:
                print "ERROR: Matrix A of derivatives does not have #vars=#equations"
            if self.jBB.shape[0] != self.jBB.shape[1]:
                print "ERROR: Matrix B of derivatives does not have #vars=#equations"
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
            from solvers.modsolvers import (PyKlein2D, MatKlein2D)
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
            
    def init_out(self):
        '''
        The final intializor section does some extra stuff after all has been done.
           
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''
        initlev = self._initlev
        # Compute the Eigenvalues of the AA matrix for inspection
        if 'jAA' in dir(self):
            self.mkeigv()

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

    '''
    # Get data method, in case model has not been loaded with database object
    def getdata(self,datafile=None,dbase=None):
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
    '''

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
        
        :param self: object instance
        :type self: dsge_inst
        
        :return self.numj: *(arr2d)* - attaches numerical model Jacobian to instance
        :return self.jdic: *(dic)* - attaches the analytical model Jacobian to instance
        :return self.numh: *(arr2d)* - attaches numerical model Hessian to instance
        :return self.hdic: *(dic)* - attaches the analytical 3D Hessian to instance
        :return self.jAA:  *(arr2d)* - attaches numerical AA matrix used in Forkleind solution method
        :return self.jBB:  *(arr2d)* - attaches numerical BB matrix used in Forkleind solution method
        '''
        mk_hessian = self._mk_hessian
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
            symdic[x] = sympycore.Symbol(x)
        locals().update(symdic)
#NOTE: don't do this?
# or don't modify _existing_ contents?
# see: http://docs.python.org/library/functions.html#locals
#also locals() is no longer a dictionary in Python 3
#NOTE: just use the dictionary itself
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
            elog = re.compile('LOG\(')
            while elog.search(str_tmp):
                ma = elog.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympycore.log('+str_tmp[poe:]
            eexp = re.compile('EXP\(')
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
        alldic = deepcopy(self.paramdic)
        alldic.update(self.sstate)
        for x in intup:
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            ssvar = float(alldic[vari+'_bar'])
            if 'log' in moddic[lookdic[vari]]:
                evalli.append(np.log(ssvar))
            else:
                evalli.append(ssvar)
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Make 2D symbolic and numeric Jacobian
        def mkjac(jrows=jrows,jcols=jcols):
            mesg = self._mesg
            rdic = dict([[x,'0'] for x in range(jcols)])
            jdic = dict([[x,rdic.copy()] for x in range(jrows)])
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
                if mesg:
                    # This will not work in Python 3.0
                    to = "INIT: Computing DSGE model's Jacobian: Equation ("+str(x+1)+"/"+str(jrows)+")..."
                    delete = "\b" * (len(to) + 1)
                    print "{0}{1}".format(delete, to),
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
                    numj[x,y] = eval(str(evalfo))
            # Take out the elements from the variable substitution equations
            lengor = len(self.nlsys_list)
            for keyo in jdicc.keys():
                if keyo > (lengor-1):
                    jdicc.pop(keyo)
            if mesg: print " DONE",
            return numj,jdic,jdicc,carry_over_dic

        # Now make 3D symbolic and numeric Hessian
        def mkhes(jrows=jrows,jcols=jcols,trans_dic=None):
            mesg = self._mesg
            rdic = dict([[x,'0'] for x in range(jcols)])
            rdic1 = dict([[x,rdic.copy()] for x in range(jcols)])
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
                if mesg:
                    # This will not work in Python 3.0
                    to = "INIT: Computing DSGE model's Hessian: Equation ("+str(x+1)+"/"+str(jrows)+")..."
                    delete = "\b" * (len(to) + 1)
                    print "{0}{1}".format(delete, to),
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
                        numh[count,z] = eval(str(evalfo))
                    count = count + 1
            # Take out the elements from the variable substitution equations
            lengor = len(self.nlsys_list)
            for keyo in hdicc.keys():
                if keyo > (lengor-1):
                    hdicc.pop(keyo)
            if mesg: print " DONE",
            return numh,hdic,hdicc

        self.numj,self.jdic,self.jdicc,carry_over_dic = mkjac()
        # To make line between Jacobian's and Hessian's computation
        if self._mesg: print
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
        A parallelized method using native Python and Sympy in oder
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
        
        # Import some configuration options for the DSGE model instance
        ncpus = copy.deepcopy(self._ncpus)
        mk_hessian = copy.deepcopy(self._mk_hessian)
        mesg = copy.deepcopy(self._mesg)
        
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
        nlsys = copy.deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + copy.deepcopy(self.nlsubsys)

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
            elog = re.compile('LOG\(')
            while elog.search(str_tmp):
                ma = elog.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympycore.log('+str_tmp[poe:]
            eexp = re.compile('EXP\(')
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
            list1 = list1 + [x[1] for x in self.vardic['other']['var']]
            list2 = list2 + [x for x in self.vardic['other']['mod']]
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
                evalli.append(np.log(alldic[vari+'_bar']))
            else:
                evalli.append(alldic[vari+'_bar'])
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Now make the 2D and 3D symbolic and numeric Jacobian and Hessian
        def mkjaheseq(lcount,func2,jcols,symdic,tmpli,paramdic,sstate,evaldic,mk_hessian):

            jdic = dict([[x,'0'] for x in range(jcols)])
            numj = numpy.matlib.zeros((1,jcols))
            if mk_hessian:
                rdic = dict([[x,'0'] for x in range(jcols)])
                hdic = dict([[x,rdic.copy()] for x in range(jcols)])
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
        # Support auto-detection of CPU cores
        if ncpus == 'auto':
            job_server = pp.Server(ppservers=ppservers)
            if mesg: print "INIT: Parallel execution started with "+str(job_server.get_ncpus())+ " CPU cores..."
        else:
            job_server = pp.Server(ncpus=ncpus,ppservers=ppservers)
            if mesg: print "INIT: Parallel execution started with "+str(job_server.get_ncpus())+ " CPU cores..."

        imports = ('numpy','numpy.matlib',)
        
        #job_server.submit(self, func, args=(), depfuncs=(), modules=(), callback=None, callbackargs=(), group='default', globals=None)
        # Submits function to the execution queue
           
        # func - function to be executed
        # args - tuple with arguments of the 'func'
        # depfuncs - tuple with functions which might be called from 'func'
        # modules - tuple with module names to import
        # callback - callback function which will be called with argument
        # list equal to callbackargs+(result,)
        # as soon as calculation is done
        # callbackargs - additional arguments for callback function
        # group - job group, is used when wait(group) is called to wait for
        # jobs in a given group to finish
        # globals - dictionary from which all modules, functions and classes
        # will be imported, for instance: globals=globals()
        
        jobs = [job_server.submit(mkjaheseq,(inputo,self.func2,jcols,symdic,tmpli,self.paramdic,self.sstate,evaldic,mk_hessian),(),imports) for inputo in inputs]
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
        mk_hessian = self._mk_hessian
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

                    vma = self.vreg(patup,vari0,False,'min')
                    vmavar = vma[2][1]+'_bar'
                    self.func2[i1] = sympycore.Symbol(vmavar)*(self.func2[i1])


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
        mk_hessian = self._mk_hessian
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
        mk_hessian = self._mk_hessian
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
        def mkjaheseq(lcount,func2,jcols,tmpli,paramdic,sstate,evaldic,mk_hessian,jdic,hdic):
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
        if self._ncpus == 'auto':
            job_server = pp.Server(ppservers=ppservers)
            if self._mesg: print "INIT: Parallel execution started with "+str(job_server.get_ncpus())+ " CPU cores..."
        else:
            job_server = pp.Server(ncpus=ncpus,ppservers=ppservers)
            if self._mesg: print "INIT: Parallel execution started with "+str(job_server.get_ncpus())+ " CPU cores..."
        imports = ('numpy','copy','numpy.matlib',)
        jobs = [job_server.submit(mkjaheseq,(input,self.func2,jcols,tmpli,self.paramdic,self.sstate,evaldic,mk_hessian,self.jdic,self.hdic),(),imports) for input in inputs]
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
            if ncpus > 1 and mk_hessian:
                self.mkjahepp()
            elif ncpus > 1 and not mk_hessian:
                self.mkjahepp()
            else:
                self.mkjahe()

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


