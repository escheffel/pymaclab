'''
.. module:: favar
   :platform: Linux
   :synopsis: The favar module contains the FAVAR class for estimating and doing further work with Factor-augmented Vector
              Autoregressions commonly used in applied macroeconometrics. It supports advanced methods such as bootstrapping
              confidence intervals including Killian's boostrap-after-bootstrap small-sample bias correction. Also CPU-intensive
              methods such as the bootstrap can be computed using Parallel Python to exploit multi-core CPUs. Pretty plotting methods
              are also included which depend on matplotlib.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''

from __future__ import division
import numpy
import scipy
# We don't need this anywhere and it should also be replaced with Pandas in the future
#from scikits import timeseries as ts
import matplotlib as mpl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
from matplotlib import pyplot as plt
from matplotlib import pylab as pyl
import matplotlib.dates as mdates
import datetime
import copy
from pymaclab.filters import hpfilter as hpf
import datetime
import pickle
from ..linalg import pca_module
import os
import glob
import sys
import pp
import time

# imports from refactor
from ..dattrans.transmeth import stdX, transX
from common import is_stable, genPCs, check_matrix, genEE, genLX, compbetta

class FAVAR:
############################################### MAIN INIT METHOD ############################################################
#############################################################################################################################
    def __init__(self,dates=None,rdates=None,data=None,freq='M',vnames=None,pnames=None,svnames=None,irfs=True,rescale=False,boot=True,plot=True,sfacs='auto',init=None,conf=None,mesg=False):
        self.confdic = conf
        self._mesg = mesg
        params = {'axes.labelsize':conf['graph_options']['labels']['label_fs'],
                  'text.fontsize': conf['graph_options']['labels']['text_fs'],
                  'legend.fontsize': conf['graph_options']['labels']['legend_fs'],
                  'xtick.labelsize': conf['graph_options']['labels']['xticks_fs'],
                  'ytick.labelsize': conf['graph_options']['labels']['yticks_fs'],
                  'text.usetex': conf['graph_options']['labels']['use_tex'],
                  'lines.linewidth': conf['graph_options']['lines']['width']}
        pyl.rcParams.update(params)        
        self.dates = [x.datetime for x in dates]
        for i1,elem in enumerate(rdates):
            for j1,elem2 in enumerate(elem):
                rdates[i1][j1] = elem2.datetime
        self.rdates = rdates
        self.orig_data = copy.deepcopy(data)
        self.tcode_dic = copy.deepcopy(self.confdic['tcode'])
        self.modname = self.confdic['modname']
        if conf['predel']: self.clean_dirs()
        if self._mesg: print '#####################################################'
        stepo=0
        stepo+=1
        if self._mesg: print str(stepo)+') Initialising FAVAR'
        # Populate the model with all of the information required
        self.populate(data=data,freq=freq,vnames=vnames,pnames=pnames,svnames=svnames,conf=conf)
        # Post-population diagnostics output printed to screen
        self.populate_diagnostics()
        if init == 1: return
        if self._mesg: print '#####################################################'
        stepo+=1
        if self._mesg: print str(stepo)+') Standardising data'
        # Standardizing and transforming data as needed
        self.transdata()
        if init == 2: return
        if self._mesg: print '#####################################################'
        stepo+=1
        print str(stepo)+') Estimate number of static and dynamic factors'
        self.firstpass_estimation(sfacs=sfacs)
        if init == 3: return
        # Plot actual and fitted values
        if plot and self._mesg: print 'Now plotting and saving graphs of actual, fitted and residual data for all series...'
        self.plot_all(plot=plot)      
        if self._mesg: print '#####################################################'
        stepo+=1
        if self._mesg: print str(stepo)+') Estimating coefficients of Factor VAR'
        self.estimate_fac_var()
        if init == 4: return
        if self._mesg: print '#####################################################'
        stepo+=1
        if self._mesg: print str(stepo)+') Creating companion matrices'
        if self._mesg: print 'Outputing some estimated coefficient matrices to a folder...'
        self.outmatrices()
        self.mkcompmatrices()
        if self._mesg:
            # Check if ff_betta_comp is stable
            if is_stable(compmat=self.ff_betta_comp): print 'The estimated factor VAR is stable and invertible !'
            else: print 'WARNING: The estimated factor VAR is not stable and not invertible !'
        if init == 5: return
        self.identification()
        if init == 6: return
        if init == 7: return
        if self.confdic['use_bbe_ident']:
            self.estimate_bbe_fac_var()
            self.redefine_post_bbe()
            if self.confdic['use_svnames']:
                self.scaleHmat()
        self.plot_factors_all(plot=plot)
        # Generate IRFs and bootstrapped confidence intervals and then plot the results, if all wanted via switch
        if init == 8: return
        self.do_irfs(irfs=irfs)
        if self.confdic['do_vdc']:
            vdcp = self.confdic['vdcp']
            self.genFEVD(self.bbe_resdic['fdata'],self.bbe_resdic['ffbettasm'],vdcp)
        self.writeFData()
        if self._mesg: print '*All done !*'

############################################## INIT/STEP PROCESS METHODS/CALLERS #############################################
##############################################################################################################################
    # This is the part of the FAVAR class which populates the instance with the model's properties
    def populate(self,data=None,freq=None,vnames=None,pnames=None,svnames=None,conf=None):
        # Get the short vnames and the pnames used in plots
        self.vnames = vnames
        self.pnames = pnames
        # Get the shocknames, which contain the size of the IRF shock
        self.svnames = svnames
        # Get the number of lags in the auto-regressors, both for variables equation and factor VAR
        self.nlags = conf['nlags']
        self.flags = conf['flags']
        # Set the number of factors to be used
        self.sfacs = conf['sfacs']
        # Set the criterion in percentages for convergence algorithms
        self.conv_crit = conf['conv_crit']
        # Attach some properties to the model instance
        self.data = data
        self.raw_data = copy.deepcopy(self.data)
        self.freq = freq
        self.regdata = self.data[self.nlags:,:]
        self.ncols = data.shape[1]
        self.nrows = data.shape[0]
        finit_pickled = pickle
        
    def populate_diagnostics(self):
        if self._mesg:
            print "This FAVAR model's name is: "+self.confdic['modname']
            print 'The dataset used is: '+self.confdic['datafile']
            print 'The chosen autoregressive order is: '+str(self.nlags)
            print 'The chosen factor autoregressive order is: '+str(self.flags)
            print 'The total number of variables in the system is: '+str(self.ncols)
            print 'The total number of observations before estimation is: '+str(self.nrows)
            print 'The frequency of the data is: '+str(self.freq)
            if self.confdic['use_bbe_ident']:
                print 'The model identifies structural shocks using: BBE'
            elif self.confdic['use_chol_ident']:
                print 'The model identifies structural shocks using: Simple Cholesky'
            else:
                print 'This model does not use any identification scheme for structural shocks'
     
    # This is the part of the FAVAR class which standardizes and transforms the data as needed    
    def transdata(self):
        # Transform raw data according to Stock and Watson's trans code encoding system
        if self._mesg: print 'Now TRANSFORMING data according to Stock and Watson transformation encoding...'
        self = transX(self)
        # Standardize the data ahead of using PC analysis, also store means and standard errors
        if self._mesg: print 'Now STANDARDIZING data ahead of principal component analysis...'
        self = stdX(self,standard=True)
        # The standardized data adjusted for lost observations due to lags in regression(s)
        self.regtdata = self.tdata[self.nlags:,:]
        self.orig_tdata = copy.deepcopy(self.tdata)
        # Prepare data blocks for BBE identification scheme if desired
        if self.confdic['use_bbe_ident']:
            self.prepSSF(tdata=self.tdata)

    # This is the part of the FAVAR class which estimates the static and dynamic factors
    def firstpass_estimation(self,stepo=None,sfacs=None):
        if self.confdic['xdata_filter']:
            if sfacs == 'auto':
                if self._mesg: print 'Using Bai-Ng IPc2 criterion to determine number of static factors...'
                self.findStatFacs()
            elif sfacs != 'auto':
                if self._mesg: print 'The number of static factors has been set EXOGENOUSLY to: sfacs='+str(sfacs)+'...'
                self.sfacs = copy.deepcopy(sfacs)
                self.BNcrit_res = copy.deepcopy(sfacs)
                self.estFBetta()
                self.ee_mat_l = copy.deepcopy(self.estresdic['ee_mat_l'])
                self = genLX(self)
                self.genLXF(fdata=self.fdata)
                self.BNcrit_pcmat_resdic = {}
                self.BNcrit_pcmat_resdic['auto_betta'] = copy.deepcopy(self.auto_betta)
                check_matrix(arraym=self.auto_betta,arrname='auto_betta',expdim=(self.ncols,self.nlags+self.sfacs))
                self.BNcrit_pcmat_resdic['pcmat'] = copy.deepcopy(self.fdata)
                self.BNcrit_pcmat_resdic['exp_var'] = copy.deepcopy(self.exp_var)
                self.exp_var = copy.deepcopy()
            if self.confdic['multicpu']: self.findDynFacs_pp()
            elif not self.confdic['multicpu']: self.findDynFacs()
            # Write the data to the harddrive so we can inspect the factors and the data
            self.writeFData()
            # Attach the found number of sfacs and re-run the first estimation method to get some matrices
            self.sfacs = self.dynindex[0]
            self.estFBetta()
            # Now that we have also bettax, create the filtered transformed data withou the additional factor regressors
            self.genFiltData()
            # Calculate also the dyanmic factor errors, using principal components, as suggested in Stock and Watson(2005)
            self.compFacErrors(resmat=self.dynresdic['ee_mat'])
        elif not self.confdic['xdata_filter']:
            self.estFonly(sfacs=sfacs)
    

    def estimate_fac_var(self,fdata=None,func=False):
        if fdata == None: fdata = copy.deepcopy(self.dynresdic['fdata'])
        self.genLFF(fdata=fdata)
        self.fbetta = compbetta(self,matd=fdata,nmatd=self.lffmat,smbetta=False,const=None,func=True)
        self.fbettasm = compbetta(self,matd=fdata,nmatd=self.lffmat,smbetta=True,const=None,func=True)
        self.mkfitres_ffvar(ymat=fdata,xmat=self.lffmat,betta=self.fbetta,flags=self.flags)
        self.idioresmat = copy.deepcopy(self.resmat_ff)
        if self.confdic['xdata_filter']:
            self = genLX(self,vnames=self.vnames,ydata=self.tdata,nlags=self.nlags)
            self.genLXF(vnames=self.vnames,xdic=self.lxmatdic,fdata=self.dynresdic['fdata'])
            self.mkfitres_auto(ymat=self.tdata,xdic=self.lxfmatdic,bettadic=self.auto_betta)

    def mkcompmatrices(self):
        # This creates relevant companion matrices which will be easier to work with later
        if self.confdic['xdata_filter']:
            self.mkCompAuto(bettaxm=self.bettax_mat)
        self.mkCompFF(bettafm=self.fbetta)
        self.mkCompFactor(deltam = self.bettaf_mat)
        # Create the VAR matrices for the DFM in VAR form
        self.deltaphi = numpy.dot(self.bettaf_mat,self.ff_betta_comp[:self.sfacs,:])
        if self.confdic['xdata_filter']:
            self.dlmat = self.auto_betta_comp[:self.ncols,:]
        # Stack all matrices together to one large one
        if self.confdic['xdata_filter']:
            self.mkCompAll()
            self.writeArray(inarray=self.varcomp)
            
    
    def plot_factors_all(self,plot=None):
        if plot:
            self.plotFactors()
            self.plot_scree()
            self.plotFactorResids(fresids=self.bbe_resdic['ee_mat'])        
        
    def plot_all(self,plot=None):
        if plot:
            self.plot_vals()
            self.plot_vals2()
            # Also plot the estimated factors
        
    def identification(self):
        # Identify structural shocks ahead of creating structural IRFs
        if self.confdic['use_bbe_ident']:
            # Find the appropriate number of factors spanning the space of the slow-moving variables
            if self.confdic['slow_facs'] == 'auto':
                self.slow_q = self.findStatFacs(self,statlimit=self.sfacs-2,
                                                vnames=self.confdic['slow_vars'],
                                                data=self.slow_tdata,func=True)
            elif type(self.confdic['slow_facs']) == type(2):
                print 'Number of factors for block of slow-moving data EXOGENOUSLY fixed at: '+str(self.confdic['slow_facs'])
                self.slow_q = copy.deepcopy(self.confdic['slow_facs'])
            # Given the above finding, now extract the factor space spanning the slow-moving variables and save it
            if self.confdic['xdata_filter']:
                self.pc_slow_res = genPCs(self,data=self.tdata_filtered_slow,sfacs=self.slow_q,func=True)[0]
            elif not self.confdic['xdata_filter']:
                self.pc_slow_res = genPCs(self,data=self.slow_tdata,sfacs=self.slow_q,func=True)[0]
            self.genBBEH2()
            #self.genBBEH()
        elif self.confdic['use_chol_ident']:
            self.genCholH()
        elif not self.confdic['use_bbe_ident'] and not self.confdic['use_chol_ident']:
            self.H_final = numpy.eye(self.dynindex[1])
            
    def estimate_bbe_fac_var(self):
        nmatd = self.genLFF(fdata=self.fdata_adj,flags=self.flags,func=True)
        ffbetta = compbetta(self,matd=self.fdata_adj,nmatd=nmatd,smbetta=False,const=None,func=True)
        ffbettasm = compbetta(self,matd=self.fdata_adj,nmatd=nmatd,smbetta=True,const=None,func=True)

        yfmat,resmat,covmat = self.mkfitres_ffvar(ymat=self.fdata_adj,xmat=nmatd,
                                                  betta=ffbetta,flags=self.flags,func=True)
        hmat = self.genCholH(facresmat=resmat,func=True)
        self.bbe_resdic = {}
        self.bbe_resdic['ffbetta'] = copy.deepcopy(ffbetta)
        self.bbe_resdic['ffbettasm'] = copy.deepcopy(ffbettasm)
        self.bbe_resdic['ee_mat'] = copy.deepcopy(resmat)
        self.bbe_resdic['ffvcmat'] = copy.deepcopy(covmat)
        self.bbe_resdic['ffmat'] = copy.deepcopy(yfmat)
        self.bbe_resdic['fdata'] = copy.deepcopy(self.fdata_adj)
        self.bbe_resdic['nmatd'] = copy.deepcopy(nmatd)
        self.fbetta = copy.deepcopy(self.bbe_resdic['ffbetta'])
        self.bbe_resdic['ffhmat'] = copy.deepcopy(hmat)
        self.writeArray(inarray=self.fdata_adj)

    def redefine_post_bbe(self):
        self.sfacs_adj = self.sfacs + 1
        self.sfacs = copy.deepcopy(self.sfacs_adj)
        sfcas = copy.deepcopy(self.sfacs)
        if self.confdic['xdata_filter']:
            bettax_mat = copy.deepcopy(self.bettax_mat)
        # Do not cut out the identified variable
        #indexo = self.vnames.index(self.confdic['separate_vars'][0])
        #vnames = copy.deepcopy(self.vnames)
        #vnames.remove(self.confdic['separate_vars'][0])
        self.vnames_adj = copy.deepcopy(self.vnames)
        #self.vnames = copy.deepcopy(self.vnames_adj)
        if self.confdic['xdata_filter']:
            xdic = genLX(self,vnames=vnames,ydata=self.tdata_adj,nlags=self.nlags,mlags=None,func=True)
            lxfmatdic = self.genLXF(vnames=vnames,xdic=xdic,fdata=self.fdata_adj,func=True)
            auto_betta,auto_ee_dic,ee_mat = genEE(self,vnames=vnames,ydata=self.tdata_adj,xdata=lxfmatdic,betta=None,hweightm=None,mlags=None,func=True)
            # Create the initial diagonal weighting matrix to control for heteroscedasticity ONLY
            diago = numpy.diag(numpy.cov(ee_mat))
            hweightm = numpy.zeros((len(diago),len(diago)))
            for diagi in range(0,len(diago),1):
                hweightm[diagi,diagi] = diago[diagi]
            hweightm = numpy.linalg.inv(hweightm)
            auto_betta,auto_ee_dic,ee_mat = genEE(self,vnames=vnames,ydata=self.tdata_adj,xdata=lxfmatdic,betta=None,hweightm=hweightm,mlags=None,func=True)
            # Produce matrix from auto_betta
            auto_betta_mat = auto_betta[vnames[0]]
            for varo in range(1,len(vnames),1):
                auto_betta_mat = numpy.vstack((auto_betta_mat,auto_betta[vnames[varo]]))        
            bettax_mat = auto_betta_mat[:,sfcas:]
            bettaf_mat = auto_betta_mat[:,:sfcas]
            yfmat,resmat,covmat = self.mkfitres_auto(ymat=self.tdata_adj,xdic=lxfmatdic,bettadic=auto_betta,mlags=None,func=True)
            ff_betta_comp = self.mkCompFF(bettafm=self.bbe_resdic['ffbetta'],func=True)
            self.bettax_mat = copy.deepcopy(bettax_mat)
            self.bettaf_mat = copy.deepcopy(bettaf_mat)
            self.ff_betta_comp = copy.deepcopy(ff_betta_comp)
        else:
            XX = numpy.dot(self.fdata_adj.T,self.fdata_adj)
            XXi = numpy.linalg.inv(XX)
            Xy = numpy.dot(self.fdata_adj.T,self.tdata_adj)
            bettaf_mat = numpy.dot(XXi,Xy)
            bettaf_mat = bettaf_mat.T
            self.bettaf_mat = copy.deepcopy(bettaf_mat)
            self.estimate_fac_var(fdata=self.fdata_adj)
            ff_betta_comp = self.mkCompFF(bettafm=self.fbetta,func=True)
            self.ff_betta_comp = copy.deepcopy(ff_betta_comp)
        self.H_final = copy.deepcopy(self.bbe_resdic['ffhmat'])
        self.tdata = copy.deepcopy(self.tdata_adj)
        self.dynresdic['fdata'] = self.fdata_adj
        if self.confdic['xdata_filter']:
            self.dynresdic['yfmat'] = yfmat
            self.dynresdic['lxfmatdic'] = lxfmatdic
            self.dynresdic['lxmatdic'] = xdic
            self.dynresdic['auto_betta'] = auto_betta
            self.dynresdic['ee_mat'] = resmat
            
    def scaleHmat(self,func=False):
        hmat = copy.deepcopy(self.H_final)
        svnames = self.svnames
        inames = self.confdic['separate_vars']
        linames = len(inames)
        tcoded = self.tcode_dic
        transd = self.trans_dic
        inames.reverse()
        for i1,namo in enumerate(inames):
            if namo in svnames.keys() and tcoded[namo] in [4,5,6,7,8]:
                eyem = numpy.identity(hmat.shape[0])
                eyem[-(i1+1),-(i1+1)] = (numpy.log(float(svnames[name]))/transd[namo]['std'])/float(hmat[-(i1+1),-(i1+1)])
            elif namo in svnames.keys() and tcoded[namo] not in [4,5,6,7,8]:
                eyem = numpy.identity(hmat.shape[0])       
                eyem[-(i1+1),-(i1+1)] = (float(svnames[namo])/transd[namo]['std'])/float(hmat[-(i1+1),-(i1+1)])
            hmat_unscaled = copy.deepcopy(hmat)
            hmat = numpy.dot(hmat,eyem)
        if not func:
            self.eyem = eyem
            self.H_final = copy.deepcopy(hmat)
            self.H_final_unscaled = copy.deepcopy(hmat_unscaled)
        else:
            return hmat
        
    def outmatrices(self):
        if self.confdic['xdata_filter']:
            self.writeArray(inarray=self.bettax_mat,fname='matrices/bettax_mat')
        self.writeArray(inarray=self.bettaf_mat,fname='matrices/bettaf_fmat')

            
    def do_irfs(self,irfs=None):
        if self.confdic['do_irfs']:
            # Create impulse responses using the BBE identification scheme
            if self._mesg: print 'Now generating IRs'
            if self.confdic['rescale']: self.mkphis(rescale=True)
            else: self.mkphis()
            #self.writeArray(self.phis[0],'phis_0')
            if self.confdic['bootstrap']:
                addstr = ''
                if self.confdic['multicpu']: addstr = '(parallel execution)'
                else: addstr = '(serial execution)'
                print 'Now bootstrapping IRs confidence intervals'+addstr
                t0 = time.time() 
                if not self.confdic['multicpu']:
                    if self.confdic['rescale']:
                        self.mkboot(rescale=True)
                    else:
                        self.mkboot(rescale=False)
                elif self.confdic['multicpu']:
                    if self.confdic['rescale']:
                        self.mkboot_pp(rescale=True)
                    else:
                        self.mkboot_pp(rescale=False)
            print round(time.time() - t0,2), "secs used to execute..."
            if irfs: print 'Now plotting and saving graphs of impulse responses...'
            if irfs and self.confdic['xdata_filter']: self.plot_xdata_irfs()
            if irfs: self.plot_facs_irfs()
            if irfs: self.plot_irfs()

#################################################### DATA PREPARATION/TRANSFORMATION METHODS #################################
##############################################################################################################################

        
    # Function to split the data into slow, separate and fast data blocks, for BBE identification    
    def prepSSF(self,tdata=None,func=False):
        if tdata == None:
            tdata = copy.deepcopy(self.tdata)
        slowli = self.confdic['slow_vars']
        fastli = self.confdic['fast_vars']
        sepli = self.confdic['separate_vars']
        slow_data = numpy.zeros((tdata.shape[0],len(slowli)))
        fast_data = numpy.zeros((tdata.shape[0],len(fastli)))
        sep_data = numpy.zeros((tdata.shape[0],len(sepli)))
        slowc=0
        fastc=0
        sepc=0
        for vname in self.vnames:
            if vname in slowli:
                slow_data[:,slowc] = tdata[:,self.vnames.index(vname)]
                slowc = slowc+1
            elif vname in fastli:
                fast_data[:,fastc] = tdata[:,self.vnames.index(vname)]
                fastc = fastc+1
            elif vname in sepli:
                sep_data[:,sepc] = tdata[:,self.vnames.index(vname)]
                sepc = sepc+1
        if not func:
            self.slow_tdata = copy.deepcopy(slow_data)
            self.sep_tdata = copy.deepcopy(sep_data)
            self.fast_tdata = copy.deepcopy(fast_data)
        elif func:
            return slow_data,sep_data,fast_data

    # Method used to create the filtered x data, based on the obtained xbetta_mat from FAVAR estimation    
    def genFiltData(self,data=None,bettax=None,nlags=None,func=False):
        if data == None: data = copy.deepcopy(self.tdata)
        if bettax == None: bettax = copy.deepcopy(self.bettax_mat)
        if nlags == None: nlags = copy.deepcopy(self.nlags)
        vnames = copy.deepcopy(self.vnames)
        xmatdic = genLX(self,vnames=vnames,ydata=data,nlags=nlags,mlags=None,func=True)
        # Now generate the corresponding regressor array
        filt_data = numpy.zeros((data.shape[0]-nlags,data.shape[1]))
        for i1,namo in enumerate(vnames):
            filt_data[:,i1] = data[nlags:,i1] - numpy.dot(xmatdic[namo],bettax[i1,:].T)
        # Also create the same datamats for slow, separate and fast data
        filt_data_slow = numpy.zeros((filt_data.shape[0],len(self.confdic['slow_vars'])))
        for i1,namo in enumerate(self.confdic['slow_vars']):
            filt_data_slow[:,i1] = copy.deepcopy(filt_data[:,vnames.index(namo)])
        filt_data_separate = numpy.zeros((filt_data.shape[0],len(self.confdic['separate_vars'])))
        filt_data_separate = copy.deepcopy(filt_data[:,vnames.index(self.confdic['separate_vars'][0])])
        if not func:
            self.tdata_filtered = copy.deepcopy(filt_data)
            self.tdata_filtered_slow = copy.deepcopy(filt_data_slow)
            self.tdata_filtered_separate = copy.deepcopy(filt_data_separate)
        elif func:
            return filt_data
        
############################################ MODEL (SHOCKS) IDENTIFICATION-RELATED METHODS ###################################
##############################################################################################################################
    def genBBEH2(self):
        # Number of factors for block of slow-moving data
        slow_q = copy.deepcopy(self.slow_q)
        # Number of factors present in all of the data
        all_q = copy.deepcopy(self.dynindex[1])
        # Matrix(array) of the extracted factors from the block of slow-moving data
        slow_facs = copy.deepcopy(self.pc_slow_res)
        # Matrix(array) of the extracted factors from all of the data
        all_facs = copy.deepcopy(self.dynresdic['fdata'])
        # This contains the filtered series of the separate variable whose shock we want to study
        if self.confdic['xdata_filter']:
            shockvar_res = copy.deepcopy(self.tdata_filtered_separate)
        elif not self.confdic['xdata_filter']:
            shockvar_res = copy.deepcopy(self.sep_tdata)
        XX = numpy.hstack((slow_facs,numpy.reshape(shockvar_res,(slow_facs.shape[0],1))))
        shock_loading = numpy.zeros((all_q,1))
        for itero in xrange(0,all_q):
            betta_temp = numpy.dot(numpy.linalg.inv(numpy.dot(XX.T,XX)),numpy.dot(XX.T,all_facs[:,itero]))
            shock_loading[itero,0] = betta_temp[-1]
        redux_space_mat = all_facs - numpy.dot(numpy.reshape(shockvar_res,(len(shockvar_res),1)),shock_loading.T)
        self.shock_loading = shock_loading
        self.redux_space_mat = redux_space_mat
        #indexo = self.vnames.index(self.confdic['separate_vars'][0])
        self.fdata_adj = numpy.hstack((redux_space_mat,numpy.reshape(shockvar_res,(len(shockvar_res),1))))
        #self.tdata_adj = numpy.hstack((self.tdata[:,:indexo],self.tdata[:,indexo+1:]))
        self.tdata_adj = copy.deepcopy(self.tdata)
        

    # Funtion to create the H matrix needed to identify structural shocks using BBE identification strategy    
    def genBBEH(self):
        slow_q = copy.deepcopy(self.slow_q)
        idioresmat = self.idioresmat
        facresmat = self.resmat_ff
        resmat = self.resmat_auto
        ee_ss_mat = numpy.zeros((resmat.shape[0],self.slow_tdata.shape[1]))
        min_row = min(facresmat.shape[0],ee_ss_mat.shape[0])
        posi=0
        for vname in self.confdic['slow_vars']:
            ee_ss_mat[:,posi] = resmat[:,self.vnames.index(vname)]
            posi=posi+1
        ee_ff_mat = numpy.zeros((resmat.shape[0],self.fast_tdata.shape[1]))
        posi=0
        for vname in self.confdic['fast_vars']:
            ee_ff_mat[:,posi] = resmat[:,self.vnames.index(vname)]
            posi=posi+1
        # Run regression of resmat on the factor shocks, then take SVD
        XX = numpy.dot(facresmat.T,facresmat)
        XXi = numpy.linalg.inv(XX)
        indi_ee = ee_ss_mat.shape[0] - min_row
        indi_fac = facresmat.shape[0] - min_row
        Xy = numpy.dot(facresmat[indi_fac:,:].T,ee_ss_mat[indi_ee:,:])
        bettass = numpy.dot(XXi,Xy)
        bettass = bettass.T
        resids = ee_ss_mat[indi_ee:,:] - numpy.dot(facresmat[indi_fac:,:],bettass.T)
        q,r = numpy.linalg.qr(facresmat)
        u,s,v = numpy.linalg.svd(numpy.dot(q.T,ee_ss_mat[indi_ee:,:]))
        sss = numpy.zeros((len(s),len(s)))
        for diago in range(0,len(s),1):
            sss[diago,diago] = s[diago]
        inner_second = numpy.dot(sss[:slow_q,:slow_q],v[:,:slow_q].T)
        inner = numpy.dot(u[:,:slow_q],inner_second)
        betta_final = numpy.dot(numpy.linalg.inv(r),inner)
        betta_final = betta_final.T
        nu,ns,nv = numpy.linalg.svd(betta_final)
        # We want the QL and not the QR factorisation, so flip original matrix on two axes before factoring
        ql,ll = numpy.linalg.qr(numpy.flipud(numpy.fliplr(betta_final)))
        qr,rr = numpy.linalg.qr(betta_final)
        # After getting QR, flip around again to get QL of original matrix
        self.ql = numpy.flipud(numpy.fliplr(ql[:,:slow_q]))
        self.ll = numpy.flipud(numpy.fliplr(ll[:slow_q,:]))
        self.qr = qr[:,:slow_q]
        self.rr = rr[:slow_q,:]
        Hs = self.rr
        # Check if this has worked
        self.testql = numpy.dot(self.ql,self.ll)
        self.testqr = numpy.dot(self.qr,self.rr)
        self.betta_final = betta_final
        self.bettass = bettass
        zeta_ss = numpy.dot(self.rr,facresmat.T).T
        ee_rr_mat = resmat[:,self.vnames.index(self.confdic['separate_vars'][0])]
        XX_fac = numpy.dot(facresmat.T,facresmat)
        XXi_fac = numpy.linalg.inv(XX_fac)
        Xz = numpy.dot(facresmat[indi_fac:,:].T,ee_rr_mat[indi_ee:])
        ZZ_ss = numpy.dot(zeta_ss.T,zeta_ss)
        ZZi_ss = numpy.linalg.inv(ZZ_ss)
        Zz = numpy.dot(zeta_ss.T,ee_rr_mat[indi_ee:])
        zeta_rr = numpy.dot(numpy.dot(XXi_fac,Xz).T,facresmat.T).T-numpy.dot(numpy.dot(ZZi_ss,Zz).T,zeta_ss.T).T
        zeta_rr = numpy.reshape(zeta_rr,(len(zeta_rr),1))
        Xzr = numpy.dot(facresmat.T,zeta_rr)
        Hr = numpy.dot(XXi_fac,Xzr).T
        ZZ_rr = numpy.dot(zeta_rr.T,zeta_rr)
        ZZi_rr = numpy.linalg.inv(ZZ_rr)
        Xf = numpy.dot(facresmat[indi_fac:,:].T,ee_ff_mat[indi_ee:,:])
        # zeta_ff = numpy.dot(numpy.dot(XXi_fac,Xf).T,facresmat.T).T-numpy.dot(numpy.dot(ZZi_ss,numpy.dot(zeta_ss.T,ee_ff_mat)).T,zeta_ss.T).T-numpy.dot(numpy.dot(ZZi_rr,numpy.dot(zeta_rr.T,ee_ff_mat)).T,zeta_rr.T).T
        # Xzf = numpy.dot(facresmat.T,zeta_ff)
        # Hf = numpy.dot(XXi_fac,Xzf).T
        Hf = numpy.cov(facresmat.T)[slow_q+1:,:]
        H_final = numpy.vstack((Hs,Hr,Hf))
        self.H_final = copy.deepcopy(H_final)

        
    # Funtion to create the H matrix needed to identify structural shocks using BBE identification strategy    
    def genCholH(self,facresmat=None,func=False):
        idioresmat = self.idioresmat
        if facresmat == None: facresmat = copy.deepcopy(self.resmat_ff)
        resmat = copy.deepcopy(self.dynresdic['ee_mat'])
        # Create the variance-covariance matrix of the residuals
        vcmat = numpy.cov(facresmat.T)
        # Create the Cholesky decomposition of vcmat
        cholmat = numpy.linalg.cholesky(vcmat)
        if not func:
            self.H_final = copy.deepcopy(cholmat)
        elif func:
            return copy.deepcopy(cholmat)
   
########################################### PRINCIPAL COMPONENTS ANALYSIS METHODS ###########################################
#############################################################################################################################       


        
################################################# VARIOUS REGRESSOR-MATRIX GENERATOR METHODS #################################
##############################################################################################################################

        
    # Function to create regressor matrix containing ONLY OWN lags AND static factors,
    # where factor coefficients are ORDERED FIRST
    def genLXF(self,vnames=None,xdic=None,fdata=None,func=False):
        # Get all the needed variables from instance
        if xdic == None:
            xdic = copy.deepcopy(self.lxmatdic)
        min_row = min(xdic[xdic.keys()[0]].shape[0],fdata.shape[0])
        if vnames == None:
            vnames = copy.deepcopy(self.vnames)
        cols = len(xdic.keys())
        indi_nmatd = xdic[xdic.keys()[0]].shape[0]-min_row
        indi_fdata = fdata.shape[0]-min_row
        lxfmatdic = {}
        for col in range(0,cols,1):
            nmatd = xdic[vnames[col]]
            nmatd = numpy.hstack((fdata[indi_fdata:,:],nmatd[indi_nmatd:,:]))
            lxfmatdic[vnames[col]] = copy.deepcopy(nmatd)
        if not func:
            self.lxfmatdic = copy.deepcopy(lxfmatdic)
        elif func:
            return lxfmatdic

        
    # Function to create a matrix of coefficients on the factor auto-regressors for each factor
    def genLFF(self,fdata=None,flags=None,func=False):
        # Get all the needed variables from instance
        if flags == None: flags = copy.deepcopy(self.flags)
        matd = fdata
        rows = matd.shape[0]
        cols = matd.shape[1]
        nmatd = numpy.zeros([rows-flags,cols*flags])
        for lag in range(0,flags,1):
            for col in range(0,cols,1):
                shift = lag*(cols-1)
                if lag == 0:
                    nmatd[:,col+shift+lag] = matd[lag:-(flags-(lag)),col]
                elif lag == flags - 1:
                    nmatd[:,col+shift+lag] = matd[lag:-1,col]
                elif lag > 0:
                    nmatd[:,col+shift+lag] = matd[lag:-(flags-(lag)),col]
        # Flip matrix left to right in order to have p(1) lags first, then re-order variables
        nmatd = numpy.fliplr(nmatd)
        for elem in range(0,flags,1):
            nmatd[:,elem*cols:(elem+1)*cols] = numpy.fliplr(nmatd[:,elem*cols:(elem+1)*cols])
        if not func:
            self.lffmat = copy.deepcopy(nmatd)
        elif func:
            return nmatd
        
    # Function to create a lagmat (w or wo the current value ordered first) with flexible lag structure
    # This is used to create flexible fdata based on larr matrix
    def genLFS(self,data=None,wcurr=True,larr=None,func=False):
        countmat = numpy.zeros((larr.shape[0],larr.shape[1]))
        for row in range(0,countmat.shape[0],1):
            for col in range(0,countmat.shape[1],1):
                if larr[row,col] > 0.0: countmat[row,col] = 1.0 
        sumarr = numpy.sum(countmat,axis=1)
        if not wcurr: maxlag = larr.shape[1]
        elif wcurr: maxlag = larr.shape[1]-1
        totcols = int(sum(sumarr))
        colvars = larr.shape[0]
        nmatd = numpy.zeros((data.shape[0]-maxlag,1))
        frows = data.shape[0]-maxlag
        tindex = data[maxlag:,:].shape[0]
        if wcurr:
            for elem in range(0,colvars,1):
                if larr[elem,0] > 0.0: nmatd = numpy.hstack((nmatd,numpy.reshape(data[maxlag:,elem],(data[maxlag:,elem].shape[0],1))))
        # Get rid of first column which was used to stack on
        nmatd = nmatd[:,1:]
        if wcurr: cmat = range(1,larr.shape[1],1)
        elif not wcurr: cmat = range(0,larr.shape[1],1)
        for i1,lagos in enumerate(cmat):
            for i2,elem in enumerate(larr[:,lagos]):
                if elem > 0.0: nmatd = numpy.hstack((nmatd,numpy.reshape(data[maxlag-lagos:(tindex+maxlag)-lagos,i2],(data[maxlag-lagos:(tindex+maxlag)-lagos,i2].shape[0],1))))
        if not wcurr: nmatd = nmatd[:,1:]
        if not func:
            self.nfdatam = nmatd
        elif func:
            return nmatd 


    # Function for creating residuals, covmatrix and fitted values, given betta and xmat for auto_reg
    def mkfitres_auto(self,ymat=None,xdic=None,bettadic=None,mlags=None,func=False):
        # Create fitted values and residuals
        vnames = self.vnames
        if bettadic == None: auto_betta = copy.deepcopy(self.auto_betta)
        else: auto_betta = bettadic
        min_row=min(xdic[xdic.keys()[0]].shape[0],ymat.shape[0])
        yindi = ymat.shape[0] - min_row
        xindi = xdic[xdic.keys()[0]].shape[0] - min_row
        yfmat = numpy.zeros((ymat[yindi:].shape[0],ymat.shape[1]))
        for i1,vname in enumerate(vnames):
            yfmat[:,i1] = numpy.dot(xdic[vname][xindi:,:],auto_betta[vname].T)
        resmat = ymat[yindi:]-yfmat
        covmat = numpy.cov(resmat.T)
        if not func:
            self.yfmat_auto = copy.deepcopy(yfmat)
            self.resmat_auto = copy.deepcopy(resmat)
            self.covmat_auto = copy.deepcopy(covmat)
        elif func:
            return yfmat,resmat,covmat
 
   
    # Function for creating residuals, covmatrix and fitted values, given betta and xmat for factor VAR
    def mkfitres_ffvar(self,ymat=None,xmat=None,betta=None,flags=None,func=False):
        # Create fitted values and residuals
        yfmat = numpy.dot(xmat,betta.T)
        resmat = ymat[flags:]-yfmat
        covmat = numpy.cov(resmat.T)
        if not func:
            self.yfmat_ff = copy.deepcopy(yfmat)
            self.resmat_ff = copy.deepcopy(resmat)
            self.covmat_ff = copy.deepcopy(covmat)
        elif func:
            return yfmat,resmat,covmat



################################################## FAVAR ESTIMATION-RELATED METHODS #########################################
#############################################################################################################################
    # Function to implement Stock and Watson's 2005 iterative Cochrane-Orcutt type estimation algorithm
    # But this is done here for GIVEN static factors and GIVEN lags on the data
    def estFBetta(self,vnames=None,data=None,func=False):
        nlags = self.nlags
        flags = self.flags
        if data == None:
            cols = self.ncols
        else:
            cols = data.shape[1]
        if vnames == None:
            vnames = self.vnames
        sfcas = self.sfacs
        flags = self.flags
        ee_crit = copy.deepcopy(self.confdic['conv_crit']) #In percent !
        if data == None:
            data = copy.deepcopy(self.tdata)
        # THIS FIRST SET OF COMMANDS CREATES THE FIRST SET OF RESULTS BEFORE ENTERING THE LOOP #
        # Create the first estimate of the bettas and also the residuals without the factors
        lxmatdic = genLX(self,vnames=vnames,ydata=data,mlags=max(nlags,flags),func=True)
        auto_betta,auto_ee_dic,ee_mat_l = genEE(self,vnames=vnames,ydata=data,xdata=lxmatdic,betta=None,hweightm=None,mlags=nlags,func=True)
        # Extract the factors from the first set of residuals
        pca_result = genPCs(self,data=ee_mat_l,func=True)
        pcmat = pca_result[0]
        exp_var = pca_result[1]
        lxfmatdic = self.genLXF(vnames=vnames,xdic=lxmatdic,fdata=pcmat,func=True)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,vnames=vnames,ydata=data,xdata=lxfmatdic,betta=None,hweightm=None,mlags=nlags,func=True)
        # Create the initial diagonal weighting matrix to control for heteroscedasticity ONLY
        diago = numpy.diag(numpy.cov(ee_mat_s))
        hweightm = numpy.zeros((len(diago),len(diago)))
        for diagi in range(0,len(diago),1):
            hweightm[diagi,diagi] = diago[diagi]
        hweightm = numpy.linalg.inv(hweightm)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,vnames=vnames,ydata=data,xdata=lxfmatdic,betta=None,hweightm=hweightm,mlags=nlags,func=True)
        # Save the initial average variance of residuals for calculation of convergence criterion value
        suma_start = copy.deepcopy(numpy.average(numpy.var(ee_mat_s,axis=0)))
        # Also copy the current sum of residuals into a variable called suma
        suma = copy.deepcopy(suma_start)
        ee_now = 100.0
        # END OF INITIALISATION LOOP #
        iter = 0
        # I take the absolute value here because in the beginning the new RSS might be bigger than the initial
        while abs(ee_now) > ee_crit:
            iter=iter+1
            # Strip out the estimated coefficients for the factors
            bettax = {}
            bettaf = {}
            for vname in vnames:
                bettax[vname] = auto_betta[vname][sfcas:]
                bettaf[vname] = auto_betta[vname][:sfcas]  
            # Get the new residuals but this time with updated regression coefficients on the lagged Xs
            lxmatdic = genLX(self,vnames=vnames,ydata=data,mlags=nlags,func=True)
            auto_betta,auto_ee_dic,ee_mat_l = genEE(self,vnames=vnames,ydata=data,xdata=lxmatdic,betta=bettax,hweightm=None,mlags=nlags,func=True)
            # Extract the factors from the current residuals
            pca_result = genPCs(self,data=ee_mat_l,func=True)
            pcmat = pca_result[0]
            exp_var = pca_result[1]
            # Build the new set of regressors with the factors, run regression and save residuals
            lxfmatdic = self.genLXF(vnames=vnames,xdic=lxmatdic,fdata=pcmat,func=True)
            auto_betta,auto_ee_dic,ee_mat_s = genEE(self,vnames=vnames,ydata=data,xdata=lxfmatdic,betta=None,hweightm=None,mlags=nlags,func=True)
            # Create the updated diagonal weighting matrix to control for heteroscedasticity ONLY
            diago = numpy.diag(numpy.cov(ee_mat_s))
            hweightm = numpy.zeros((len(diago),len(diago)))
            for diagi in range(0,len(diago),1):
                hweightm[diagi,diagi] = diago[diagi]
            hweightm = numpy.linalg.inv(hweightm)
            auto_betta,auto_ee_dic,ee_mat_s = genEE(self,vnames=vnames,ydata=data,xdata=lxfmatdic,betta=None,hweightm=hweightm,mlags=nlags,func=True)
            # Save the sum of residuals to compare with the next round
            sume = numpy.average(numpy.var(ee_mat_s,axis=0))
            ee_now=((suma-sume)/suma_start)*100.0
            # Do not allow negative numbers to converge or make progress, as RSS needs to be non-increasing
            if str(ee_now)[0] == '-':
                ee_now=100.0
            suma = copy.deepcopy(sume)         
        # Produce matrix from auto_betta
        auto_betta_mat = auto_betta[vnames[0]]
        for varo in range(1,cols,1):
            auto_betta_mat = numpy.vstack((auto_betta_mat,auto_betta[vnames[varo]]))
        # Copy the final results into a dictionary
        estresdic = {}
        estresdic['auto_betta'] = copy.deepcopy(auto_betta)
        estresdic['auto_ee_dic'] = copy.deepcopy(auto_ee_dic)
        estresdic['ee_mat'] = copy.deepcopy(ee_mat_s)
        estresdic['pcmat'] = copy.deepcopy(pcmat)
        estresdic['exp_var'] = copy.deepcopy(exp_var)
        estresdic['hweightm'] = copy.deepcopy(hweightm)
        estresdic['conv_crit'] = copy.deepcopy(ee_crit)
        estresdic['last_crit'] = copy.deepcopy(ee_now)
        estresdic['max_iters'] = copy.deepcopy(iter)
        estresdic['ee_mat_l'] = copy.deepcopy(ee_mat_l)
        if not func:
            self.auto_betta_mat = copy.deepcopy(auto_betta_mat)
            self.bettax_mat = copy.deepcopy(self.auto_betta_mat[:,sfcas:])
            self.bettaf_mat = copy.deepcopy(self.auto_betta_mat[:,:sfcas])
            self.auto_betta = copy.deepcopy(estresdic['auto_betta'])
            self.auto_ee_mat = copy.deepcopy(estresdic['ee_mat'])
            self.fdata = estresdic['pcmat']
            self.exp_var = estresdic['exp_var']
            self.hweightm = estresdic['hweightm']
            self.regfdata = estresdic['pcmat'][:,:self.sfacs]
            self.estresdic = copy.deepcopy(estresdic)
        elif func:
            return estresdic
        
    def estFonly(self,data=None,vnames=None,sfacs=None,func=False):
        flags = self.flags
        if data == None:
            cols = self.ncols
        else:
            cols = data.shape[1]
        if vnames == None:
            vnames = self.vnames
        if sfacs == None: sfcas = self.sfacs
        flags = self.flags
        if data == None:
            data = copy.deepcopy(self.tdata)
        # Extract the factors from the first set of residuals
        pca_result = genPCs(self,data=data,func=True)
        pcmat = pca_result[0]
        exp_var = pca_result[1]
        self.exp_var = copy.deepcopy(exp_var)
        self.fdata = copy.deepcopy(pcmat)
        XX = numpy.dot(pcmat.T,pcmat)
        XXi = numpy.linalg.inv(XX)
        Xy = numpy.dot(pcmat.T,data)
        betta = numpy.dot(XXi,Xy)
        betta = betta.T
        self.bettaf_mat = copy.deepcopy(betta)
        self.dynresdic = {}
        self.dynresdic['fdata'] = copy.deepcopy(pcmat)
        self.dynresdic['yfmat'] = copy.deepcopy(numpy.dot(betta,pcmat.T).T)
        self.dynresdic['ee_mat'] = copy.deepcopy(data - self.dynresdic['yfmat'])
        self.dynindex = []
        self.dynindex.append(pcmat.shape[1])
        self.dynindex.append(pcmat.shape[1])
        self.sfacs = pcmat.shape[1]
    
    # This is the IC(p2) criterion from Ng and Bai 2002 Econometrica    
    def compBNcrit(self,resmat=None,nobs=None,tnobs=None,nfacs=None,func=False):
        avgvar = numpy.average(numpy.var(resmat,axis=0))
        Cnt = min(numpy.sqrt(nobs),numpy.sqrt(tnobs))
        ICp2 = numpy.log(avgvar)+nfacs*((nobs+tnobs)/(nobs*tnobs))*numpy.log(Cnt**2)
        if not func:
            self.icp2_crit = copy.deepcopy(ICp2)
        elif func:
            return ICp2
     
    # Function to find the static factors using the NB criterion   
    def findStatFacs(self,crit=None,statlimit=None,data=None,vnames=None,func=False):
        limswitch = False
        if crit != None:
            orig_crit = copy.deepcopy(self.conv_crit)
            self.conv_crit = crit
        if statlimit != None:
            statend = statlimit
        cols = self.ncols
        if statlimit == None:
            statend = cols
        if vnames == None:
            vnames = self.vnames
        if data != None:
            rows = data.shape[0]
            cols = data.shape[1]
        elif data == None:
            rows = self.nrows
            cols = self.ncols
        orig_sfacs = copy.deepcopy(self.sfacs)
        sfac_ops = [x for x in range(1,statend+1,1)]
        ipc2_crit = []
        estres_li = []
        for sfaco in sfac_ops:
            print 'Now estimating model for '+str(sfaco)+' factor representation...'
            self.sfacs = copy.deepcopy(sfaco)
            if data == None:
                estres_li.append(self.estFBetta(vnames=vnames,func=True))
            else:
                estres_li.append(self.estFBetta(vnames=vnames,data=data,func=True))
            # Keep length of recursive storage for results efficiently short
            if len(estres_li) == 3: estres_li = estres_li[1:]
            ipc2_crit.append(self.compBNcrit(resmat=estres_li[-1]['ee_mat'],nobs=cols,tnobs=rows,nfacs=self.sfacs,func=True))
            print 'Found ICp2 criterion at value of: '+str(ipc2_crit[-1])
            if len(ipc2_crit) > 1 and ipc2_crit[-2] < ipc2_crit[-1]: break
            elif sfaco == statend:
                limswitch=True
                break
        if ipc2_crit.index(min(ipc2_crit))+1 == 1:
                print 'WARNING: The was only ONE static factor in the model, please adjust model parameters !!!'
        print 'The model appears to be best represented using: '+str(ipc2_crit.index(min(ipc2_crit))+1)+' STATIC factors.'
        if not func:
            self.sfacs = copy.deepcopy(orig_sfacs)
            if crit != None: self.conv_crit = copy.deepcopy(orig_crit)
            self.BNcrit_res = copy.deepcopy(ipc2_crit.index(min(ipc2_crit))+1)
            if not limswitch:
                self.BNcrit_pcmat_resdic = copy.deepcopy(estres_li[0])
                # Copy also the final result of the unfactored residuals, we need this later...
                self.ee_mat_l = copy.deepcopy(estres_li[0]['ee_mat_l'])
            elif limswitch:
                self.ee_mat_l = copy.deepcopy(estres_li[1]['ee_mat_l'])
                self.BNcrit_pcmat_resdic = copy.deepcopy(estres_li[1])
        elif func:
            self.sfacs = copy.deepcopy(orig_sfacs)
            if crit != None: self.conv_crit = copy.deepcopy(orig_crit)
            return copy.deepcopy(ipc2_crit.index(min(ipc2_crit))+1)

    # Function to find the dynamic factors using the NB criterion
    def findDynFacs(self,sfacs=None,crit=None,func=None):
        cols = self.ncols
        nlags = self.nlags
        if sfacs == None:
            statfacs = copy.deepcopy(self.BNcrit_res)
        else:
            statfacs = sfacs
        data = copy.deepcopy(self.tdata)
        auto_betta = copy.deepcopy(self.BNcrit_pcmat_resdic['auto_betta'])
        pcmat = copy.deepcopy(self.BNcrit_pcmat_resdic['pcmat'])
        exp_var = copy.deepcopy(self.BNcrit_pcmat_resdic['exp_var'])
        vnames = self.vnames
        if crit != None:
            orig_crit = copy.deepcopy(self.conv_crit)
            self.conv_crit = crit
        orig_sfacs = copy.deepcopy(self.sfacs)
        orig_flags = copy.deepcopy(self.flags)
        critmat = numpy.zeros((statfacs,statfacs))
        ssemat = numpy.zeros((statfacs,statfacs))
        fitmat = {}
        for row in range(0,critmat.shape[0],1):
            for col in range(0,critmat.shape[1],1):
                critmat[row,col] = numpy.nan
        for stato in range(1,statfacs+1,1):
            fitmat[str(stato)]={}
            # print 'Now computing combinations for '+str(stato)+' static factors...'
            for dyno in range(1,stato+1,1):
                min_lagos = int(numpy.trunc(stato/dyno))
                rest_lagos = stato-dyno*min_lagos
                if rest_lagos == 0:
                    larr = numpy.ones((dyno,min_lagos))
                    end_lagos = min_lagos
                elif rest_lagos > 0:
                    larr = numpy.ones((dyno,min_lagos+1))
                    larr[:,-1] = 0.0
                    end_lagos = min_lagos+1
                    for elem in range(0,rest_lagos,1):
                        larr[elem,-1] = 1.0
                lxmatdic = genLX(self,ydata=data,nlags=nlags,mlags=nlags,func=True)
                if dyno == stato+1:
                    fdata = self.genLFS(data=pcmat,wcurr=True,larr=larr,func=True)
                else:
                    fdata = self.genLFS(data=pcmat[:,:dyno],wcurr=True,larr=larr,func=True)
                lxfmatdic = self.genLXF(xdic=lxmatdic,fdata=fdata,func=True)
                auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=None,mlags=max(nlags,end_lagos+rest_lagos),func=True)
                # Create the updated diagonal weighting matrix to control for heteroscedasticity ONLY
                diago = numpy.diag(numpy.cov(ee_mat_s))
                hweightm = numpy.zeros((len(diago),len(diago)))
                for diagi in range(0,len(diago),1):
                    hweightm[diagi,diagi] = diago[diagi]
                hweightm = numpy.linalg.inv(hweightm)
                auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=hweightm,mlags=max(nlags,end_lagos),func=True)
                yfmat,resmat,covmat = self.mkfitres_auto(ymat=data,xdic=lxfmatdic,bettadic=auto_betta,mlags=max(nlags,end_lagos),func=True)
                critmat[dyno-1,stato-1] = self.compBNcrit(resmat=ee_mat_s,nobs=self.ncols,tnobs=self.nrows,nfacs=statfacs,func=True)
                ssemat[dyno-1,stato-1] = numpy.average(numpy.std(ee_mat_s,axis=0))
                ssemat[dyno-1,stato-1] = numpy.round(ssemat[dyno-1,stato-1],3)
                fitmat[str(stato)][str(dyno)] = yfmat
        self.dynmat = copy.deepcopy(critmat)
        self.dynmat_round = numpy.round(self.dynmat,3)
        self.ssemat = copy.deepcopy(ssemat)
        self.fitmat = copy.deepcopy(fitmat)
        minval = numpy.nanmin(self.dynmat)
        minindi = []
        for col in range(0,self.dynmat.shape[1],1):
            for row in range(0,self.dynmat.shape[0],1):
                if self.dynmat[row,col] == minval:
                    minindi.append(row+1)
                    minindi.append(col+1)
        self.dynindex = copy.deepcopy(minindi)
        # One more round with the found r and q to get the results
        self.dynresdic = {}
        stato = self.dynindex[0]
        dyno = self.dynindex[1]
        min_lagos = int(numpy.trunc(stato/dyno))
        rest_lagos = stato-dyno*min_lagos
        if rest_lagos == 0:
            larr = numpy.ones((dyno,min_lagos))
            end_lagos = min_lagos
        elif rest_lagos > 0:
            larr = numpy.ones((dyno,min_lagos+1))
            larr[:,-1] = 0.0
            end_lagos = min_lagos+1
            for elem in range(0,rest_lagos,1):
                larr[elem,-1] = 1.0
        lxmatdic = genLX(self,ydata=data,nlags=nlags,mlags=nlags,func=True)
        self.dynresdic['lxmatdic'] = copy.deepcopy(lxmatdic)
        self.dynresdic['exp_var'] = copy.deepcopy(exp_var)
        if dyno == stato+1:
            fdata = self.genLFS(data=pcmat,wcurr=True,larr=larr,func=True)
        else:
            fdata = self.genLFS(data=pcmat[:,:dyno],wcurr=True,larr=larr,func=True)
        self.dynresdic['fdata'] = copy.deepcopy(fdata)
        lxfmatdic = self.genLXF(xdic=lxmatdic,fdata=fdata,func=True)
        self.dynresdic['lxfmatdic'] = copy.deepcopy(lxfmatdic)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=None,mlags=max(nlags,end_lagos+rest_lagos),func=True)
        # Create the updated diagonal weighting matrix to control for heteroscedasticity ONLY
        diago = numpy.diag(numpy.cov(ee_mat_s))
        hweightm = numpy.zeros((len(diago),len(diago)))
        for diagi in range(0,len(diago),1):
            hweightm[diagi,diagi] = diago[diagi]
        hweightm = numpy.linalg.inv(hweightm)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=hweightm,mlags=max(nlags,end_lagos),func=True)
        self.dynresdic['auto_betta'] = copy.deepcopy(auto_betta)
        self.dynresdic['ee_mat'] = copy.deepcopy(ee_mat_s)
        self.dynresdic['ee_mat_l'] = copy.deepcopy(ee_ma)
        yfmat,resmat,covmat = self.mkfitres_auto(ymat=data,xdic=lxfmatdic,bettadic=auto_betta,mlags=max(nlags,end_lagos),func=True)
        self.dynresdic['yfmat'] = copy.deepcopy(yfmat)
        self.dynresdic['covmat'] = copy.deepcopy(covmat)
        print 'The model appears to be best represented by '+str(self.dynindex[0])+' static and '+str(self.dynindex[1])+' dynamic factors'
        if self.dynindex[0] == self.dynindex[1]:
            print 'Since in this model q=r, the model is best represented as a static factor model with '+str(self.dynindex[0])+' static factors'


    def refac_pp2(self,stato,dyno,data,nlags,pcmat,statfacs):
        min_lagos = int(numpy.trunc(stato/dyno))
        rest_lagos = stato-dyno*min_lagos
        if rest_lagos == 0:
            larr = numpy.ones((dyno,min_lagos))
            end_lagos = min_lagos
        elif rest_lagos > 0:
            larr = numpy.ones((dyno,min_lagos+1))
            larr[:,-1] = 0.0
            end_lagos = min_lagos+1
            for elem in range(0,rest_lagos,1):
                larr[elem,-1] = 1.0
        lxmatdic = genLX(self,ydata=data,nlags=nlags,mlags=nlags,func=True)
        if dyno == stato+1:
            fdata = self.genLFS(data=pcmat,wcurr=True,larr=larr,func=True)
        else:
            fdata = self.genLFS(data=pcmat[:,:dyno],wcurr=True,larr=larr,func=True)
        lxfmatdic = self.genLXF(xdic=lxmatdic,fdata=fdata,func=True)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=None,mlags=max(nlags,end_lagos+rest_lagos),func=True)
        # Create the updated diagonal weighting matrix to control for heteroscedasticity ONLY
        diago = numpy.diag(numpy.cov(ee_mat_s))
        hweightm = numpy.zeros((len(diago),len(diago)))
        for diagi in range(0,len(diago),1):
            hweightm[diagi,diagi] = diago[diagi]
        hweightm = numpy.linalg.inv(hweightm)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=hweightm,mlags=max(nlags,end_lagos),func=True)
        yfmat,resmat,covmat = self.mkfitres_auto(ymat=data,xdic=lxfmatdic,bettadic=auto_betta,mlags=max(nlags,end_lagos),func=True)
        critval = self.compBNcrit(resmat=ee_mat_s,nobs=self.ncols,tnobs=self.nrows,nfacs=statfacs,func=True)
        return critval


    # Function to find the dynamic factors using the NB criterion
    def findDynFacs_pp(self,sfacs=None,crit=None,func=None):
        # Configure and start the parallel Python jobserver
        ppservers = ()
        if self.confdic['ncpus'] == 'auto':
            job_server = pp.Server(ncpus='autodetect',ppservers=ppservers)
        else:
            job_server = pp.Server(ncpus=self.confdic['ncpus'],ppservers=ppservers)
        cols = self.ncols
        nlags = self.nlags
        if sfacs == None:
            statfacs = copy.deepcopy(self.BNcrit_res)
        else:
            statfacs = sfacs
        data = copy.deepcopy(self.tdata)
        auto_betta = copy.deepcopy(self.BNcrit_pcmat_resdic['auto_betta'])
        pcmat = copy.deepcopy(self.BNcrit_pcmat_resdic['pcmat'])
        exp_var = copy.deepcopy(self.BNcrit_pcmat_resdic['exp_var'])
        vnames = self.vnames
        if crit != None:
            orig_crit = copy.deepcopy(self.conv_crit)
            self.conv_crit = crit
        orig_sfacs = copy.deepcopy(self.sfacs)
        orig_flags = copy.deepcopy(self.flags)
        critmat = numpy.zeros((statfacs,statfacs))
        joblist = []
        for row in range(0,critmat.shape[0],1):
            for col in range(0,critmat.shape[1],1):
                critmat[row,col] = numpy.nan
        for stato in range(1,statfacs+1,1):
            for dyno in range(1,stato+1,1):
                joblist.append([stato,dyno,job_server.submit(self.refac_pp2,(stato,dyno,data,nlags,pcmat,statfacs),modules=("copy","numpy","scipy"))])
        for job in joblist:
            print 'Now computing combinations for '+str(job[0])+' static factors... (using '+str(job[1])+' dynamic factors)'
            critmat[job[0]-1,job[1]-1] = job[2]()
        self.dynmat = copy.deepcopy(critmat)
        self.dynmat_round = numpy.round(self.dynmat,3)
        minval = numpy.nanmin(self.dynmat)
        minindi = []
        for col in range(0,self.dynmat.shape[1],1):
            for row in range(0,self.dynmat.shape[0],1):
                if self.dynmat[row,col] == minval:
                    minindi.append(row+1)
                    minindi.append(col+1)
        self.dynindex = copy.deepcopy(minindi)
        # One more round with the found r and q to get the results
        self.dynresdic = {}
        stato = self.dynindex[0]
        dyno = self.dynindex[1]
        min_lagos = int(numpy.trunc(stato/dyno))
        rest_lagos = stato-dyno*min_lagos
        if rest_lagos == 0:
            larr = numpy.ones((dyno,min_lagos))
            end_lagos = min_lagos
        elif rest_lagos > 0:
            larr = numpy.ones((dyno,min_lagos+1))
            larr[:,-1] = 0.0
            end_lagos = min_lagos+1
            for elem in range(0,rest_lagos,1):
                larr[elem,-1] = 1.0
        lxmatdic = genLX(self,ydata=data,nlags=nlags,mlags=nlags,func=True)
        self.dynresdic['lxmatdic'] = copy.deepcopy(lxmatdic)
        if dyno == stato+1:
            fdata = self.genLFS(data=pcmat,wcurr=True,larr=larr,func=True)
        else:
            fdata = self.genLFS(data=pcmat[:,:dyno],wcurr=True,larr=larr,func=True)
        self.dynresdic['fdata'] = copy.deepcopy(fdata)
        lxfmatdic = self.genLXF(xdic=lxmatdic,fdata=fdata,func=True)
        self.dynresdic['lxfmatdic'] = copy.deepcopy(lxfmatdic)
        self.dynresdic['exp_var'] = copy.deepcopy(exp_var)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=None,mlags=max(nlags,end_lagos+rest_lagos),func=True)
        # Create the updated diagonal weighting matrix to control for heteroscedasticity ONLY
        diago = numpy.diag(numpy.cov(ee_mat_s))
        hweightm = numpy.zeros((len(diago),len(diago)))
        for diagi in range(0,len(diago),1):
            hweightm[diagi,diagi] = diago[diagi]
        hweightm = numpy.linalg.inv(hweightm)
        auto_betta,auto_ee_dic,ee_mat_s = genEE(self,ydata=data,xdata=lxfmatdic,betta=None,hweightm=hweightm,mlags=max(nlags,end_lagos),func=True)
        self.dynresdic['auto_betta'] = copy.deepcopy(auto_betta)
        self.dynresdic['ee_mat'] = copy.deepcopy(ee_mat_s)
        yfmat,resmat,covmat = self.mkfitres_auto(ymat=data,xdic=lxfmatdic,bettadic=auto_betta,mlags=max(nlags,end_lagos),func=True)
        self.dynresdic['yfmat'] = copy.deepcopy(yfmat)
        self.dynresdic['covmat'] = copy.deepcopy(covmat)
        print 'The model appears to be best represented by '+str(self.dynindex[0])+' static and '+str(self.dynindex[1])+' dynamic factors'
        if self.dynindex[0] == self.dynindex[1]:
            print 'Since in this model q=r, the model is best represented as a static factor model with '+str(self.dynindex[0])+' static factors'

    # Compute the factor innovations and corresponding covmat using the easy method of applying
    # principal components analysis to the x-errors from Stock and Watson 2005
    # Also, compute the corresponding G matrix
    def compFacErrors(self,resmat=None,func=False):
        sfacs = self.sfacs
        if resmat == None:
            resmat = copy.deepcopy(self.dynresdic['ee_mat'])
        pcm = pca_module.PCA_svd(resmat,standardize=False)[0][:,:sfacs]
        pcmm = numpy.zeros((len(pcm),len(pcm[0])))
        for factor in range(0,len(pcm),1):
            pcmm[factor,:] = pcm[factor]
        if not func:
            self.facresmat = copy.deepcopy(pcmm)
            self.facres_covmat = copy.deepcopy(numpy.cov(pcmm.T))
            # Also solve the residual idiosyncratic shocks matrix
            self.idioresmat = copy.deepcopy(numpy.dot(pcmm,self.bettaf_mat.T)-resmat)
        elif func:
            return pcmm

####################################################### COMPANION MATRICES GENERATORS ########################################
##############################################################################################################################
    # Function to build the companion matrix corresponding to the data autoregressors
    def mkCompAuto(self,bettaxm=None,func=False):
        cols = bettaxm.shape[0]
        nlags = self.nlags
        eyes = numpy.eye(cols)
        compmat = numpy.zeros(((cols*nlags),(cols*nlags)))
        for colo in range(0,cols,1):
            for lago in range(0,nlags,1):
                compmat[colo,lago*cols+colo] = bettaxm[colo,lago]
        for lago in range(0,(nlags-1),1):
            compmat[(cols*(lago+1)):(cols*(lago+2)),(cols*lago):(cols*(lago+1))] = eyes
        if not func:
            self.auto_betta_comp = copy.deepcopy(compmat)
        elif func:
            return compmat
        
    # Function to build the companion matrix corresponding to the factors in the VAR for the x data
    def mkCompFactor(self,deltam=None,func=False):
        if deltam == None:
            deltam = copy.deepcopy(self.bettaf_mat)
        cols = deltam.shape[0]
        nlags = self.nlags
        sfacs = self.bettaf_mat.shape[1]
        dynos = self.dynindex[1]
        minflags = int(dynos/sfacs)
        restlags = int(sfacs-dynos*minflags)
        eyes = numpy.eye(sfacs)
        compmat = numpy.zeros(((cols*sfacs),sfacs))
        for faco in range(0,sfacs,1):
            compmat[(cols*faco):(cols*(faco+1)),:] = deltam
        if not func:
            self.auto_fbetta_comp = copy.deepcopy(compmat)
        elif func:
            return compmat

    # Function to build the companion matrix corresponding to the VAR in the factors    
    def mkCompFF(self,bettafm=None,func=False):
        cols = bettafm.shape[0]
        flags = self.flags
        sfacs = bettafm.shape[0]
        eyes = numpy.eye(sfacs)
        compmat = numpy.zeros(((sfacs*flags),(sfacs*flags)))
        compmat[:sfacs,:] = bettafm
        for lago in range(0,(flags-1),1):
            compmat[(sfacs*(lago+1)):(sfacs*(lago+2)),(sfacs*lago):(sfacs*(lago+1))] = eyes
        if not func:
            self.ff_betta_comp = copy.deepcopy(compmat)
        elif func:
            return compmat
    
    # Function to build the FAVAR restricted companion matrix, to get a stacked FAVAR(1)    
    def mkCompAll(self,ff_betta_comp=None,deltaphi=None,dlmat=None,func=False):
        if ff_betta_comp == None: ff_betta_comp = copy.deepcopy(self.ff_betta_comp)
        if deltaphi == None: deltaphi = copy.deepcopy(self.deltaphi)
        if dlmat == None: dlmat = copy.deepcopy(self.dlmat)
        cols = dlmat.shape[0]
        nlags = self.nlags
        dfacs = self.bettaf_mat.shape[1]
        flags = copy.deepcopy(self.flags)
        compmat = numpy.zeros((dfacs*flags+cols*nlags,dfacs*flags+cols*nlags))
        compmat[-dfacs*flags:,:dfacs*flags] = copy.deepcopy(ff_betta_comp)
        compmat[:cols,:dfacs*flags] = copy.deepcopy(deltaphi)
        compmat[:cols,dfacs*flags:] = copy.deepcopy(dlmat)
        eyes = numpy.eye(cols)
        for lago in range(0,nlags-1,1):
            compmat[dfacs*flags+cols*(lago+1):dfacs*flags+cols*(lago+2),cols*lago:cols*(lago+1)] = eyes
        if not func:
            self.varcomp = copy.deepcopy(compmat)
        else:
            return compmat
        
################################################# IRFs AND BOOTSTRAPPING-RELATED METHODS ##################################
###########################################################################################################################
    def mkphis(self,betta=None,maxphi=None,cholmat=None,cholimp=None,compmats=None,rescale=False,func=False):
        if maxphi == None: maxphi = self.confdic['maxphi']
        if cholimp == None : cholimp = self.confdic['cholimp']
        if cholmat == None: cholmat = copy.deepcopy(self.H_final)
        loglist = self.confdic['loglist']
        vnames = self.vnames
        filter_bool = self.confdic['xdata_filter']
        const = self.confdic['const_term']
        transdic = copy.deepcopy(self.trans_dic)
        lags = self.nlags
        flags = self.flags
        # Now try the new version using the three matrices separately and then multiply them
        if compmats == None:
            #auto_betta_comp = copy.deepcopy(self.auto_betta_comp)
            bettaf_mat = copy.deepcopy(self.bettaf_mat)
            ff_betta_comp = copy.deepcopy(self.ff_betta_comp)
            if filter_bool: bettax_mat = copy.deepcopy(self.bettax_mat)
        else:
            #auto_betta_comp = compmats[0]
            bettaf_mat = compmats[1]
            ff_betta_comp = compmats[2]
            if filter_bool: bettax_mat = compmats[3]
        # Number of dynamic factors found by previous method
        df = bettaf_mat.shape[1]
        if filter_bool: k = bettax_mat.shape[0]
        else: k = bettaf_mat.shape[0]
        # Matrices for the factor autoregressive coefficients of order flag
        phis_compta2 = numpy.zeros((maxphi+1,df,df))
        phis_compta2[0] = numpy.eye(df)
        # Matrices for the x-data autoregressive coefficients of order nlag
        phis_compta = numpy.zeros((maxphi+1,k,df))
        phis_compta[0] = numpy.ones((phis_compta.shape[1],phis_compta.shape[2]))
        # Matrices for total multiplicative effect of all of the series, this has dim ncols*sfacs
        phis_compta3 = numpy.zeros((maxphi+1,k,df))
        phis_compta4 = numpy.zeros((lags,k,df))
        if filter_bool:
            for j in xrange(0,lags,1):
                phis_compta4[j] = numpy.tile(numpy.reshape(bettax_mat[:,j],(k,1)),(1,df))
        # recursively compute Phi matrices
        if filter_bool:
            # phis_compta contains the autoregressive matrices for the auto-regressive elements of the x data only
            for i in xrange(1, maxphi+1):
                for j in xrange(1,i+1):
                    if j > lags: break
                    phis_compta[i] += phis_compta[i-j]*numpy.tile(numpy.reshape(bettax_mat[:,j-1],(k,1)),(1,df))
        # phis_compta2 contains the autoregressive matrices for the auto-regressive elements of the factors only
        for i in xrange(1, maxphi+1):
            for j in xrange(1,i+1):
                if j > flags: break
                phis_compta2[i] += numpy.dot(phis_compta2[i-j],ff_betta_comp[:df,(j-1)*df:j*df])
        # phis_compta2c contains the autoregressive matrices for the auto-regressive elements of the factors only, with cholmat
        phis_compta2c = numpy.array([list(numpy.dot(tmat, cholmat)) for tmat in phis_compta2])
        if filter_bool:
            # In phis_compta3, everything is multiplied up to get the full dynamic response
            for kk in xrange(0,maxphi+1):
                for j in xrange(0,kk):
                    if j > lags-1: break
                    phis_compta3[kk] += phis_compta4[j]*numpy.dot(bettaf_mat,phis_compta2[kk-j]) 
        else:
            # In phis_compta3, everything is multiplied up to get the full dynamic response
            for kk in xrange(0,maxphi+1):
                phis_compta3[kk] = numpy.dot(bettaf_mat,phis_compta2[kk])            
        if filter_bool:           
            # Also create the orthogonalized equivalents
            phisc = numpy.zeros((phis_compta3.shape[0],phis_compta3.shape[1],phis_compta3.shape[2]))
            for kk in xrange(0,maxphi+1):
                for j in xrange(0,kk):
                    if j > lags-1: break
                    phisc[kk] += phis_compta4[j]*numpy.dot(bettaf_mat,phis_compta2c[kk-j])
        else:
            # Also create the orthogonalized equivalents
            phisc = numpy.zeros((phis_compta3.shape[0],phis_compta3.shape[1],phis_compta3.shape[2]))
            for kk in xrange(0,maxphi+1):
                phisc[kk] = numpy.dot(bettaf_mat,phis_compta2c[kk])
        phisc_unscaled = copy.deepcopy(phisc)
        if rescale:
            for i1,namo in enumerate(vnames):
                    phisc[:,i1] = phisc[:,i1]*transdic[namo]['std']
        if not func:
            self.phis_compta2c = phis_compta2c # This is for the factors only, but with cholmat
            self.phis_compta2 = phis_compta2 # This is for the factors only, without cholmat
            self.phis_compta = phis_compta
            self.phis_compta3 = phis_compta3 # This is the final response without the cholmat
            self.phisc = phisc # This is the final response with the cholmat, but scaled
            self.phisc_unscaled = copy.deepcopy(phisc_unscaled) # This is the final response with cholmat and scaled
        elif func:
            return phis_compta2,phisc,phisc_unscaled,phis_compta2c

    def mkboot(self,bdraw=None,yfmat=None,resiarray=None,rescale=False):
        # Set instance into boot mode
        self.boot = True
        # Get all of the paramters
        if bdraw == None: bdraw = self.confdic['bdraw']
        breplace = self.confdic['breplace']
        maxphi = self.confdic['maxphi']
        signif = self.confdic['signif']
        irfp = self.confdic['irfp']
        const_term = self.confdic['const_term']
        resmat = self.dynresdic['ee_mat']
        if self.confdic['use_bbe_ident']: resmat_facs = self.bbe_resdic['ee_mat']
        else: resmat_facs = copy.deepcopy(self.resmat_ff)
        vnames = self.vnames
        cols = self.ncols
        dfacs = self.dynindex[1]
        flags = self.flags
        if yfmat == None: yfmat = self.dynresdic['yfmat']
        nlags = self.nlags
        matd = self.tdata
        # Produce bootstrap confidence intervals
        barray = []
        betarray = []
        rarray = []
        resarray = []
        nmatarray = []
        betarray_n = []
        for bodraw in range(0,bdraw,1):
            if resiarray == None:
                resmat_n = copy.deepcopy(resmat)
                if breplace:
                    matd_n = numpy.zeros((yfmat.shape[0],yfmat.shape[1]))
                    for elem in range(0,yfmat.shape[0],1):
                        numpy.random.shuffle(resmat_n)
                        matd_n[elem,:] = yfmat[elem,:] + resmat_n[0,:]
                    restmp = matd_n - yfmat
                    resarray.append(restmp)
                else:
                    numpy.random.shuffle(resmat_n)
                    resarray.append(resmat_n)
                    matd_n = yfmat + resmat_n
            elif resiarray != None:
                matd_n = yfmat + resiarray[bodraw]
                resarray.append(resiarray[bodraw])
            # Add actual obs from the original real data to make up for lost observations from estimation
            matd_n = numpy.vstack((matd[:nlags,:],matd_n))
            self.boot_matd = copy.deepcopy(matd_n)
            rarray.append(copy.deepcopy(matd_n))
            rows_n,cols_n = matd_n.shape
            nmatd_n = copy.deepcopy(self.dynresdic['lxfmatdic'])
            nmatarray.append(nmatd_n)
            betta_n = genEE(self,ydata=matd_n,xdata=nmatd_n,func=True)[0]
            # Produce matrix from auto_betta
            auto_betta_mat = betta_n[vnames[0]]
            for varo in range(1,cols,1):
                auto_betta_mat = numpy.vstack((auto_betta_mat,betta_n[vnames[varo]]))
            betarray_n.append(auto_betta_mat)
            bettax_mat = copy.deepcopy(auto_betta_mat[:,dfacs:])
            bettaf_mat = copy.deepcopy(auto_betta_mat[:,:dfacs])
            # This creates relevant companion matrices which will be easier to work with later
            auto_betta_comp = self.mkCompAuto(bettaxm=bettax_mat,func=True)
            ff_betta_comp = self.mkCompFF(bettafm=self.fbetta,func=True)
            auto_fbetta_comp = self.mkCompFactor(deltam = bettaf_mat,func=True)
            # Create the VAR matrices for the DFM in VAR form
            deltaphi = numpy.dot(bettaf_mat,ff_betta_comp)
            dlmat = auto_betta_comp[:self.ncols,:]
            # Stack all matrices together to one large one
            varcomp_n = self.mkCompAll(ff_betta_comp=ff_betta_comp,deltaphi=deltaphi,dlmat=dlmat,func=True)            
            betarray.append(varcomp_n)
            yfmat_n,resmat_n,vcmat_n = self.mkfitres_auto(ymat=matd_n,xdic=nmatd_n,bettadic=betta_n,func=True)
            compmats = [auto_betta_comp,bettaf_mat,ff_betta_comp,bettax_mat]
            phis_n,phisc_n = self.mkphis(betta=varcomp_n,maxphi=maxphi,compmats=compmats,rescale=rescale,func=True)
            barray.append(phisc_n)
        nmatarray = numpy.array(nmatarray)
        resarray = numpy.array(resarray)
        barray = numpy.array(barray)
        betarray = numpy.array(betarray)
        betarray_n = numpy.array(betarray_n)
        rarray = numpy.array(rarray)
        if type(signif) == type(10.0):
            lower = scipy.stats.scoreatpercentile(barray,per=(signif/2.0))
            upper = scipy.stats.scoreatpercentile(barray,per=(100.0-(signif/2.0)))
            confint = [lower,upper]
        elif signif == 'mix':
            lower_outer = scipy.stats.scoreatpercentile(barray,per=5.0)
            lower_inner = scipy.stats.scoreatpercentile(barray,per=(34.0/2.0))
            upper_inner = scipy.stats.scoreatpercentile(barray,per=100.0-(34.0/2.0))
            upper_outer = scipy.stats.scoreatpercentile(barray,per=95.0)
            confint = [lower_outer,lower_inner,upper_inner,upper_outer]
        self.rarray = rarray
        self.betarray = betarray
        self.betarray_n = betarray_n
        self.barray = barray
        self.confint = numpy.array(confint)[:,:irfp+1]
        self.resarray = resarray
        self.nmatarray = nmatarray
        # Switch out of boot mode
        self.boot = False
        

    def refac_killian_pp(self,intup):
        resmat = intup[0]
        yfmat = intup[1]
        matd = intup[2]
        nlags = intup[3]
        maxphi = intup[4]
        breplace = intup[5]
        rescale = intup[6]
        dfacs = intup[7]
        fbetta = intup[8]
        if self.confdic['use_bbe_ident']:
            vnames = copy.deepcopy(self.vnames_adj)
        else:
            vnames = copy.deepcopy(self.vnames)
        # Get all of the matrices needed to compute new actual values
        if self.confdic['xdata_filter']: bettax_mat = copy.deepcopy(self.bettax_mat)
        bettaf_mat = copy.deepcopy(self.bettaf_mat)
        ffbetta = copy.deepcopy(fbetta)
        flags = int(ffbetta.shape[1]/dfacs)
        # Also get the fitted factor series from the one-order lag VAR in the factors
        if self.confdic['use_bbe_ident']: ffmat = self.bbe_resdic['ffmat']
        else: ffmat = self.fdata[flags:,:] - self.resmat_ff
        # This is the composite matrix in the x-data equation(s) linking them to the lagged factors
        phimat = numpy.dot(bettaf_mat,ffbetta)
        if self.confdic['use_bbe_ident']: resmat_facs_n = copy.deepcopy(self.bbe_resdic['ee_mat'])
        else: resmat_facs_n = copy.deepcopy(self.resmat_ff)
        if self.confdic['use_bbe_ident']: fdata = copy.deepcopy(self.fdata_adj)
        else: fdata = copy.deepcopy(self.fdata)
        # Also create data series of fdata with appropriate number of lags
        fdata_lagged = self.genLFF(fdata=fdata,flags=flags,func=True)
        if self.confdic['xdata_filter']: xdic = genLX(self,vnames=vnames,ydata=matd,nlags=nlags,mlags=None,func=True)
        resmat_n = copy.deepcopy(resmat)
        cols = yfmat.shape[1]
        fcols = fdata.shape[1]
        frows = fdata.shape[0]
        # Reshuffle the data innovations
        # numpy.random.shuffle(resmat_n)
        # Reshuffle the factor residuals
        numpy.random.shuffle(resmat_facs_n)
        # This part updates the VAR(p) coefficient matrix in the factors (due to sampling uncertainty)
        init_vals = fdata[:flags,:]
        fdata_n = numpy.zeros((frows,fcols))
        fdata_n[:flags,:] = copy.deepcopy(init_vals)
        for obso in xrange(0,(frows-flags),1):
            for j in xrange(1,flags+1):
                fdata_n[flags+obso,:] += numpy.dot(fdata_n[flags+obso-j,:],ffbetta[:,(j-1)*fcols:j*fcols].T)
            fdata_n[flags+obso,:] += resmat_facs_n[obso,:]
        ffbetta_n = pymaclab.stats.common.compbetta(self,matd=fdata_n,nmatd=fdata_lagged,smbetta=False,const=None,func=True)
        # This creates the updated composite matrix in the x-data equations conditional on the lagged factors
        # should be (k*df)*(df,df*flags)=(k,df*flags)
        phimat_n = numpy.dot(bettaf_mat,ffbetta_n)
        matd_n = yfmat + resmat_n
        if self.confdic['xdata_filter']:
            auto_fitted = numpy.zeros((matd_n.shape[0],matd_n.shape[1]))
            for i1,namo in enumerate(vnames):
                auto_fitted[:,i1] = numpy.dot(xdic[namo],bettax_mat[i1,:].T)
            # Need to drop first 'flags' observation on the two arrays as we employ the lagged coefficients of the factors
            auto_fitted = auto_fitted[flags:,:]
            resmat_n = resmat_n[flags:,:]
            # This creates a new matrix of actual (resamples) values of the x data
            # Should be (n*k) = (n*k) + (n,df*flags)*(k,df*flags).T
            matd_n = auto_fitted + numpy.dot(fdata_lagged,phimat_n.T) + resmat_n
            # Add actual obs from the original real data to make up for lost observations from estimation
            matd_n = numpy.vstack((matd[:nlags,:],matd_n))
            self.boot_matd = copy.deepcopy(matd_n)
            rows_n,cols_n = matd_n.shape
            nmatd_n = copy.deepcopy(self.dynresdic['lxfmatdic'])
            betta_n = genEE(self,ydata=matd_n,xdata=nmatd_n,func=True)[0]
            # Produce matrix from auto_betta
            auto_betta_mat = betta_n[vnames[0]]
            for varo in range(1,cols,1):
                auto_betta_mat = numpy.vstack((auto_betta_mat,betta_n[vnames[varo]]))
                bettax_mat = copy.deepcopy(auto_betta_mat[:,dfacs:])
                bettaf_mat = copy.deepcopy(auto_betta_mat[:,:dfacs])
            # This creates relevant companion matrices which will be easier to work with later
            auto_betta_comp = self.mkCompAuto(bettaxm=bettax_mat,func=True)
            ff_betta_comp = self.mkCompFF(bettafm=ffbetta_n,func=True)
            auto_fbetta_comp = self.mkCompFactor(deltam = bettaf_mat,func=True)
            # Create the VAR matrices for the DFM in VAR form
            deltaphi = numpy.dot(bettaf_mat,ff_betta_comp[:dfacs,:])
            dlmat = auto_betta_comp[:cols,:]
            # Stack all matrices together to one large one
            varcomp_n = self.mkCompAll(ff_betta_comp=ff_betta_comp,deltaphi=deltaphi,dlmat=dlmat,func=True)
            yfmat_n,resmat_n,vcmat_n = self.mkfitres_auto(ymat=matd_n,xdic=nmatd_n,bettadic=betta_n,func=True)
            compmats = [auto_betta_comp,bettaf_mat,ff_betta_comp,bettax_mat]
        else:
            nmatdf_n = self.genLFF(fdata=fdata_n,func=True)
            bettaf_n = pymaclab.stats.common.compbetta(self,matd=fdata_n,nmatd=nmatdf_n,smbetta=False,const=None,func=True)
            ff_betta_comp = self.mkCompFF(bettafm=bettaf_n,func=True)
            ffmat,fresmat,fcovmat = self.mkfitres_ffvar(ymat=fdata_n,xmat=nmatdf_n,betta=bettaf_n,flags=self.flags,func=True)
            idioresmat_n = copy.deepcopy(fresmat)
            compmats = [None,bettaf_mat,ff_betta_comp,None]            
        return bettaf_n,ff_betta_comp


    def refac_pp(self,intup):
        resmat = intup[0]
        yfmat = intup[1]
        matd = intup[2]
        nlags = intup[3]
        maxphi = intup[4]
        breplace = intup[5]
        rescale = intup[6]
        dfacs = intup[7]
        fbetta = intup[8]
        if self.confdic['use_bbe_ident']:
            vnames = copy.deepcopy(self.vnames_adj)
        else:
            vnames = copy.deepcopy(self.vnames)
        # Get all of the matrices needed to compute new actual values
        if self.confdic['xdata_filter']: bettax_mat = copy.deepcopy(self.bettax_mat)
        bettaf_mat = copy.deepcopy(self.bettaf_mat)
        ffbetta = copy.deepcopy(fbetta)
        flags = int(ffbetta.shape[1]/dfacs)
        # Also get the fitted factor series from the one-order lag VAR in the factors
        if self.confdic['use_bbe_ident']: ffmat = self.bbe_resdic['ffmat']
        else: ffmat = self.fdata[flags:,:] - self.resmat_ff
        # This is the composite matrix in the x-data equation(s) linking them to the lagged factors
        phimat = numpy.dot(bettaf_mat,ffbetta)
        if self.confdic['use_bbe_ident']: resmat_facs_n = copy.deepcopy(self.bbe_resdic['ee_mat'])
        else: resmat_facs_n = copy.deepcopy(self.resmat_ff)
        if self.confdic['use_bbe_ident']: fdata = copy.deepcopy(self.fdata_adj)
        else: fdata = copy.deepcopy(self.fdata)
        # Also create data series of fdata with appropriate number of lags
        fdata_lagged = self.genLFF(fdata=fdata,flags=flags,func=True)
        if self.confdic['xdata_filter']: xdic = genLX(self,vnames=vnames,ydata=matd,nlags=nlags,mlags=None,func=True)
        resmat_n = copy.deepcopy(resmat)
        cols = yfmat.shape[1]
        fcols = fdata.shape[1]
        frows = fdata.shape[0]
        # Reshuffle the data innovations
        # numpy.random.shuffle(resmat_n)
        # Reshuffle the factor residuals
        numpy.random.shuffle(resmat_facs_n)
        # This part updates the VAR(p) coefficient matrix in the factors (due to sampling uncertainty)
        init_vals = fdata[:flags,:]
        fdata_n = numpy.zeros((frows,fcols))
        fdata_n[:flags,:] = copy.deepcopy(init_vals)
        for obso in xrange(0,(frows-flags),1):
            for j in xrange(1,flags+1):
                fdata_n[flags+obso,:] += numpy.dot(fdata_n[flags+obso-j,:],ffbetta[:,(j-1)*fcols:j*fcols].T)
            fdata_n[flags+obso,:] += resmat_facs_n[obso,:]
        ffbetta_n = pymaclab.stats.common.compbetta(self,matd=fdata_n,nmatd=fdata_lagged,smbetta=False,const=None,func=True)
        # This creates the updated composite matrix in the x-data equations conditional on the lagged factors
        # should be (k*df)*(df,df*flags)=(k,df*flags)
        phimat_n = numpy.dot(bettaf_mat,ffbetta_n)
        matd_n = yfmat + resmat_n
        if self.confdic['xdata_filter']:
            auto_fitted = numpy.zeros((matd_n.shape[0],matd_n.shape[1]))
            for i1,namo in enumerate(vnames):
                auto_fitted[:,i1] = numpy.dot(xdic[namo],bettax_mat[i1,:].T)
            # Need to drop first 'flags' observation on the two arrays as we employ the lagged coefficients of the factors
            auto_fitted = auto_fitted[flags:,:]
            resmat_n = resmat_n[flags:,:]
            # This creates a new matrix of actual (resamples) values of the x data
            # Should be (n*k) = (n*k) + (n,df*flags)*(k,df*flags).T
            matd_n = auto_fitted + numpy.dot(fdata_lagged,phimat_n.T) + resmat_n
            # Add actual obs from the original real data to make up for lost observations from estimation
            matd_n = numpy.vstack((matd[:nlags,:],matd_n))
            self.boot_matd = copy.deepcopy(matd_n)
            rows_n,cols_n = matd_n.shape
            nmatd_n = copy.deepcopy(self.dynresdic['lxfmatdic'])
            betta_n = genEE(self,ydata=matd_n,xdata=nmatd_n,func=True)[0]
            # Produce matrix from auto_betta
            auto_betta_mat = betta_n[vnames[0]]
            for varo in range(1,cols,1):
                auto_betta_mat = numpy.vstack((auto_betta_mat,betta_n[vnames[varo]]))
                bettax_mat = copy.deepcopy(auto_betta_mat[:,dfacs:])
                bettaf_mat = copy.deepcopy(auto_betta_mat[:,:dfacs])
            # This creates relevant companion matrices which will be easier to work with later
            auto_betta_comp = self.mkCompAuto(bettaxm=bettax_mat,func=True)
            ff_betta_comp = self.mkCompFF(bettafm=ffbetta_n,func=True)
            auto_fbetta_comp = self.mkCompFactor(deltam = bettaf_mat,func=True)
            # Create the VAR matrices for the DFM in VAR form
            deltaphi = numpy.dot(bettaf_mat,ff_betta_comp[:dfacs,:])
            dlmat = auto_betta_comp[:cols,:]
            # Stack all matrices together to one large one
            varcomp_n = self.mkCompAll(ff_betta_comp=ff_betta_comp,deltaphi=deltaphi,dlmat=dlmat,func=True)
            yfmat_n,resmat_n,vcmat_n = self.mkfitres_auto(ymat=matd_n,xdic=nmatd_n,bettadic=betta_n,func=True)
            compmats = [auto_betta_comp,bettaf_mat,ff_betta_comp,bettax_mat]
        else:
            nmatdf_n = self.genLFF(fdata=fdata_n,func=True)
            bettaf_n = pymaclab.stats.common.compbetta(self,matd=fdata_n,nmatd=nmatdf_n,smbetta=False,const=None,func=True)
            ff_betta_comp = self.mkCompFF(bettafm=bettaf_n,func=True)
            ffmat,fresmat,fcovmat = self.mkfitres_ffvar(ymat=fdata_n,xmat=nmatdf_n,betta=bettaf_n,flags=self.flags,func=True)
            idioresmat_n = copy.deepcopy(fresmat)
            compmats = [None,bettaf_mat,ff_betta_comp,None]
        phis_n,phisc_n,phisc_unscaled_n,phiscf_n = self.mkphis(betta=None,maxphi=maxphi,compmats=compmats,rescale=rescale,func=True)
        return phisc_n,phiscf_n


    def mkboot_pp(self,bdraw=None,yfmat=None,resiarray=None,rescale=False):
        # Set instance into boot mode
        self.boot = True
        # Get all of the paramters
        if bdraw == None: bdraw = self.confdic['bdraw']
        breplace = self.confdic['breplace']
        maxphi = self.confdic['maxphi']
        signif = self.confdic['signif']
        irfp = self.confdic['irfp']
        const_term = self.confdic['const_term']
        resmat = self.dynresdic['ee_mat']
        fbetta = copy.deepcopy(self.fbetta)
        if self.confdic['use_bbe_ident']: resmat_facs = self.bbe_resdic['ee_mat']
        else: resmat_facs = copy.deepcopy(self.resmat_ff)
        bdraw_chunk = self.confdic['bdraw_chunk']
        vnames = self.vnames
        dfacs = resmat_facs.shape[1]
        flags = self.flags
        if yfmat == None: yfmat = self.dynresdic['yfmat']
        cols = yfmat.shape[1]
        nlags = self.nlags
        matd = self.tdata
        # Produce bootstrap confidence intervals, here using parallel python to speed up computation
        barray = []
        ppservers = ()
        if self.confdic['ncpus'] == 'auto':
            job_server = pp.Server(ncpus='autodetect',ppservers=ppservers)
        else:
            job_server = pp.Server(ncpus=self.confdic['ncpus'],ppservers=ppservers)
        split_bdraw = int(bdraw/bdraw_chunk)
        split_mod = bdraw%bdraw_chunk
        # Do the first bootstrap run a la Killian 1998 to reduce the small sample bias of the IRFs, if wanted
        if self.confdic['killian_bsinbs']:
            print "Calculating bias correction using Killian's(1998) bs-in-bs method"
            barray_biasred = []
            # Note: Only pass the betta coefficients of the dynamic responses, not constant or time-trend
            inputs_biasred = ((resmat,yfmat,matd,nlags,maxphi,fbetta),)*bdraw_chunk
            for elem in xrange(0,split_bdraw):
                inputs = ((resmat,yfmat,matd,nlags,maxphi,breplace,rescale,dfacs,fbetta),)*bdraw_chunk
                barray_jobs = [job_server.submit(self.refac_killian_pp,(input,),modules=("copy","numpy","scipy","pymaclab.stats.common")) for input in inputs]
                barray_biasred = barray_biasred + [x() for x in barray_jobs]
            if split_mod != 0:
                inputs = ((resmat,yfmat,matd,nlags,maxphi,breplace,rescale,dfacs,fbetta),)*split_mod
                barray_jobs = [job_server.submit(self.refac_killian_pp,(input,),modules=("copy","numpy","scipy","pymaclab.stats.common")) for input in inputs]
                barray_biasred = barray_biasred + [x() for x in barray_jobs]
            barray_biasred = numpy.array([x[0] for x in barray_biasred])
            self.barray_biasred = copy.deepcopy(barray_biasred)
            biasred_term = numpy.mean(barray_biasred-fbetta,axis=0)
            # Make the bias in value and percentage available for inspection
            self.killbias_ff = copy.deepcopy(biasred_term)
            self.killbiasperc_ff = copy.deepcopy((biasred_term/fbetta)*100.0)
            self.betta_kcorr_ff = copy.deepcopy(fbetta + biasred_term)
            betta_corrected = fbetta + biasred_term
            # Lets keep the original residuals, comment the following out if you want the new ones
            '''
            nmatd = copy.deepcopy(self.nmatd)
            yfmat = numpy.dot(nmatd,betta_corrected.T)
            resmat = matd[nlags:,:] - yfmat
            # Need to centre residuals prior to using them in the bootstrap procedure
            resmat = resmat - numpy.mean(resmat,axis=0)
            '''
            # Pass the corrected betta on to the boostrap procedure
            fbetta = copy.deepcopy(betta_corrected)
            # Recalculation with new betta
            #nmatdf = self.genLFF(fdata=self.fdata_adj,flags=self.flags,func=True)
            #yfmat,resmat,covmat = self.mkfitres_ffvar(ymat=self.fdata_adj,
            #                                        xmat=nmatdf,
            #                                        betta=betta_corrected,
            #                                        flags=self.flags,
            #                                        func=True)
            #self.covmat_ff_kcorr = copy.deepcopy(covmat)
            #self.resmat_ff_kcorr = copy.deepcopy(resmat)
            #self.yfmat_ff_kcorr = copy.deepcopy(yfmat)            
            #yfmat_auto = numpy.dot(nmatdf,numpy.dot(self.bettaf_mat,betta_corrected).T)
            #resmat_auto = self.tdata[self.flags:,:] - yfmat_auto
            #covmat_auto = numpy.cov(resmat_auto.T)
            #self.covmat_auto_kcorr = copy.deepcopy(covmat_auto)
            #self.H_final = numpy.linalg.cholesky(self.covmat_ff_kcorr)
            # Needs to be called before new impulses are created, because colmat needs to be updated
            self.betta_uncorr_ff = copy.deepcopy(self.fbetta)
            self.phisc_uncorr_ff = copy.deepcopy(self.phisc)
            ff_betta_comp = self.mkCompFF(bettafm=betta_corrected,func=True)
            compmats = [None,self.bettaf_mat,ff_betta_comp,None]
            phisf_kcorr,phisc_kcorr,phisc_kcorr_unscaled,phiscf_kcorr = self.mkphis(betta=self.betta_kcorr_ff,compmats=compmats,cholmat=self.H_final,rescale=rescale,func=True)
            self.phisc_kcorr = copy.deepcopy(phisc_kcorr)
            self.phisc_uncorr = copy.deepcopy(self.phisc)
            self.phisc_kcorr_unscaled = copy.deepcopy(phisc_kcorr_unscaled)
            self.phiscf_kcorr = copy.deepcopy(phiscf_kcorr)
            self.phisf_kcorr = copy.deepcopy(phisf_kcorr)
            self.phisc = copy.deepcopy(self.phisc_kcorr)
            print 'Calculation of bias correction done...'
        barray = []
        barrayf = []
        for elem in xrange(0,split_bdraw):
            inputs = ((resmat,yfmat,matd,nlags,maxphi,breplace,rescale,dfacs,fbetta),)*bdraw_chunk
            barray_jobs = [job_server.submit(self.refac_pp,(input,),modules=("copy","numpy","scipy","pymaclab.stats.common")) for input in inputs]
            barray = barray + [x()[0] for x in barray_jobs]
            barrayf = barrayf + [x()[1] for x in barray_jobs]
        if split_mod != 0:
            inputs = ((resmat,yfmat,matd,nlags,maxphi,breplace,rescale,dfacs,fbetta),)*split_mod
            barray_jobs = [job_server.submit(self.refac_pp,(input,),modules=("copy","numpy","scipy")) for input in inputs]
            barray = barray + [x()[0] for x in barray_jobs]
            barrayf = barrayf + [x()[1] for x in barray_jobs]
        barray = numpy.array(barray)
        barrayf = numpy.array(barrayf)
        if type(signif) == type(10.0):
            lower = scipy.stats.scoreatpercentile(barray,per=(signif/2.0))
            upper = scipy.stats.scoreatpercentile(barray,per=(100.0-(signif/2.0)))
            lowerf = scipy.stats.scoreatpercentile(barrayf,per=(signif/2.0))
            upperf = scipy.stats.scoreatpercentile(barrayf,per=(100.0-(signif/2.0)))
            confint = [lower,upper]
            confintf = [lowerf,upperf]
        elif signif == 'mix':
            lower_outer = scipy.stats.scoreatpercentile(barray,per=5.0)
            lower_outerf = scipy.stats.scoreatpercentile(barrayf,per=5.0)
            lower_inner = scipy.stats.scoreatpercentile(barray,per=(34.0/2.0))
            lower_innerf = scipy.stats.scoreatpercentile(barrayf,per=(34.0/2.0))
            upper_inner = scipy.stats.scoreatpercentile(barray,per=100.0-(34.0/2.0))
            upper_innerf = scipy.stats.scoreatpercentile(barrayf,per=100.0-(34.0/2.0))
            upper_outer = scipy.stats.scoreatpercentile(barray,per=95.0)
            upper_outerf = scipy.stats.scoreatpercentile(barrayf,per=95.0)
            confint = [lower_outer,lower_inner,upper_inner,upper_outer]
            confintf = [lower_outerf,lower_innerf,upper_innerf,upper_outerf]
        self.barray = copy.deepcopy(barray)
        self.barrayf = copy.deepcopy(barrayf)
        self.confint = numpy.array(confint)[:,:irfp+1]
        self.confintf = numpy.array(confintf)[:,:irfp+1]
        # Switch out of boot mode
        self.boot = False
    
    def genFEVD(self,y,coefs,steps,func=False):
        p = len(coefs)
        k = len(coefs[0])
        bettaf = copy.deepcopy(self.bettaf_mat)
        vnames = self.vnames
        
        # initial value
        forcs = numpy.zeros((steps, k))
        # h=0 forecast should be latest observation
        # forcs[0] = y[-1]
    
        # make indices easier to think about
        for h in xrange(1, steps + 1):
            # y_t(h) = intercept + sum_1^p A_i y_t_(h-i)
            f = forcs[h - 1]
            for i in xrange(1, p + 1):
                # slightly hackish
                if h - i <= 0:
                    # e.g. when h=1, h-1 = 0, which is y[-1]
                    prior_y = y[h - i - 1]
                else:
                    # e.g. when h=2, h-1=1, which is forcs[0]
                    prior_y = forcs[h - i - 1]
    
                # i=1 is coefs[0]
                f = f + numpy.dot(coefs[i - 1], prior_y)
            forcs[h - 1] = f
        self.forcs = forcs

        ma_coefs = copy.deepcopy([numpy.dot(bettaf,tmat) for tmat in self.phisf_kcorr])
        # ma_coefs = copy.deepcopy(self.phis)
        sigma_uf = self.covmat_ff
        #sigma_u = self.covmat_auto
        k = len(sigma_uf)
        k2 = len(vnames)
        forc_covs = numpy.zeros((steps, k2, k2))
        prior = numpy.zeros((k2, k2))
        for h in xrange(steps):
            # Sigma(h) = Sigma(h-1) + Phi Sig_u Phi'
            phi = ma_coefs[h]
            #inner_one = numpy.dot(bettaf,phi)
            #sigo1 = numpy.dot(resmat,bettaf)
            #sigo = numpy.dot(sigo1.T,sigo1)
            #lefto = numpy.dot(bettaf,phi)
            var = numpy.dot(numpy.dot(phi,sigma_uf),phi.T)
            forc_covs[h] = prior = prior + var
        self.forc_covs = forc_covs
            
        periods = steps
        neqs = self.fdata_adj.shape[1]
        
        if 'phisc_kcorr_unscaled' in dir(self):
            orthirfs = self.phiscf_kcorr
        elif 'phisc_unscaled' in dir(self):
            orthirfs = self.phis_compta2c
        else:
            orthirfs = self.phisc_kcorr
        orthirfs = self.phisc_kcorr_unscaled

        # cumulative impulse responses
        irfs = (orthirfs[:periods] ** 2).cumsum(axis=0)
        self.forc_covs_irfs = copy.deepcopy(irfs)

        rng = range(neqs)
        rng2 = range(k2)
        mse = forc_covs[:, rng2, rng2]
        self.mse = copy.deepcopy(mse)

        # lag x equation x component
        fevd = numpy.empty_like(irfs)

        for i in range(periods):
            fevd[i] = (irfs[i].T / mse[i]).T

        # Express in percentage terms
        fevd = fevd*100.0

        # switch to equation x lag x component
        self.decomp = fevd.swapaxes(0, 1)

################################################# GRAPHING/PLOTTING METHODS ###############################################
###########################################################################################################################
    def plotFactors(self,fdata=None):
        years    = mdates.YearLocator()   # every year
        months   = mdates.MonthLocator()  # every month
        yearsFmt = mdates.DateFormatter('%Y')        
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'factors' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/factors')
        dato = copy.deepcopy(self.dates)
        plotrec = self.confdic['plot_rec']
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']        
        rdates = self.rdates
        if 'rec_colour' in self.confdic.keys(): reccol = self.confdic['rec_colour']
        if 'rec_transp' in self.confdic.keys(): rectrans = self.confdic['rec_transp']
        flags = self.flags
        if fdata == None:
            if 'fdata_adj' in dir(self):
                fdata = copy.deepcopy(self.fdata_adj)
            else:
                fdata = copy.deepcopy(self.fdata)
        for fac_num,faco in enumerate(fdata.T):
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(1,1,1)
            if glopt['xlabel']: ax1.set_xlabel('Time',fontsize=glopt['label_fs'])
            for label in ax1.get_xticklabels():
                label.set_rotation(55)
            if gridopt: ax1.grid()
            ax1.plot(dato,faco,color='black')
            if glopt['title']: ax1.set_title('Estimated Values for FACTOR '+str(fac_num+1),fontsize=glopt['title_fs'])
            fig1.savefig('../graphs/'+modname+'/factors/'+'FACTOR_'+str(fac_num+1)+'.eps',bbox_inches='tight')
        # Also plot unadjusted factors before BBE identification
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'factors' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/factors')
        if 'unadj' not in os.listdir('../graphs/'+modname+'/factors/'):
            os.mkdir('../graphs/'+modname+'/factors/unadj') 
        for fac_num,faco in enumerate(self.fdata.T):
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(1,1,1)
            if glopt['xlabel']: ax1.set_xlabel('Time')
            for label in ax1.get_xticklabels():
                label.set_rotation(55)
            if gridopt: ax1.grid()
            ax1.plot(dato,faco,color='black')
            if plotrec:
                for reco in rdates:
                    plt.fill([reco[0],reco[1],reco[1],reco[0]],
                             [ax1.get_ylim()[0],ax1.get_ylim()[0],ax1.get_ylim()[1],ax1.get_ylim()[1]],
                             reccol, alpha=rectrans, edgecolor=reccol)
            if glopt['title']: ax1.set_title('Estimated Values for FACTOR '+str(fac_num+1),fontsize=glopt['title_fs'])
            fig1.savefig('../graphs/'+modname+'/factors/unadj/'+'FACTOR_'+str(fac_num+1)+'.eps',bbox_inches='tight')
        # Also plot some factor sums
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1,1,1)
        if glopt['xlabel']: ax1.set_xlabel('Time',fontsize=glopt['label_fs'])
        for label in ax1.get_xticklabels():
            label.set_rotation(55)
        if gridopt: ax1.grid()
        ax1.plot(dato,numpy.sum(self.fdata[:,0:4],axis=1).T,color='black')
        if plotrec:
            for reco in rdates:
                plt.fill([reco[0],reco[1],reco[1],reco[0]],
                         [ax1.get_ylim()[0],ax1.get_ylim()[0],ax1.get_ylim()[1],ax1.get_ylim()[1]],
                         reccol, alpha=rectrans, edgecolor=reccol)
        if glopt['title']: ax1.set_title('Sum of Factors 1-4',fontsize=glopt['title_fs'])
        fig1.savefig('../graphs/'+modname+'/factors/unadj/'+'FACTOR_1to4.eps',bbox_inches='tight')  
        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1,1,1)
        if glopt['xlabel']: ax1.set_xlabel('Time',fontsize=glopt['label_fs'])
        for label in ax1.get_xticklabels():
            label.set_rotation(55)
        if gridopt: ax1.grid()
        ax1.plot(dato,numpy.sum(self.fdata[:,4:],axis=1).T,color='black')
        if plotrec:
            for reco in rdates:
                plt.fill([reco[0],reco[1],reco[1],reco[0]],
                         [ax1.get_ylim()[0],ax1.get_ylim()[0],ax1.get_ylim()[1],ax1.get_ylim()[1]],
                         reccol, alpha=rectrans, edgecolor=reccol)
        if glopt['title']: ax1.set_title('Sum of Factors 5-12',fontsize=glopt['title_fs'])
        fig1.savefig('../graphs/'+modname+'/factors/unadj/'+'FACTOR_5to12.eps',bbox_inches='tight')
            
    def plotFactorResids(self,fresids=None):
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'factors' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/factors')
        dato = copy.deepcopy(self.dates)
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']
        if fresids == None:
            fresids = self.resmat_ff
        for fac_num,faco in enumerate(fresids.T):
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(1,1,1)
            if gridopt: ax1.grid()
            for label in ax1.get_xticklabels():
                label.set_rotation(55)
            ax1.plot(dato[-len(faco):],faco,color='black')
            if glopt['title']: ax1.set_title('Estimated Residuals for FACTOR '+str(fac_num+1),fontsize=glopt['title_fs'])
            fig1.savefig('../graphs/'+modname+'/factors/'+'FACTOR_RESIDS_'+str(fac_num+1)+'.eps',bbox_inches='tight')     

    def plot_scree(self):
        from matplotlib.ticker import MultipleLocator, FormatStrFormatter
        majorLocator   = MultipleLocator(1)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(1)
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']        
        expvar = self.exp_var
        modname = self.modname
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'factors' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/factors')
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1,1,1)
        if gridopt: ax1.grid()
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)        
        if glopt['xlabel']: ax1.set_xlabel('Factors',fontsize=glopt['label_fs'])
        if glopt['ylabel']: ax1.set_ylabel('Explained Variation',fontsize=glopt['label_fs'])
        ax1.set_xlim(1,len(expvar))
        ax1.plot(numpy.arange(1,len(expvar)+1,1),expvar*100,'ko-')
        if glopt['title']: ax1.set_title('Scree Plot of explained variation due to Factors',fontsize=glopt['title_fs'])
        fig1.savefig('../graphs/'+modname+'/factors/'+'SCREE.eps',bbox_inches='tight')
     
        expvarcum=[]
        for i1,elem in enumerate(expvar):
            if i1 == 0:
                expvarcum.append(elem)
            else:
                expvarcum.append(elem+expvarcum[i1-1])
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1,1,1)
        if gridopt: ax1.grid()
        majorLocator   = MultipleLocator(1)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(1) 
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(10) 
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)        
        if glopt['xlabel']: ax1.set_xlabel('Factors',fontsize=glopt['label_fs'])
        if glopt['ylabel']: ax1.set_ylabel('Explained Cumulative Variation',fontsize=glopt['label_fs'])
        ax1.set_xlim(1,len(expvarcum))
        ax1.set_ylim(0,100)
        ax1.plot(numpy.arange(1,len(expvarcum)+1,1),[x*100.0 for x in expvarcum],'ko-')
        if glopt['title']: ax1.set_title('Scree Plot of explained cumulative variation due to Factors',fontsize=glopt['title_fs'])
        fig1.savefig('../graphs/'+modname+'/factors/'+'CUMSCREE.eps',bbox_inches='tight')

    
    def plot_vals(self):
        # Graph some vars and dump them in graphs directory
        dato = copy.deepcopy(self.dates)
        vnames = self.vnames
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']          
        matd = self.orig_tdata
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'series_fitted' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/series_fitted')
        yfmat = self.dynresdic['yfmat']
        resmat = self.dynresdic['ee_mat']
        for col,name in enumerate(vnames):
            fig1 = plt.figure()
            pyl.subplots_adjust(hspace=0.075)
            ax1 = fig1.add_subplot(311)
            if gridopt: ax1.grid()
            ax1.plot(matd[-yfmat[:,col].shape[0]:,col],color='black')
            labelo = ax1.get_xticklabels()
            pyl.setp(labelo, visible=False)            
            if glopt['title'] and not glopt['title_top_only']: ax1.set_title('Actual Values for '+name,fontsize=glopt['title_fs'])
            elif glopt['title'] and glopt['title_top_only']: ax1.set_title('Actual, Fitted and Residual Values for '+name,fontsize=glopt['title_fs'])            
            ax2 = fig1.add_subplot(312)             
            if gridopt: ax2.grid()
            ax2.plot(yfmat[:,col],color='black')
            labelo = ax2.get_xticklabels()
            pyl.setp(labelo, visible=False)            
            if glopt['title'] and not glopt['title_top_only']: ax2.set_title('Fitted Values for '+name,fontsize=glopt['title_fs'])
            ax3 = fig1.add_subplot(313)
            if gridopt: ax3.grid()
            for label in ax3.get_xticklabels():
                label.set_rotation(55)
            ax3.plot(dato[-len(resmat[:,col]):],resmat[:,col],color='black')
            if glopt['title'] and not glopt['title_top_only']: ax3.set_title('Residuals for '+name,fontsize=glopt['title_fs'])           
            fig1.savefig('../graphs/'+modname+'/series_fitted/'+name+'.eps',bbox_inches='tight')
        plt.close('all')

    def plot_vals2(self):
        # Graph some vars and dump them in graphs directory, here the original, transformed and standardized
        dato = copy.deepcopy(self.dates)
        vnames = self.vnames
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']          
        data = copy.deepcopy(self.orig_data[self.max_shift:,:])
        matd = copy.deepcopy(self.orig_tdata)
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'series' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/series')
        for col,name in enumerate(vnames):
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(211)
            pyl.subplots_adjust(hspace=0.075)
            if gridopt: ax1.grid()
            ax1.plot(dato[-len(data[:,col]):],data[:,col],color='black')
            labelo = ax1.get_xticklabels()
            pyl.setp(labelo, visible=False)
            if glopt['title'] and not glopt['title_top_only']: ax1.set_title('Raw data for '+name,fontsize=glopt['title_fs'])
            elif glopt['title'] and glopt['title_top_only']: ax1.set_title('Raw and transformed data for '+name,fontsize=glopt['title_fs'])
            ax2 = fig1.add_subplot(212)
            if gridopt: ax2.grid()
            for label in ax2.get_xticklabels():
                label.set_rotation(55)
            ax2.plot(dato[-len(matd[:,col]):],matd[:,col],color='black')
            if glopt['title'] and not glopt['title_top_only']: ax2.set_title('Transformed data for '+name,fontsize=glopt['title_fs'])            
            fig1.savefig('../graphs/'+modname+'/'+'series/'+name+'.eps',bbox_inches='tight')
        plt.close('all')

    def plot_irfs(self):
        # Graph some IRFs (Cholesky)
        irf_options = self.confdic['graph_options']['irfs']
        outer_colour = irf_options['outer_colour']
        inner_colour = irf_options['inner_colour']
        line_colour = irf_options['line_colour']
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']          
        vnames = self.vnames
        pnames = self.pnames
        use_svnames = self.confdic['use_svnames']
        svnames = self.confdic['svnames']
        irfp = self.confdic['irfp']
        cols = self.bettaf_mat.shape[0]
        dfacs = self.bettaf_mat.shape[1]
        phisc = self.phisc
        if self.confdic['bootstrap']:
            confint = self.confint
        loglist = self.confdic['loglist']
        translog = self.confdic['translog']
        tcode_dic = copy.deepcopy(self.tcode_dic)
        hzline = self.confdic['hzline']
        rescale = self.confdic['rescale']
        irf_anot = self.confdic['irf_anot']
        hzlwidth = self.confdic['hzline_width']
        bootstrap = self.confdic['bootstrap']
        kbsinbs = self.confdic['killian_bsinbs']
        irfident = self.confdic['irf_ident_graph']
        if kbsinbs: phisc_uncorr = self.phisc_uncorr
        graph_biased = self.confdic['graph_options']['irfs']['biased_line']
        graph_biased_type = self.confdic['graph_options']['irfs']['biased_line_type']
        param_bootstrap = self.confdic['param_bootstrap']
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'shock_irfs' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/shock_irfs')
        signif = self.confdic['signif']
        for shockv,sname in enumerate(range(dfacs)):
            #fig1 = plt.figure()
            if irfident and shockv != self.sfacs-1: continue
            for respv,rname in enumerate(vnames):
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(1,1,1)
                if glopt['xlabel']: ax2.set_xlabel('Time periods after shock',fontsize=glopt['label_fs'])
                plt.xlim((0,irfp))
                #ax = fig1.add_subplot(cols,1,respv+1)
                #ax.set_xlabel('Lags')
                plt.xlim((0,irfp))
                if tcode_dic[rname] in [4,5,6,7,8] and rescale and translog:
                    #ax.set_ylabel('\% Deviation')
                    if glopt['ylabel']: ax2.set_ylabel('Percentage Deviation',fontsize=glopt['label_fs'])
                elif tcode_dic[rname] not in [4,5,6,7,8] and rescale:
                    if glopt['ylabel']: ax2.set_ylabel('Absolute Deviation',fontsize=glopt['label_fs'])
                else:
                    #ax.set_ylabel('Abs. Deviation')
                    if glopt['ylabel']: ax2.set_ylabel('Normalized Deviation',fontsize=glopt['label_fs'])
                if hzline:
                    ax2.hlines(numpy.array([0,]*irfp),xmin=0,xmax=irfp,linestyles='solid',linewidth=hzlwidth,colors='k')
                    #ax.hlines(numpy.array([0,]*irfp),xmin=0,xmax=irfp,linestyles='solid',linewidth=hzlwidth,colors='k')
                if self.tcode_dic[rname] in [4,5,6,7,8] and translog:
                    percen = 100.0
                else:
                    percen = 1.0
                # When using the non-parametric bootstrap the central (mean) IR has to be the average of all IRs
                if not bootstrap:
                    ax2.plot(self.phisc[:,respv,shockv].flatten()*percen,color=line_colour)
                elif bootstrap:
                    ax2.plot(self.phisc_kcorr[:,respv,shockv].flatten()*percen,color=line_colour)
                if kbsinbs and graph_biased:
                    ax2.plot(self.phisc_uncorr[:,respv,shockv].flatten()*percen,graph_biased_type)
                if gridopt:ax2.grid()
                #ax.plot(phisc[:irfp+1,respv,shockv].flatten()*percen,color='black')
                #ax.grid()
                if bootstrap:
                    if type(signif) == type(10.0):
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
                        #ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*shocko,confint[1,:,respv,shockv].flatten()*percen, facecolor='#D3B0B0', edgecolor='#D3B0B0')
                    elif signif == 'mix':
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[1,:,respv,shockv].flatten()*percen,confint[2,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[2,:,respv,shockv].flatten()*percen,confint[3,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
                        #ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor='#EED6D6', edgecolor='#EED6D6')
                        #ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[1,:,respv,shockv].flatten()*percen,confint[2,:,respv,shockv].flatten()*percen, facecolor='#D3B0B0', edgecolor='#D3B0B0')
                        #ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[2,:,respv,shockv].flatten()*percen,confint[3,:,respv,shockv].flatten()*percen, facecolor='#EED6D6', edgecolor='#EED6D6')
                if self.confdic['use_bbe_ident'] and (sname+1) == self.fdata_adj.shape[1]:
                    if glopt['title']: ax2.set_title('Response of '+pnames[rname]+' to shock on '+self.confdic['separate_vars'][0],fontsize=glopt['title_fs'])
                else:
                    if glopt['title']: ax2.set_title('Response of '+pnames[rname]+' to shock on FACTOR '+str(sname+1),fontsize=glopt['title_fs'])        
                #ax.set_title(r'$\textrm{Response of '+pnames[rname]+' to shock on '+pnames[sname]+'}$')
                if irf_anot:
                    sigli = [cmp(x,0) for x in self.phisc[:20,respv,shockv].flatten()]
                    posc = sigli.count(1)
                    nevc = sigli.count(-1)
                    if nevc > posc: extr = 'min'
                    elif posc > nevc: extr = 'max'
                    elif posc == nevc: extr = 'max'
                    if extr == 'max': val = max(self.phisc[:20,respv,shockv].flatten())
                    elif extr == 'min': val = min(self.phisc[:20,respv,shockv].flatten())        
                    pos = list(self.phisc[:20,respv,shockv].flatten()).index(val)
                    if tcode_dic[rname] in [4,5,6,7,8]:
                        anval = numpy.round(self.trans_dic[rname]['std']*val*100.0,3)
                    else:
                        anval = numpy.round(self.trans_dic[rname]['std']*val,3)        
                    ax2.annotate(str(anval), xy=(pos, val), xytext=(-15,15), 
                                textcoords='offset points', ha='center', va='bottom', color='k')                
                fig2.savefig('../graphs/'+modname+'/'+'shock_irfs/irf_'+str(sname)+'_'+rname+'.eps',bbox_inches='tight')
                plt.close()
        plt.close('all')
            #fig1.savefig('../graphs/'+modname+'/'+'irf_'+str(sname)+'.png')
               
    def plot_facs_irfs(self):
        # Graph some IRFs (Cholesky)
        vnames = self.vnames
        pnames = self.pnames
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']          
        irfp = self.confdic['irfp']
        cols = self.bettaf_mat.shape[0]
        phisc = copy.deepcopy(self.phis_compta2c)
        dfacs = phisc[0].shape[0]
        if self.confdic['bootstrap']:
            confint = copy.deepcopy(self.confintf)
        loglist = self.confdic['loglist']
        translog = self.confdic['translog']
        hzline = self.confdic['hzline']
        bootstrap = self.confdic['bootstrap']
        param_bootstrap = self.confdic['param_bootstrap']
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'factor_irfs' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/factor_irfs')
        signif = self.confdic['signif']
        for shockv,sname in enumerate(range(dfacs)):
            for respv,rname in enumerate(range(dfacs)):
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(1,1,1)
                if glopt['xlabel']: ax2.set_xlabel('Time periods after shock')
                if rname in loglist and translog:
                    if glopt['ylabel']: ax2.set_ylabel('\% Deviation')
                else:
                    if glopt['ylabel']: ax2.set_ylabel('Abs. Deviation')
                if hzline:
                    ax2.hlines(numpy.array([0,]*irfp),xmin=0,xmax=irfp,linestyles='solid',linewidth=0.2,colors='k')
                if rname in loglist and translog:
                    percen = 100.0
                else:
                    percen = 1.0
                # When using the non-parametric bootstrap the central (mean) IR has to be the average of all IRs
                if not bootstrap:
                    ax2.plot(phisc[:irfp+1,respv,shockv].flatten()*percen,color='black')
                elif bootstrap and not param_bootstrap:
                    ax2.plot(numpy.mean(self.barrayf,axis=0)[:irfp+1,respv,shockv]*percen,color='black')
                elif bootstrap and param_bootstrap:
                    ax2.plot(phisc[:irfp+1,respv,shockv].flatten()*percen,color='black')
                plt.xlim((0,irfp))                
                if gridopt: ax2.grid()
                if bootstrap:
                    if type(signif) == type(10.0):
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor='#D3B0B0', edgecolor='#D3B0B0')
                    elif signif == 'mix':
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor='#EED6D6', edgecolor='#EED6D6')
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[1,:,respv,shockv].flatten()*percen,confint[2,:,respv,shockv].flatten()*percen, facecolor='#D3B0B0', edgecolor='#D3B0B0')
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[2,:,respv,shockv].flatten()*percen,confint[3,:,respv,shockv].flatten()*percen, facecolor='#EED6D6', edgecolor='#EED6D6')
                fig2.savefig('../graphs/'+modname+'/'+'factor_irfs/irf_'+str(sname)+'_'+str(rname)+'.eps',bbox_inches='tight')
                plt.close()
        plt.close('all')


    def plot_xdata_irfs(self):
        # Graph some IRFs (Cholesky)
        vnames = self.vnames
        pnames = self.pnames
        glopt = self.confdic['graph_options']['labels']
        gridopt = self.confdic['graph_options']['grid']          
        irfp = self.confdic['xirfp']
        cols = self.bettax_mat.shape[0]
        phisc = copy.deepcopy(self.phis_compta)
        dfacs = phisc[0].shape[0]
        if self.confdic['bootstrap']:
            confint = copy.deepcopy(self.confintf)
        loglist = self.confdic['loglist']
        translog = self.confdic['translog']
        hzline = self.confdic['hzline']
        bootstrap = self.confdic['bootstrap']
        param_bootstrap = self.confdic['param_bootstrap']
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'xdata_irfs' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/xdata_irfs')
        signif = self.confdic['signif']
        for i1,vnamo in enumerate(vnames):
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(1,1,1)
            if glopt['xlabel']: ax2.set_xlabel('Lags')
            plt.xlim((0,irfp))
            plt.xlim((0,irfp))
            if vnamo in loglist and translog:
                if glopt['ylabel']: ax2.set_ylabel('\% Deviation')
            else:
                if glopt['ylabel']: ax2.set_ylabel('Abs. Deviation')
            if hzline:
                ax2.hlines(numpy.array([0,]*irfp),xmin=0,xmax=irfp,linestyles='solid',linewidth=0.2,colors='k')
            if vnamo in loglist and translog:
                percen = 100.0
            else:
                percen = 1.0
            # When using the non-parametric bootstrap the central (mean) IR has to be the average of all IRs
            '''
            if not bootstrap:
            '''
            ax2.plot(phisc[:irfp+1,i1].flatten()*percen,color='black')
            '''
            elif bootstrap and not param_bootstrap:
                ax2.plot(numpy.mean(self.barrayf,axis=0)[:irfp+1,vnamo]*percen,color='black')
            elif bootstrap and param_bootstrap:
                ax2.plot(phisc[:irfp+1,vnamo].flatten()*percen,color='black')
            '''
            if gridopt: ax2.grid()
            '''
            if bootstrap:
                if type(signif) == type(10.0):
                    ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor='#D3B0B0', edgecolor='#D3B0B0')
                elif signif == 'mix':
                    ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor='#EED6D6', edgecolor='#EED6D6')
                    ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[1,:,respv,shockv].flatten()*percen,confint[2,:,respv,shockv].flatten()*percen, facecolor='#D3B0B0', edgecolor='#D3B0B0')
                    ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[2,:,respv,shockv].flatten()*percen,confint[3,:,respv,shockv].flatten()*percen, facecolor='#EED6D6', edgecolor='#EED6D6')
            '''        
            fig2.savefig('../graphs/'+modname+'/'+'xdata_irfs/irf_'+vnamo+'.eps',bbox_inches='tight')
            plt.close()
        plt.close('all')

################################################# HELPER METHODS (TOOLS) ##################################################
###########################################################################################################################
 

    # This is a function which will write the computed companion matrix into a csv file
    def writeArray(self,inarray=None,fname=None):
        modname = self.confdic['modname']
        # Check if folder exists
        if 'matrices' not in os.listdir('../graphs/'+modname):
            os.mkdir('../graphs/'+modname+'/matrices')
        if inarray == None: varcomp = self.varcomp
        else: varcomp = inarray
        rows = []
        for row in range(0,varcomp.shape[0],1):
            tmp_str = ''
            tmp_str = tmp_str+str(varcomp[row][0])+';'
            for numo in varcomp[row][1:]:
                tmp_str = tmp_str+str(numpy.round(numo,8))+';'
            rows.append(tmp_str[:-1])
        if fname == None:
            if modname+'_compmat.csv' in os.listdir(os.getcwd()):
                os.remove(os.path.join(os.getcwd(),modname+'_compmat.csv'))
            csvfile = open(modname+'_compmat.csv',"wb")
        else:
            if fname+'.csv' in os.listdir(os.getcwd()):
                os.remove(os.path.join(os.getcwd(),fname+'.csv'))
            csvfile = open(os.path.join(os.getcwd(),'../graphs/'+modname+'/'+fname+'.csv'),"wb")            
        for row in rows:
            csvfile.writelines(row)
            csvfile.write('\n')
        csvfile.close()

     # This is a function which will write the found factors and the original data into a csv file
    def writeFData(self):
        modname = self.confdic['modname']
        # Check if folder exists
        if 'matrices' not in os.listdir('../graphs/'+modname):
            os.mkdir('../graphs/'+modname+'/matrices')        
        fdata = self.fdata_adj
        vnames = self.vnames
        fnames = []
        for fname in range(0,fdata.shape[1],1):
            fnames.append('Factor'+str(fname+1))
        all_names = fnames + vnames
        min_rows = fdata.shape[0]
        data = self.tdata
        max_rows = data.shape[0]
        diff_rows = max_rows - min_rows
        data = data[diff_rows:,:]
        all_data = numpy.hstack((fdata,data))
        header = ''
        for vname in all_names:
            header = header+'"'+vname+'";'
        header = header[:-1]
        rows = []
        for row in range(0,all_data.shape[0],1):
            tmp_str = ''
            tmp_str = tmp_str+str(all_data[row][0])+';'
            for numo in all_data[row][1:]:
                tmp_str = tmp_str+str(numpy.round(numo,8))+';'
            rows.append(tmp_str[:-1])
        if modname+'.csv' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modname+'.csv'))
        csvfile = open(modname+'.csv',"wb")
        csvfile.writelines(header)
        csvfile.write('\n')
        for row in rows:
            csvfile.writelines(row)
            csvfile.write('\n')
        csvfile.close()
        

    # Method used to clean all directories of their contents    
    def clean_dirs(self,modpath=None):
        if modpath == None: base_path = '../graphs/'+self.modname
        else: base_path = modpath
        if 'xdata_irfs' in os.listdir(base_path):
            seriesxx = glob.glob(base_path+'/xdata_irfs/*')
            for f in seriesxx:
                os.remove(f)
        if 'factor_irfs' in os.listdir(base_path):
            seriesff = glob.glob(base_path+'/factor_irfs/*')
            for f in seriesff:
                os.remove(f)
        if 'shock_irfs' in os.listdir(base_path):
            seriesf = glob.glob(base_path+'/shock_irfs/*')
            for f in seriesf:
                os.remove(f)
        if 'series' in os.listdir(base_path):
            seriess = glob.glob(base_path+'/series/*')
            for f in seriess:
                os.remove(f)
        if 'series_fitted' in os.listdir(base_path):
            seriesfit = glob.glob(base_path+'/series_fitted/*')
            for f in seriesfit:
                os.remove(f)
        if 'factors' in os.listdir(base_path):
            seriesfac = glob.glob(base_path+'/factors/*')
            # Remove the subfolder from glob list as we don't want to delete this
            seriesfac.remove('../graphs/'+self.modname+'/factors/unadj')
            for f in seriesfac:
                os.remove(f)
        if 'matrices' in os.listdir(base_path):
            seriesmat = glob.glob(base_path+'/matrices/*')
            for f in seriesmat:
                os.remove(f)