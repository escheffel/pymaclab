'''
.. module:: var
   :platform: Linux
   :synopsis: The var module contains the VAR class for estimating and doing further work with Vector Autoregressions commonly
              used in applied macroeconometrics. It supports advanced methods such as bootstrapping confidence intervals including
              Killian's boostrap-after-bootstrap small-sample bias correction. Also CPU-intensive methods such as the bootstrap can
              be computed using Parallel Python to exploit multi-core CPUs. Pretty plotting methods are also included which depend
              on matplotlib.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
import numpy
import scipy
import scipy.stats
# We don't need this anywhere and it should also be replaced by Pandas in the future
#from scikits import timeseries as ts
import matplotlib as mpl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
from matplotlib import pyplot as plt
from matplotlib import pylab as pyl
import copy
import datetime
import pickle
import pp
import time
import os

# Re-factored imports
from ..dattrans.transmeth import stdX, transX
from .common import genXX, compbetta

# Switch off runtime warnings here
import warnings
warnings.filterwarnings("ignore")


class VAR:

    def __init__(self,data=None,vnames=None,pnames=None,svnames=None,irfs=True,boot=True,plot=True,conf=None,mesg=False):
        print '#####################################################'
        stepo=0
        stepo+=1
        self._mesg = mesg
        if self._mesg: print str(stepo)+') Initialising VAR'
        self.confdic = conf
        params = {'axes.labelsize':conf['graph_options']['labels']['label_fs'],
                  'text.fontsize': conf['graph_options']['labels']['text_fs'],
                  'legend.fontsize': conf['graph_options']['labels']['legend_fs'],
                  'xtick.labelsize': conf['graph_options']['labels']['xticks_fs'],
                  'ytick.labelsize': conf['graph_options']['labels']['yticks_fs'],
                  'text.usetex': conf['graph_options']['labels']['use_tex'],
                  'lines.linewidth': conf['graph_options']['lines']['width']}
        pyl.rcParams.update(params)        
        self.vnames = vnames
        self.pnames = pnames
        self.svnames = svnames
        self.const_term = conf['const_term']
        self.nlags = conf['nlags']
        self.freq = conf['freq']
        self.transperc_dic = conf['transperc_dic']
        self.tcode_dic = copy.deepcopy(conf['tcode'])
        self.data = data
        self.ncols = data.shape[1]
        self.nrows = data.shape[0]
        # Create irf_dic for graphs
        self.irf_graphs_dic = {}
        for sname in vnames:
            self.irf_graphs_dic[sname] = {}
            for rname in vnames:
                self.irf_graphs_dic[sname][rname] = None
        ############################################################
        if self._mesg:
            print "This VAR model's name is: "+self.confdic['modname']
            print 'The dataset used is: '+self.confdic['datafile']
            print 'The chosen autoregressive order is: '+str(self.nlags)
            print 'The total number of variables in the system is: '+str(self.ncols)
            print 'The total number of observations before estimation is: '+str(self.nrows)
            if conf['const_term'] == None:
                print 'The model was estimated using neither a constant nor a time-trend'
            elif conf['const_term'] == 'cc':
                print 'The model was estimated using only a constant'
            elif conf['const_term'] == 'tt':
                print 'The model was estimated using only a time-trend'
            elif conf['const_term'] == 'ct':
                print 'The model was estimated using both a constant and a time-trend'
            elif conf['const_term'] == 'ctt':
                print 'The model was estimated using both a constant, a linear, and a quadratic time trend'
            if conf['translog']:
                print "The model's impulse responses are in % changes for vars in logs"
            elif not conf['translog']:
                print "The model's impulse responses are in abs. log dev. for vars in logs"
            if conf['use_svnames']:
                print "The model's impulse responses are based on a customized shock matrix"
            elif not conf['use_svnames']:
                print "The model's impulse responses are based on the identity shock matrix"
            print 'The frequency of the data is: '+str(self.freq)
            print '#####################################################'
        stepo+=1
        if self._mesg: print str(stepo)+') Standardising data'
        # Transform raw data according to Stock and Watson's trans code encoding system
        if self._mesg: print 'Now TRANSFORMING data according to Stock and Watson transformation encoding...'
        self = transX(self,data=self.data)
        if self.confdic['trans_data']:
            # Standardize the data if wanted, also store means and standard errors
            if self._mesg: print 'Now STANDARDIZING data ahead of estimation procedure...'
            self = stdX(self,data=self.data,standard=True)
            self.data = copy.deepcopy(self.tdata)
        else:
            self = stdX(self,data=self.data,standard=False)
        if self._mesg: print '#####################################################'
        stepo+=1
        if self._mesg: print str(stepo)+') Estimate coefficient matrix of the VAR'
        self = genXX(self)
        compbetta(self)
        if self.confdic['const_term'] != None:
            self.compbettax()
        if self._mesg: print '#####################################################'
        stepo+=1
        if self._mesg: print str(stepo)+') Creating companion matrices'
        self.mkCompMat()
        self.mkfitres()
        self.mksigma()
        if irfs:
            if self._mesg: print 'Now generating impulse responses...'
            self.mkphis()       
        if irfs and boot:
            addstr = ''
            if self.confdic['multicpu']: addstr = '(parallel execution)'
            else: addstr = '(serial execution)'
            if self._mesg: print 'Now generating bootstrap confidence intervals'+addstr+'...'
            t0 = time.time()            
            if not self.confdic['multicpu']: self.mkboot()
            elif self.confdic['multicpu']: self.mkboot_pp()
            if self._mesg: print round(time.time() - t0,2), "secs used to execute..."
        if self.confdic['do_vdc']:
            vdcp = self.confdic['vdcp']
            self.genFEVD(self.data,self.bettasm,vdcp) 
        if plot:
            if self._mesg: print 'Now plotting and saving graphs of actual, fitted and residual data for all series...'
            self.plot_vals()
        if irfs:
            if self._mesg: print 'Now plotting and saving graphs of impulse responses...'
            self.plot_irfs()
        if self._mesg: print '*All done !*'
   

    # Computes coefficient matrix but with a column of ones in the Y matrix.
    # This is useful for when there is a constant in the coefficient matrix
    def compbettax(self,matd=None,nmatd=None,smbetta=True,const='standard',func=False):
        # Get all the needed variables from instance
        if matd == None: matd = self.data
        if nmatd == None: nmatd = self.nmatd
        if const == 'standard': const_term = self.confdic['const_term']
        cols = matd.shape[1]
        if const_term == None:
            nlags = int(nmatd.shape[1]/cols)
        elif const_term == 'cc' or const_term == 'tt':
            nlags = int((nmatd.shape[1]-1)/cols)
            # Add column of ones or time index to the ymat
            if const_term == 'cc': matd = numpy.hstack((numpy.ones((matd.shape[0],1)),matd))
            elif const_term == 'tt': matd = numpy.hstack((numpy.arange(0,matd.shape[0],1.0).reshape(len(numpy.arange(matd.shape[0])),1),matd))
        elif const_term == 'ct':
            # Add column of ones and column of time index to ymat, so two extra columns
            nlags = int((nmatd.shape[1]-2)/cols)
            const = numpy.ones((matd.shape[0],1))
            timet = numpy.arange(0,matd.shape[0],1.0).reshape(len(numpy.arange(matd.shape[0])),1)
            matd = numpy.hstack((const,timet,matd))
        elif const_term == 'ctt':
            # Add column of ones and column of time index to ymat, so two extra columns
            nlags = int((nmatd.shape[1]-3)/cols)
            const = numpy.ones((matd.shape[0],1))
            timet = numpy.arange(0,matd.shape[0],1.0).reshape(len(numpy.arange(matd.shape[0])),1)
            timet2 = timet**2
            matd = numpy.hstack((const,timet,timet2,matd))
        fmatd = nmatd
        XX = numpy.dot(fmatd.T,fmatd)
        XXi = numpy.linalg.inv(XX)
        if const_term == None:
            betta = numpy.zeros((cols,nlags*cols))
        elif const_term == 'tt' or const_term == 'cc':
            betta = numpy.zeros((cols+1,1+nlags*cols))
        elif const_term == 'ct':
            betta = numpy.zeros((cols+2,2+nlags*cols))
        elif const_term == 'ctt':
            betta = numpy.zeros((cols+3,3+nlags*cols))
        if const_term == 'cc' or const_term == 'tt': fcols = cols+1
        elif const_term == 'ct': fcols = cols+2
        elif const_term == 'ctt': fcols = cols+3
        elif const_term == None: fcols = cols
        for vari in range(0,fcols,1):
            ymat=matd[nlags:,vari]
            Xy = numpy.dot(fmatd.T,ymat)
            betta[vari,:] = numpy.dot(XXi,Xy)
        if smbetta:
            # Produce the betta format from statsmodels library so we can use some of their code
            # Exclude the trend or constant from bettasm, as opposed to in betta
            bettasm = numpy.zeros((nlags,cols,cols))
            for lago in range(0,nlags,1):
                for varos in range(0,cols,1):
                    if const_term == 'tt' or const_term == 'cc':
                        bettasm[lago][varos,:] = betta[:,1:][varos][cols*lago:cols*(lago+1)]
                    elif const_term == 'ct':
                        bettasm[lago][varos,:] = betta[:,2:][varos][cols*lago:cols*(lago+1)]
                    elif const_term == 'ctt':
                        bettasm[lago][varos,:] = betta[:,3:][varos][cols*lago:cols*(lago+1)]
                    elif const_term == None:
                        bettasm[lago][varos,:] = betta[varos][cols*lago:cols*(lago+1)]
        if not func:
            if smbetta:
                self.betta_one_sm = bettasm
                self.betta_one = betta
            elif not smbetta:
                self.betta_one = betta
        elif func:
            if smbetta:
                return bettasm
            elif not smbetta:
                return betta

            
    def mkCompMat(self,betta=None,betta_one=None,func=False):
        if betta == None: betta = copy.deepcopy(self.betta)
        if betta_one == None and 'betta_one' in dir(self):
            betta_one = copy.deepcopy(self.betta_one)
        cols = self.ncols
        lags = self.nlags
        const = self.confdic['const_term']
        if const == None:
            compmat = numpy.zeros((cols*lags,cols*lags))
            eyem = numpy.eye(cols)
            compmat[:cols,:] = betta
            for lago in range(1,lags,1):
                compmat[cols*lago:cols*(lago+1),cols*(lago-1):cols*lago] = eyem
        elif const == 'cc' or const == 'tt':
            compmat = numpy.zeros((1+cols*lags,1+cols*lags))
            compmat[0,0] = 1.0
            compmat[0,1:] = 0.0
            compmat[1:,0] = 0.0
            eyem = numpy.eye(cols)
            compmat[1:cols+1,:] = betta_one[1:,:]
            for lago in range(1,lags,1):
                compmat[1+cols*lago:cols*(lago+1)+1,cols*(lago-1):cols*lago] = eyem
        elif const == 'ct':
            compmat = numpy.zeros((2+cols*lags,2+cols*lags))
            compmat[0,0] = 1.0
            compmat[1,1] = 1.0
            compmat[0,2:] = 0.0
            compmat[2:,0] = 0.0
            eyem = numpy.eye(cols)
            compmat[2:cols+2,:] = betta_one[2:,:]
            for lago in range(1,lags,1):
                compmat[2+cols*lago:cols*(lago+1)+2,cols*(lago-1):cols*lago] = eyem
        elif const == 'ctt':
            compmat = numpy.zeros((3+cols*lags,3+cols*lags))
            compmat[0,0] = 1.0
            compmat[1,1] = 1.0
            compmat[0,3:] = 0.0
            compmat[3:,0] = 0.0
            eyem = numpy.eye(cols)
            compmat[3:cols+3,:] = betta_one[3:,:]
            for lago in range(1,lags,1):
                compmat[3+cols*lago:cols*(lago+1)+3,cols*(lago-1):cols*lago] = eyem
        if not func:
            self.compmat = copy.deepcopy(compmat)
            # Also create a companion matrix whic strips out the constant and time trend effects, only focus on short-run behaviour
            if const == None:
                self.compmatr = compmat
            elif const == 'cc' or const == 'tt':
                self.compmatr = compmat[1:,1:]
                self.compmatr[cols:,:] = 0.0
                for lago in range(1,lags,1):
                    self.compmatr[cols*lago:cols*(lago+1),cols*(lago-1):cols*lago] = eyem             
            elif const == 'ct':
                self.compmatr = compmat[2:,2:]
                self.compmatr[cols:,:] = 0.0
                for lago in range(1,lags,1):
                    self.compmatr[cols*lago:cols*(lago+1),cols*(lago-1):cols*lago] = eyem
            elif const == 'ctt':
                self.compmatr = compmat[3:,3:]
                self.compmatr[cols:,:] = 0.0
                for lago in range(1,lags,1):
                    self.compmatr[cols*lago:cols*(lago+1),cols*(lago-1):cols*lago] = eyem
        elif func:
            # Also create a companion matrix whic strips out the constant and time trend effects, only focus on short-run behaviour
            if const == None:
                compmatr = compmat
            elif const == 'cc' or const == 'tt':
                compmatr = compmat[1:,1:]
                compmatr[cols:,:] = 0.0
                for lago in range(1,lags,1):
                    compmatr[cols*lago:cols*(lago+1),cols*(lago-1):cols*lago] = eyem                            
            elif const == 'ct':
                compmatr = compmat[2:,2:]
                compmatr[cols:,:] = 0.0
                for lago in range(1,lags,1):
                    compmatr[cols*lago:cols*(lago+1),cols*(lago-1):cols*lago] = eyem
            elif const == 'ctt':
                compmatr = compmat[3:,3:]
                compmatr[cols:,:] = 0.0
                for lago in range(1,lags,1):
                    compmatr[cols*lago:cols*(lago+1),cols*(lago-1):cols*lago] = eyem
            return compmatr

    def mkfitres(self,matd=None,nmatd=None,betta=None,nlags=None,func=False):
        # Create fitted values and residuals
        if matd == None: matd = self.data
        if nmatd == None: nmatd = self.nmatd
        if betta == None: betta = self.betta
        if nlags == None: nlags = self.nlags
        yfmat = numpy.dot(nmatd,betta.T)
        resmat = matd[nlags:]-yfmat
        covmat = numpy.cov(resmat.T)
        if not func:
            self.yfmat = yfmat
            self.resmat = resmat
            self.covmat = covmat
        elif func:
            return yfmat,resmat,covmat
    
    def mksigma(self,resmat=None,func=False):
        if resmat == None: resmat = self.resmat
        vnames = self.vnames
        svnames = self.svnames
        transd = self.tcode_dic
        # Create the variance-covariance matrix of the residuals
        vcmat = numpy.cov(resmat.T)
        # Create the Cholesky decomposition of vcmat
        cholmat = numpy.linalg.cholesky(vcmat)
        # This is the step which scales up the initial shock, if so configured in svname
        eyem = numpy.eye(cholmat.shape[0])
        if self.confdic['use_svnames']:
            for i1,name in enumerate(vnames):
                if transd[name] in [4,5,6,7,8]:
                    if svnames[name] != 1.0:
                        adjo = numpy.log(float(svnames[name]))/float(cholmat[i1,i1])
                        eyem[i1,i1] = eyem[i1,i1]*adjo
                    else: eyem[i1,i1] = 1.0
                else:
                    if svnames[name] != 1.0:
                        adjo = float(svnames[name])/float(cholmat[i1,i1])
                        eyem[i1,i1] = eyem[i1,i1]*adjo
                    else:
                        eyem[i1,i1] = 1.0
        cholmat_scaled = numpy.dot(cholmat,eyem)
        if not func:
            self.vcmat = vcmat
            self.cholmat_scaled = cholmat_scaled
            self.cholmat = cholmat_scaled
            self.cholmat_unscaled = cholmat
            self.eyem = eyem
        elif func:
            return vcmat,cholmat_scaled,eyem

    def mkphis(self,bettasm=None,betta=None,compmatr=None,maxphi=None,cholmat=None,cholimp=None,func=False):
        if bettasm == None: bettasm = self.bettasm
        if maxphi == None: maxphi = self.confdic['maxphi']
        if cholimp == None : cholimp = self.confdic['cholimp']
        if cholmat == None: cholmat = self.cholmat
        if compmatr == None: compmatr = copy.deepcopy(self.compmatr)
        if betta == None: betta = copy.deepcopy(self.betta)
        nlags = self.nlags
        ncols = self.ncols
        const_term = self.confdic['const_term']
        bettar = betta[:,-ncols*nlags:]
        loglist = self.confdic['loglist']
        vnames = self.vnames
        const = self.confdic['const_term']
        k = copy.deepcopy(self.ncols)
        p = copy.deepcopy(self.nlags)
        transd = self.trans_dic
        # Now try the version using the companion matrix
        betta_comp = compmatr
        phis_comp = numpy.zeros((maxphi+1,k,k))
        bigmat_temp = copy.deepcopy(betta_comp)
        phis_comp[0] = numpy.eye(k)       
        # recursively compute Phi matrices, but using the first-order companion matrix
        for i in xrange(1, maxphi + 1):
            '''
            phis_comp[i] = bigmat_temp[:k,:k]
            bigmat_temp = numpy.dot(bigmat_temp,betta_comp)
            '''
            for j in xrange(1,i+1):
                if j > nlags: break
                phis_comp[i] += numpy.dot(phis_comp[i-j],bettar[:,(j-1)*k:j*k])
        self.phis_comp = copy.deepcopy(phis_comp)
        # Also create the orthogonalized equivalents
        if cholimp:
            phisc = numpy.array([list(numpy.dot(tmat, cholmat)) for tmat in phis_comp[1:]])
            phisc[0] = numpy.eye(phisc.shape[2])
            for i1,name in enumerate(vnames):
                if transd[name] in [4,5,6,7,8]:
                    phisc[0][i1,i1] = 0.01
        else:
            phisc = numpy.array([list(numpy.dot(tmat, cholmat)) for tmat in phis_comp])
        if not func:
            self.phis = copy.deepcopy(phis_comp)
            self.phisc = copy.deepcopy(phisc)
        elif func:
            return phis_comp,phisc

    def mkboot(self,bdraw=None,yfmat=None,resiarray=None):
        # Get all of the paramters
        if bdraw == None: bdraw = self.confdic['bdraw']
        maxphi = self.confdic['maxphi']
        signif = self.confdic['signif']
        irfp = self.confdic['irfp']
        if resiarray == None: resmat_n = copy.deepcopy(self.resmat)
        # Need to centre residuals prior to using them in the bootstrap procedure
        resmat_n = resmat_n - numpy.mean(resmat_n,axis=0)        
        if yfmat == None: yfmat = self.yfmat
        nlags = copy.deepcopy(self.nlags)
        nrows = copy.deepcopy(self.nrows)
        ncols = copy.deepcopy(self.ncols)
        matd = copy.deepcopy(self.data)
        const = self.confdic['const_term']
        shifto = 0
        if const == 'cc' or const == 'tt':
            shifto = 1
        elif const == 'ct':
            shifto = 2
        elif const == 'ctt':
            shifto = 3
        betta = copy.deepcopy(self.betta)

        # Produce bootstrap confidence intervals
        barray = []
        betarray = []
        rarray = []
        resarray = []
        nmatarray = []
        # Do the first bootstrap run a la Killian 1998 to reduce the small sample bias of the IRFs, if wanted
        if self.confdic['killian_bsinbs']:
            self.phis_uncorr = copy.deepcopy(self.phis)
            print "Calculating bias correction using Killian's(1998) bs-in-bs method"
            kill_array = []
            for bodraw in range(0,bdraw,1):
                init_vals = self.data[:nlags,:]
                matd_n = numpy.zeros((nrows,ncols))
                matd_n[:nlags,:] = copy.deepcopy(init_vals)                
                numpy.random.shuffle(resmat_n)
                for obso in xrange(0,(nrows-nlags),1):
                    for j in xrange(1,nlags+1):
                        # Producing the fitted values period-by-period
                        if const == None:
                            matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,(j-1)*ncols:j*ncols].T)
                        elif const == 'cc':
                            matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)
                        elif const == 'tt':
                            matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)+betta[:,0].T
                        elif const == 'ct':
                            matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,2+(j-1)*ncols:2+j*ncols].T)+betta[:,1].T                    
                    # Then finally adding the reshuffled residuals period-by-period
                    matd_n[nlags+obso,:] += resmat_n[obso,:]
                nmatd_n = genXX(self,matd=matd_n,func=True)
                betta_n = compbetta(self,matd=matd_n,nmatd=nmatd_n,smbetta=False,func=True)
                kill_array.append(betta_n)
            kill_array = numpy.array(kill_array)
            biasred_betta = numpy.mean(kill_array,axis=0)
            # Now calculated bias-corrected matrix of estimated coefficients and build alternative yfmat
            biasred_term = biasred_betta - self.betta
            # Make the bias in value and percentage available for inspection
            self.killbias = copy.deepcopy(biasred_term)
            self.killbiasperc = copy.deepcopy((biasred_term/self.betta)*100.0)
            self.betta_kcorr = copy.deepcopy(self.betta + biasred_term)
            betta_corrected = self.betta + biasred_term
            '''
            nmatd = copy.deepcopy(self.nmatd)
            yfmat = numpy.dot(nmatd,betta_corrected.T)
            resmat_n = matd[nlags:,:] - yfmat
            # Need to centre residuals prior to using them in the bootstrap procedure
            resmat_n = resmat_n - numpy.mean(resmat_n,axis=0)
            '''
            self.betta_uncorr = copy.deepcopy(self.betta)
            betta = copy.deepcopy(self.betta_kcorr)
            if 'betta_one' in dir(self):
                self.betta_one_uncorr = copy.deepcopy(self.betta_one)
                betta_one_kcorr = copy.deepcopy(self.betta_one)
                if self.confdic['const_term'] == 'cc' or self.confdic['const_term'] == 'tt':
                    betta_one_kcorr[1:,:] = copy.deepcopy(self.betta_kcorr)
                elif self.confdic['const_term'] == 'ct':
                    betta_one_kcorr[2:,:] = copy.deepcopy(self.betta_kcorr)
                elif self.confdic['const_term'] == 'ctt':
                    betta_one_kcorr[3:,:] = copy.deepcopy(self.betta_kcorr)
                compmatr = self.mkCompMat(betta=self.betta_kcorr,betta_one=betta_one_kcorr,func=True)
            else:
                compmatr = self.mkCompMat(betta=self.betta_kcorr,func=True)
            self.phisc_uncorr = copy.deepcopy(self.phisc)
            phis_kcorr,phisc_kcorr = self.mkphis(betta=self.betta_kcorr,compmatr=compmatr,func=True)
            self.phisc_kcorr = copy.deepcopy(phisc_kcorr)
            self.phis_kcorr = copy.deepcopy(phis_kcorr)
            self.phisc_uncorr = copy.deepcopy(self.phisc)
            self.phisc = copy.deepcopy(self.phisc_kcorr)     
            print 'Calculation of bias correction done...'
        print 'Now computing bootstraps'            
        for bodraw in range(0,bdraw,1):
            numpy.random.shuffle(resmat_n)
            init_vals = self.data[:nlags,:]
            matd_n = numpy.zeros((nrows,ncols))
            matd_n[:nlags,:] = copy.deepcopy(init_vals)
            for obso in xrange(0,(nrows-nlags),1):
                for j in xrange(1,nlags+1):
                    if const == None:
                        matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,(j-1)*ncols:j*ncols].T)
                    elif const == 'cc':
                        matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)
                    elif const == 'tt':
                        matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)+betta[:,0].T
                    elif const == 'ct':
                        matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,2+(j-1)*ncols:2+j*ncols].T)+betta[:,1].T
                matd_n[nlags+obso,:] += resmat_n[obso,:]
            rows_n,cols_n = matd_n.shape
            nmatd_n = genXX(self,matd=matd_n,func=True)
            betta_n = compbetta(self,matd=matd_n,nmatd=nmatd_n,smbetta=False,func=True)
            betta_one_n = self.compbettax(matd=matd_n,nmatd=nmatd_n,smbetta=False,func=True)
            compmatr_n = self.mkCompMat(betta=betta_n,betta_one=betta_one_n,func=True)
            # Use different names here as we don't want to overwrite the original resmat_n
            yfmat_nn,resmat_nn,covmat_nn = self.mkfitres(matd=matd_n,nmatd=nmatd_n,betta=betta_n,func=True)
            vcmat_nn,cholmat_nn,eyem_nn = self.mksigma(resmat=resmat_nn,func=True)
            # Use the original cholmat
            cholmat_n = copy.deepcopy(self.cholmat)
            phis_n,phisc_n = self.mkphis(compmatr=compmatr_n,betta=betta_n,maxphi=maxphi,cholmat=cholmat_n,func=True)
            barray.append(phisc_n)
        print 'Calculation of bootstraps done...'
        nmatarray = numpy.array(nmatarray)
        resarray = numpy.array(resarray)
        barray = numpy.array(barray)
        betarray = numpy.array(betarray)
        rarray = numpy.array(rarray)
        if type(signif) == type(10.0):
            lower = scipy.stats.scoreatpercentile(barray,per=((100.0-signif*100.0)/2.0))
            upper = scipy.stats.scoreatpercentile(barray,per=signif*100.0+((100.0-signif*100.0)/2.0))
            confint = [lower,upper]
        elif type(signif) == type([10.0,10.0]):
            lower_outer = scipy.stats.scoreatpercentile(barray,per=(100.0-signif[1]*100.0)/2.0)
            lower_inner = scipy.stats.scoreatpercentile(barray,per=(100.0-signif[0]*100.0)/2.0)
            upper_inner = scipy.stats.scoreatpercentile(barray,per=signif[0]*100.0+((100.0-signif[0]*100.0)/2.0))
            upper_outer = scipy.stats.scoreatpercentile(barray,per=signif[1]*100.0+((100.0-signif[1]*100.0)/2.0))
            confint = [lower_outer,lower_inner,upper_inner,upper_outer]
        self.rarray = rarray
        self.betarray = betarray
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
        betta = intup[5]
        nrows = self.nrows
        ncols = self.ncols
        nlags = self.nlags
        const = self.const_term
        resmat_n = copy.deepcopy(resmat)
        numpy.random.shuffle(resmat_n)
        init_vals = self.data[:nlags,:]
        matd_n = numpy.zeros((nrows,ncols))
        matd_n[:nlags,:] = copy.deepcopy(init_vals)
        numpy.random.shuffle(resmat_n)
        for obso in xrange(0,(nrows-nlags),1):
            for j in xrange(1,nlags+1):
                if const == None:
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,(j-1)*ncols:j*ncols].T)
                elif const == 'cc':
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)
                elif const == 'tt':
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)+betta[:,0].T
                elif const == 'ct':
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,2+(j-1)*ncols:2+j*ncols].T)+betta[:,1].T                    
            matd_n[nlags+obso,:] += resmat_n[obso,:]
        nmatd_n = pymaclab.stats.common.genXX(self,matd=matd_n,func=True)
        betta_n = pymaclab.stats.common.compbetta(self,matd=matd_n,nmatd=nmatd_n,smbetta=False,func=True)
        compbetta_n = self.mkCompMat(betta=betta_n,func=True)
        eigs = numpy.linalg.eigvals(compbetta_n)
        sbool = (numpy.abs(eigs) <= 1).all()
        if not sbool: print 'NOT STABLE !'
        return betta_n

    def refac_pp(self,intup):
        resmat = intup[0]
        yfmat = intup[1]
        matd = intup[2]
        nlags = intup[3]
        maxphi = intup[4]
        betta = intup[5]
        const = self.const_term
        nlags = self.nlags
        nrows = self.nrows
        ncols = self.ncols
        resmat_n = copy.deepcopy(resmat)
        numpy.random.shuffle(resmat_n)
        init_vals = self.data[:nlags,:]
        matd_n = numpy.zeros((nrows,ncols))
        matd_n[:nlags,:] = copy.deepcopy(init_vals)
        for obso in xrange(0,(nrows-nlags),1):
            for j in xrange(1,nlags+1):
                if const == None:
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,(j-1)*ncols:j*ncols].T)
                elif const == 'cc':
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)
                elif const == 'tt':
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,1+(j-1)*ncols:1+j*ncols].T)+betta[:,0].T
                elif const == 'ct':
                    matd_n[nlags+obso,:] += numpy.dot(matd_n[nlags+obso-j,:],betta[:,2+(j-1)*ncols:2+j*ncols].T)+betta[:,1].T                    
            matd_n[nlags+obso,:] += resmat_n[obso,:]
        nmatd_n = pymaclab.stats.common.genXX(self,matd=matd_n,func=True)
        betta_n = pymaclab.stats.common.compbetta(self,matd=matd_n,nmatd=nmatd_n,smbetta=False,func=True)
        betta_one_n = self.compbettax(matd=matd_n,nmatd=nmatd_n,smbetta=False,func=True)
        compmatr_n = self.mkCompMat(betta=betta_n,betta_one=betta_one_n,func=True)
        # Lets try here to actually keep the Cholesky shock matrix fixed at the original one...
        yfmat_n,resmat_n,covmat_n = self.mkfitres(matd=matd_n,nmatd=nmatd_n,betta=betta_n,func=True)
        vcmat_n,cholmat_n,eyem_n = self.mksigma(resmat=resmat_n,func=True)
        cholmat_n = copy.deepcopy(self.cholmat)
        phis_n,phisc_n = self.mkphis(compmatr=compmatr_n,betta=betta_n,maxphi=maxphi,cholmat=cholmat_n,func=True)
        return phisc_n

    def mkboot_pp(self,bdraw=None,yfmat=None,resiarray=None):
        # Set instance into boot mode
        self.boot = True
        # Get all of the paramters
        if bdraw == None: bdraw = self.confdic['bdraw']
        maxphi = self.confdic['maxphi']
        signif = self.confdic['signif']
        irfp = self.confdic['irfp']
        const_term = copy.deepcopy(self.confdic['const_term'])
        betta = copy.deepcopy(self.betta)
        resmat = copy.deepcopy(self.resmat)
        ncols = copy.deepcopy(self.ncols)
        nlags = copy.deepcopy(self.nlags)
        # Need to centre the residuals prior to using them in boostrapping
        resmat = resmat - numpy.mean(resmat,axis=0)
        if yfmat == None: yfmat = copy.deepcopy(self.yfmat)
        matd = copy.deepcopy(self.data)
        # Produce bootstrap confidence intervals, here using parallel python to speed up computation
        barray = []
        ppservers = ()
        if self.confdic['ncpus'] == 'auto':
            job_server = pp.Server(ncpus='autodetect',ppservers=ppservers)
        else:
            job_server = pp.Server(ncpus=self.confdic['ncpus'],ppservers=ppservers)
        # Do the first bootstrap run a la Killian 1998 to reduce the small sample bias of the IRFs, if wanted
        if self.confdic['killian_bsinbs']:
            self.phis_uncorr = copy.deepcopy(self.phis)
            print "Calculating bias correction using Killian's(1998) bs-in-bs method"
            barray_biasred = []
            # Note: Only pass the betta coefficients of the dynamic responses, not constant or time-trend
            inputs_biasred = ((resmat,yfmat,matd,nlags,maxphi,betta),)*bdraw
            barray_biasred_jobs = [job_server.submit(self.refac_killian_pp,(input,),modules=("copy","numpy","scipy","pymaclab.stats.common",)) for input in inputs_biasred]
            barray_biasred = [x() for x in barray_biasred_jobs]
            barray_biasred = numpy.array(barray_biasred)
            biasred_term = numpy.mean(barray_biasred-self.betta,axis=0)
            # Make the bias in value and percentage available for inspection
            self.killbias = copy.deepcopy(biasred_term)
            self.killbiasperc = copy.deepcopy((biasred_term/self.betta)*100.0)
            self.betta_kcorr = copy.deepcopy(self.betta + biasred_term)
            betta_corrected = self.betta + biasred_term
            # Pass the corrected betta on to the boostrap procedure
            betta = copy.deepcopy(betta_corrected)
            self.betta_uncorr = copy.deepcopy(self.betta)
            if 'betta_one' in dir(self):
                self.betta_one_uncorr = copy.deepcopy(self.betta_one)
                betta_one_kcorr = copy.deepcopy(self.betta_one)
                if self.confdic['const_term'] == 'cc' or self.confdic['const_term'] == 'tt':
                    betta_one_kcorr[1:,:] = copy.deepcopy(self.betta_kcorr)
                elif self.confdic['const_term'] == 'ct':
                    betta_one_kcorr[2:,:] = copy.deepcopy(self.betta_kcorr)
                elif self.confdic['const_term'] == 'ctt':
                    betta_one_kcorr[3:,:] = copy.deepcopy(self.betta_kcorr)
                compmatr = self.mkCompMat(betta=self.betta_kcorr,betta_one=betta_one_kcorr,func=True)
            else:
                compmatr = self.mkCompMat(betta=self.betta_kcorr,func=True)
            self.phisc_uncorr = copy.deepcopy(self.phisc)
            # Recalculation with new betta
            #yfmat,resmat,covmat = self.mkfitres(betta=self.betta_kcorr,func=True)
            #self.covmat_kcorr = copy.deepcopy(covmat)
            #self.covmat = copy.deepcopy(self.covmat_kcorr)
            #self.resmat_kcorr = copy.deepcopy(resmat)
            #self.resmat = copy.deepcopy(self.resmat_kcorr)
            #self.yfmat_kcorr = copy.deepcopy(yfmat)
            #self.yfmat = copy.deepcopy(self.yfmat_kcorr)
            # Needs to be called before new impulses are created, because colmat needs to be updated
            #self.cholmat_kcorr = self.mksigma(resmat=self.resmat_kcorr,func=True)[1]
            #self.cholmat_uncorr = copy.deepcopy(self.cholmat)
            #self.cholmat = copy.deepcopy(self.cholmat_kcorr)
            phis_kcorr,phisc_kcorr = self.mkphis(betta=self.betta_kcorr,compmatr=compmatr,func=True)
            self.phisc_kcorr = copy.deepcopy(phisc_kcorr)
            self.phis_kcorr = copy.deepcopy(phis_kcorr)
            self.phisc_uncorr = copy.deepcopy(self.phisc)
            self.phisc = copy.deepcopy(self.phisc_kcorr)
            print 'Calculation of bias correction done...'
        print 'Now computing bootstraps'
        # Note: Again, pass only the betta coefficients for dynamic responses, no constant or time-trend
        inputs = ((resmat,yfmat,matd,nlags,maxphi,betta),)*bdraw
        barray_jobs = [job_server.submit(self.refac_pp,(input,),modules=("copy","numpy","scipy","pymaclab.stats.common")) for input in inputs]
        barray = [x() for x in barray_jobs]
        print 'Calculation of bootstraps done...'
        barray = numpy.array(barray)
        if type(signif) == type(10.0):
            lower = scipy.stats.scoreatpercentile(barray,per=((100.0-signif*100.0)/2.0))
            upper = scipy.stats.scoreatpercentile(barray,per=signif*100.0+((100.0-signif*100.0)/2.0))
            confint = [lower,upper]
        elif type(signif) == type([10.0,10.0]):
            lower_outer = scipy.stats.scoreatpercentile(barray,per=(100.0-signif[1]*100.0)/2.0)
            lower_inner = scipy.stats.scoreatpercentile(barray,per=(100.0-signif[0]*100.0)/2.0)
            upper_inner = scipy.stats.scoreatpercentile(barray,per=signif[0]*100.0+((100.0-signif[0]*100.0)/2.0))
            upper_outer = scipy.stats.scoreatpercentile(barray,per=signif[1]*100.0+((100.0-signif[1]*100.0)/2.0))
            confint = [lower_outer,lower_inner,upper_inner,upper_outer]
        self.barray = barray
        self.confint = numpy.array(confint)[:,:irfp+1]
        # Switch out of boot mode
        self.boot = False
        
    def genFEVD(self,y,coefs,steps,func=False):
        p = len(coefs)
        k = len(coefs[0])
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
            
        self.forcs = copy.deepcopy(forcs)

        if 'phis_kcorr' in dir(self): ma_coefs = copy.deepcopy(self.phis_kcorr[:steps])
        else: ma_coefs = copy.deepcopy(self.phis[:steps])
        sigma_u = copy.deepcopy(self.covmat)
        k = len(sigma_u)
        forc_covs = numpy.zeros((steps, k, k))

        prior = numpy.zeros((k, k))
        for h in xrange(steps):
            # Sigma(h) = Sigma(h-1) + Phi Sig_u Phi'
            phi = ma_coefs[h]
            var = numpy.dot(numpy.dot(phi,sigma_u),phi.T)
            forc_covs[h] = prior = prior + var
            
        self.forc_covs = copy.deepcopy(forc_covs)
            
        periods = steps
        neqs = self.data.shape[1]
        vnames = self.vnames
        
        if 'phisc_kcorr' in dir(self): orthirfs = copy.deepcopy(self.phisc_kcorr)
        else: orthirfs = copy.deepcopy(self.phisc)

        # cumulative impulse responses
        irfs = (orthirfs[:periods] ** 2).cumsum(axis=0)
        self.forc_covs_irfs = copy.deepcopy(irfs)

        rng = range(neqs)
        mse = forc_covs[:, rng, rng]

        # lag x equation x component
        fevd = numpy.empty_like(irfs)

        for i in range(periods):
            fevd[i] = (irfs[i].T / mse[i]).T
        
        # Express in percentage terms
        fevd = fevd*100.0

        # switch to equation x lag x component
        self.decomp = fevd.swapaxes(0, 1)


    def plot_vals(self):
        # Graph some vars and dump them in graphs directory
        vnames = self.vnames
        matd = self.data
        nlags = self.nlags
        modname = self.confdic['modname']
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        yfmat = self.yfmat
        resmat = self.resmat
        tcode_dic = self.confdic['tcode']
        for col,name in enumerate(vnames):
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(3,1,1)
            ax1.grid()
            if tcode_dic[name] in [4,5,6,7,8,16]:
                ax1.plot(matd[nlags:,col]*100.0,color='black')
            else:
                ax1.plot(matd[nlags:,col],color='black')
            ax1.set_title('Actual Values for '+name)
            ax2 = fig1.add_subplot(3,1,2)
            ax2.grid()
            if tcode_dic[name] in [4,5,6,7,8,16]:
                ax2.plot(yfmat[:,col]*100.0,color='black')
            else:
                ax2.plot(yfmat[:,col],color='black')
            ax2.set_title('Fitted Values for '+name)
            ax3 = fig1.add_subplot(3,1,3)
            ax3.grid()
            if tcode_dic[name] in [4,5,6,7,8,16]:
                ax3.plot(resmat[:,col]*100.0,color='black')
            else:
                ax3.plot(resmat[:,col],color='black')
            ax3.set_title('Residuals for '+name)
            if self.confdic['graph_options']['save']['format']['eps']:
                fig1.savefig('../graphs/'+modname+'/'+name+'.eps',bbox_inches='tight')
            if self.confdic['graph_options']['save']['format']['pdf']:
                fig1.savefig('../graphs/'+modname+'/'+name+'.pdf',bbox_inches='tight')
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
        irfp = self.confdic['irfp']
        cols = self.ncols
        phisc = self.phisc
        confint = self.confint
        loglist = self.confdic['loglist']
        tcode_dic = self.confdic['tcode']
        translog = self.confdic['translog']
        transperc = self.confdic['transperc']
        transperc_dic = self.confdic['transperc_dic']
        transdic = self.trans_dic
        hzline = self.confdic['hzline']
        hzlwidth = self.confdic['hzline_width']
        bootstrap = self.confdic['bootstrap']
        kbsinbs = self.confdic['killian_bsinbs']
        if kbsinbs: phisc_uncorr = self.phisc_uncorr
        graph_biased = self.confdic['graph_options']['irfs']['biased_line']
        graph_biased_type = self.confdic['graph_options']['irfs']['biased_line_type']
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'shock_irfs' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/shock_irfs') 
        signif = self.confdic['signif']
        for shockv,sname in enumerate(vnames):
            fig1 = plt.figure()
            for respv,rname in enumerate(vnames):
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(1,1,1)
                if glopt['xlabel']: ax2.set_xlabel('Time periods after shock',fontsize=glopt['xlabel_fs'])
                ax = fig1.add_subplot(cols,1,respv+1)
                if glopt['xlabel']: ax.set_xlabel('Time periods after shock',fontsize=glopt['xlabel_fs'])
                if tcode_dic[rname] in [4,5,6,7,8,16] and translog and not transperc_dic[rname]:
                    if glopt['ylabel']: ax.set_ylabel(r'Percentage Deviation',fontsize=glopt['ylabel_fs'])
                    if glopt['ylabel']: ax2.set_ylabel(r'Percentage Deviation',fontsize=glopt['ylabel_fs'])
                else:
                    if glopt['ylabel']: ax.set_ylabel(r'Absolute Deviation',fontsize=glopt['ylabel_fs'])
                    if glopt['ylabel']: ax2.set_ylabel(r'Absolute Deviation',fontsize=glopt['ylabel_fs'])
                if hzline:
                    ax2.hlines(numpy.array([0,]*irfp),xmin=0,xmax=irfp,linestyles='solid',linewidth=hzlwidth,colors='k')
                    ax.hlines(numpy.array([0,]*irfp),xmin=0,xmax=irfp,linestyles='solid',linewidth=hzlwidth,colors='k')
                if tcode_dic[rname] in [4,5,6,7,8,16] and translog and not transperc_dic[rname]:
                    percen = 100.0
                elif transperc and transperc_dic[rname]:
                    percen = numpy.exp(transdic[vnames[respv]]['mean'])
                else:
                    percen = 1.0
                # When using the non-parametric bootstrap the central (mean) IR has to be the average of all IRs
                if not bootstrap:
                    ax2.plot(phisc[:irfp+1,respv,shockv].flatten()*percen,color=line_colour)
                elif bootstrap:
                    ax2.plot(phisc[:irfp+1,respv,shockv].flatten()*percen,color=line_colour)
                if kbsinbs and graph_biased:
                    ax2.plot(phisc_uncorr[:irfp+1,respv,shockv].flatten()*percen,graph_biased_type)
                plt.xlim((0,irfp))
                if gridopt: ax2.grid()
                ax.plot(phisc[:irfp+1,respv,shockv].flatten()*percen,color=line_colour)
                if gridopt: ax.grid()
                if bootstrap:
                    if type(signif) == type(10.0):
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
                        ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
                    elif type(signif) == type([10.0,10.0]):
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[1,:,respv,shockv].flatten()*percen,confint[2,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
                        ax2.fill_between([xx for xx in range(0,irfp+1,1)],confint[2,:,respv,shockv].flatten()*percen,confint[3,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
                        ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
                        ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[1,:,respv,shockv].flatten()*percen,confint[2,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
                        ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[2,:,respv,shockv].flatten()*percen,confint[3,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
                if glopt['title']: ax2.set_title('Response of '+pnames[rname]+' to shock on '+pnames[sname],fontsize=glopt['title_fs'])        
                if glopt['title']: ax.set_title('Response of '+pnames[rname]+' to shock on '+pnames[sname],fontsize=glopt['title_fs'])
                if self.confdic['graph_options']['save']['format']['eps']:
                    fig2.savefig('../graphs/'+modname+'/'+'shock_irfs/irf_'+sname+'_'+rname+'.eps',bbox_inches='tight')
                if self.confdic['graph_options']['save']['format']['pdf']:
                    fig2.savefig('../graphs/'+modname+'/'+'shock_irfs/irf_'+sname+'_'+rname+'.pdf',bbox_inches='tight')
                plt.close()
            if self.confdic['graph_options']['save']['format']['eps']:
                fig1.savefig('../graphs/'+modname+'/'+'irf_'+sname+'.eps',bbox_inches='tight')
            if self.confdic['graph_options']['save']['format']['pdf']:
                fig1.savefig('../graphs/'+modname+'/'+'irf_'+sname+'.pdf',bbox_inches='tight')
            plt.close()
        plt.close('all')


    def plot_one_irf(self,sname=None,rname=None,save_hd=True):       
        irf_options = self.confdic['graph_options']['irfs']
        outer_colour = irf_options['outer_colour']
        inner_colour = irf_options['inner_colour']
        line_colour = irf_options['line_colour']
        vnames = self.vnames
        shockv = vnames.index(sname)
        respv = vnames.index(rname)
        pnames = self.pnames
        irfp = self.confdic['irfp']
        cols = self.ncols
        phisc = self.phisc
        confint = self.confint
        loglist = self.confdic['loglist']
        tcode_dic = self.confdic['tcode']
        translog = self.confdic['translog']
        transperc = self.confdic['transperc']
        transperc_dic = self.confdic['transperc_dic']
        transdic = self.trans_dic
        hzline = self.confdic['hzline']
        bootstrap = self.confdic['bootstrap']
        modname = self.confdic['modname']
        # Check if folder(s) needs to be created
        if modname not in os.listdir('../graphs/'):
            os.mkdir('../graphs/'+modname)
        if 'shock_irfs' not in os.listdir('../graphs/'+modname+'/'):
            os.mkdir('../graphs/'+modname+'/shock_irfs') 
        signif = self.confdic['signif']
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('Lags')
        plt.xlim((0,irfp))
        ax.set_xlabel('Lags')
        plt.xlim((0,irfp))
        if tcode_dic[rname] in [4,5,6,7,8,16] and translog and not transperc_dic[rname]:
            ax.set_ylabel('\% Deviation')
        else:
            ax.set_ylabel('Abs. Deviation')
        if hzline:
            ax.hlines(numpy.array([0,]*irfp),xmin=0,xmax=irfp,linestyles='solid',linewidth=0.2,colors='k')
        if tcode_dic[rname] in [4,5,6,7,8,16] and translog and not transperc_dic[rname]:
            percen = 100.0
        elif transperc and transperc_dic[rname]:
            percen = numpy.exp(transdic[vnames[respv]]['mean'])
        else:
            percen = 1.0
            ax.plot(phisc[:irfp+1,respv,shockv].flatten()*percen,color=line_colour)
        ax.grid()
        if bootstrap:
            if type(signif) == type(10.0):
                ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
            elif type(signif) == type([10.0,10.0]):
                ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[0,:,respv,shockv].flatten()*percen,confint[1,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
                ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[1,:,respv,shockv].flatten()*percen,confint[2,:,respv,shockv].flatten()*percen, facecolor=outer_colour, edgecolor=outer_colour)
                ax.fill_between([xx for xx in range(0,irfp+1,1)],confint[2,:,respv,shockv].flatten()*percen,confint[3,:,respv,shockv].flatten()*percen, facecolor=inner_colour, edgecolor=inner_colour)
        ax.set_title(r'$\textrm{Response of '+pnames[rname]+' to shock on '+pnames[sname]+'}$')        
        if save_hd:
            if self.confdic['graph_options']['save']['format']['eps']:
                fig.savefig('../graphs/'+modname+'/'+'shock_irfs/irf_'+sname+'_'+rname+'.eps',bbox_inches='tight')
            if self.confdic['graph_options']['save']['format']['pdf']:
                fig.savefig('../graphs/'+modname+'/'+'shock_irfs/irf_'+sname+'_'+rname+'.pdf',bbox_inches='tight')
        else:
            self.irf_graphs_dic[sname][rname] = fig     
        plt.close('all')