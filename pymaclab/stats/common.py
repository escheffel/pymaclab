import numpy
import copy
from ..linalg import pca_module

def genXX(self,matd=None,nlags=None,const='standard',func=False):
    # Get all the needed variables from instance
    if nlags == None: nlags = copy.deepcopy(self.nlags)
    if matd == None: matd = copy.deepcopy(self.data)
    if const == 'standard': const = copy.deepcopy(self.confdic['const_term'])
    rows = matd.shape[0]
    cols = matd.shape[1]
    nmatd = numpy.zeros([rows-nlags,cols*nlags])
    for lag in range(0,nlags,1):
        for col in range(0,cols,1):
            shift = lag*(cols-1)
            if lag == 0:
                nmatd[:,col+shift+lag] = matd[lag:-(nlags-(lag)),col]
            elif lag == nlags - 1:
                nmatd[:,col+shift+lag] = matd[lag:-1,col]
            elif lag > 0:
                nmatd[:,col+shift+lag] = matd[lag:-(nlags-(lag)),col]
    # Flip matrix left to right in order to have p(1) lags first, then re-order variables
    nmatd = numpy.fliplr(nmatd)
    for elem in range(0,nlags,1):
        nmatd[:,elem*cols:(elem+1)*cols] = numpy.fliplr(nmatd[:,elem*cols:(elem+1)*cols])
    # Add the colum for the constants
    if const == 'tt':
        constar = numpy.reshape(numpy.array([x for x in range(0,nmatd.shape[0],1)]),(nmatd.shape[0],1))
        nmatd = numpy.hstack((constar,nmatd))
    elif const == 'cc':
        constar = numpy.reshape(numpy.array([1,]*nmatd.shape[0]),(nmatd.shape[0],1))
        nmatd = numpy.hstack((constar,nmatd))
    elif const == 'ct':
        constar = numpy.reshape(numpy.array([1,]*nmatd.shape[0]),(nmatd.shape[0],1))
        consttrend = numpy.reshape(numpy.array([x for x in range(0,nmatd.shape[0],1)]),(nmatd.shape[0],1))
        nmatd = numpy.hstack((constar,consttrend,nmatd))
    elif const == 'ctt':
        constar = numpy.reshape(numpy.array([1,]*nmatd.shape[0]),(nmatd.shape[0],1))
        consttrend = numpy.reshape(numpy.array([x for x in range(0,nmatd.shape[0],1)]),(nmatd.shape[0],1))
        consttrend2 = consttrend**2
        nmatd = numpy.hstack((constar,consttrend,consttrend2,nmatd))
    if not func:
        self.nmatd = nmatd
        return self
    elif func:
        return nmatd
    

def is_stable(compmat=None,verbose=False):
    """
Determine stability of VAR(p) system by examining the eigenvalues of the
VAR(1) representation

Parameters
----------
compmat : ndarray (p x k,p x k)

Returns
-------
is_stable : bool
"""
    eigs = numpy.linalg.eigvals(compmat)

    if verbose:
        print 'Eigenvalues of VAR(1) rep'
        for val in numpy.abs(eigs):
            print val

    return (numpy.abs(eigs) <= 1).all()


def genPCs(self,data=None,sfacs=None,func=False):
    '''
    Using PCA to extract the factors from some array of data
    This utility function uses non-stacked vectors or data, s stands for simple
    '''
    if sfacs == None: npc = copy.deepcopy(self.sfacs)
    else: npc = sfacs
    cols = self.ncols
    nlags = self.nlags
    pcm = pca_module.PCA_svd(data,standardize=False)[0][:,:npc]
    exp_var = pca_module.PCA_svd(data,standardize=False)[2][:npc]
    pcmm = numpy.zeros((len(pcm),len(pcm[0])))
    for factor in range(0,len(pcm),1):
        pcmm[factor,:] = pcm[factor]
    if not func:
        self.fdata = copy.deepcopy(pcmm)
        self.exp_var = copy.deepcopy(exp_var)
    elif func:
        return pcmm,exp_var



def genPC(self,data=None,func=False):
    '''
    Using PCA to extract the factors from some array of data
    This utility function uses the stacked data instead
    '''
    npc = copy.deepcopy(self.sfacs)
    cols = self.ncols
    flags = self.flags
    nlags = self.nlags
    stackX = copy.deepcopy(data)
    pcm = numpy.zeros((stackX.shape[0],npc*flags))
    for lago in range(0,flags,1):
        pcm[:,lago*npc:(lago+1)*npc] = pca_module.PCA_svd(stackX[:,lago*cols:(lago+1)*cols],standardize=False)[0][:,:npc]
    if not func:
        self.fdata = copy.deepcopy(pcm)
    elif func:
        return pcm

def compbetta(self,matd=None,nmatd=None,smbetta=True,const='standard',func=False):
    # Get all the needed variables from instance
    if matd == None: matd = self.data
    if nmatd == None: nmatd = self.nmatd
    cols = matd.shape[1]
    if const == 'standard': const = self.confdic['const_term']
    if const == None:
        nlags = int(nmatd.shape[1]/cols)
    elif const == 'cc':
        nlags = int((nmatd.shape[1]-1)/cols)
    elif const == 'tt':
        nlags = int((nmatd.shape[1]-1)/cols)
    elif const == 'ct':
        nlags = int((nmatd.shape[1]-2)/cols)
    elif const == 'ctt':
        nlags = int((nmatd.shape[1]-3)/cols)
    XX = numpy.dot(nmatd.T,nmatd)
    XXi = numpy.linalg.inv(XX)
    if const == None:
        betta = numpy.zeros((cols,nlags*cols))
    elif const == 'tt' or const == 'cc':
        betta = numpy.zeros((cols,1+nlags*cols))
    elif const == 'ct':
        betta = numpy.zeros((cols,2+nlags*cols))
    elif const == 'ctt':
        betta = numpy.zeros((cols,3+nlags*cols))
    for vari in range(0,cols,1):
        ymat=matd[nlags:,vari]
        Xy = numpy.dot(nmatd.T,ymat)
        betta[vari,:] = numpy.dot(XXi,Xy)
    if smbetta:
        # Produce the betta format from statsmodels library so we can use some of their code
        # Exclude the trend or constant from bettasm, as opposed to in betta
        bettasm = numpy.zeros((nlags,cols,cols))
        for lago in range(0,nlags,1):
            for varos in range(0,cols,1):
                if const == 'tt' or const == 'cc':
                    bettasm[lago][varos,:] = betta[:,1:][varos][cols*lago:cols*(lago+1)]
                elif const == None:
                    bettasm[lago][varos,:] = betta[varos][cols*lago:cols*(lago+1)]
                elif const == 'ct':
                    bettasm[lago][varos,:] = betta[:,2:][varos][cols*lago:cols*(lago+1)]
                elif const == 'ctt':
                    bettasm[lago][varos,:] = betta[:,3:][varos][cols*lago:cols*(lago+1)]
    if not func:
        if smbetta:
            self.bettasm = bettasm
            self.betta = betta
            return self
        elif not smbetta:
            self.betta = betta
            return self
    elif func:
        if smbetta:
            return bettasm
        elif not smbetta:
            return betta

def check_matrix(arraym=None,arrname=None,expdim=None):
    if type(arraym) == type(numpy.array([[1,2,3],[1,2,3]])):
        rows = arraym.shape[0]
        cols = arraym.shape[1]
    elif type(arraym) == type({}):
        cols = len(arraym[arraym.keys()[0]])
        rows = len(arraym.keys())
    print '%%%%%%%%%%%%%%% MATRIX DIM CHECK %%%%%%%%%%%%%%%%%%'
    print 'Matrix ##'+arrname+'## has dim ('+str(rows)+','+str(cols)+')'
    if expdim != None:
        print 'The expected dimension was: ('+str(expdim[0])+','+str(expdim[1])+')'
        if (rows,cols) == expdim: print 'Matrix dim check PASSED !'
    print '%%%%%%%%%%%%%%% DIM CHECK END %%%%%%%%%%%%%%%%%%%%%%'

  
def genStackX(self,data=None,flags=None,func=False):
    '''
    Utility function used to take data X (which may be n-dimensional on horizontal axis)
    and then to create a stacked vector based on flags 
    '''
    # Get all the needed variables from instance
    if flags == None:
        flags = copy.deepcopy(self.flags)
    matd = copy.deepcopy(data)
    rows = matd.shape[0]
    cols = matd.shape[1]
    vnames = copy.deepcopy(self.vnames)
    nmatd = numpy.zeros([rows-flags,cols*flags])
    for lag in range(0,flags,1):
        for col in range(0,cols,1):
            shift = lag*(cols-1)
            if lag == 0:
                nmatd[:,col+shift+lag] = matd[lag:-(flags-lag),col]
            elif lag == flags - 1:
                nmatd[:,col+shift+lag] = matd[lag:-1,col]
            elif lag > 0:
                nmatd[:,col+shift+lag] = matd[lag:-(flags-lag),col]
    # Flip matrix left to right in order to have p(1) lags first, then re-order variables
    nmatd = numpy.fliplr(nmatd)
    for elem in range(0,flags,1):
        nmatd[:,elem*cols:(elem+1)*cols] = numpy.fliplr(nmatd[:,elem*cols:(elem+1)*cols])
    # Add original data back to lags in first columns to get stacked X with current and lagged
    nmatd = numpy.hstack((matd[flags:,:],nmatd))
    if not func:
        self.stackmat = copy.deepcopy(nmatd)
        return self
    elif func:
        return nmatd


  
def genStackX2(self,data=None,flags=None,fleads=None,func=False):
    '''
    Utility function used to take data X (which may be n-dimensional on horizontal axis)
    and then to create a stacked vector based on flags WORK IN PROGRESS 
    '''
    # Get all the needed variables from instance
    if flags == None:
        flags = copy.deepcopy(self.flags)
    if fleads == None:
        fleads = 0
    matd = copy.deepcopy(data)
    rows = matd.shape[0]
    cols = matd.shape[1]
    vnames = copy.deepcopy(self.vnames)
    nmatd = numpy.zeros([rows-(flags+fleads),cols*flags])
    # Do this only if leads are also required !
    if fleads != None:
        nmatdf = numpy.zeros([rows-(flags+fleads),cols*fleads])
        for lead in range(0,fleads,1):
            for col in range(0,cols,1):
                shift = lead*(cols-1)
                if lead == 0:
                    nmatdf[:,col+shift+lead] = matd[flags+lead:-fleads,col]
                elif lead == fleads - 1:
                    nmatdf[:,col+shift+lead] = matd[flags+lead:-1,col]
                elif lead > 0:
                    nmatdf[:,col+shift+lead] = matd[flags+lead:-(fleads-lead),col]
    for lag in range(0,flags,1):
        for col in range(0,cols,1):
            shift = lag*(cols-1)
            if lag == 0:
                nmatd[:,col+shift+lag] = matd[lag:-(fleads+flags-lag),col]
            elif lag == flags - 1:
                nmatd[:,col+shift+lag] = matd[lag:-(1+fleads),col]
            elif lag > 0:
                nmatd[:,col+shift+lag] = matd[lag:-(fleads+flags-lag),col]
    # Flip matrix left to right in order to have p(1) lags first, then re-order variables
    nmatd = numpy.fliplr(nmatd)
    for elem in range(0,flags,1):
        nmatd[:,elem*cols:(elem+1)*cols] = numpy.fliplr(nmatd[:,elem*cols:(elem+1)*cols])
    # Add original data back to lags in first columns to get stacked X with current and lagged
    nmatd = numpy.hstack((matd[flags:-fleads,:],nmatd))
    # Stack also the leading terms if wanted !
    if fleads != None:
        nmatd = numpy.hstack((nmatd,nmatdf))
    # Also make a dictionary for each variable
    stackdic = {}
    for namo in vnames:
        tmpmat = numpy.zeros((rows-flags-fleads,flags+1+fleads))
        # Save the current-period value
        tmpmat[:,flags+1] = data[flags:fleads,vnames.index(namo)]
        for i1,lago in enumerate(xrange(flags,0,1)):
            tmpmat[:,i1] = nmatd[:,:]
            
    if not func:
        self.stackmat = copy.deepcopy(nmatd)
        return self
    elif func:
        return nmatd



def genEE(self,vnames=None,ydata=None,xdata=None,betta=None,hweightm=None,mlags=None,func=False):
    '''
    Function to generate equations residuals after controlling for own lags, used to extract factors
    Can also be used to estimate bettas, also optionally uses weighting matrix hweightm
    '''
    data = ydata
    cols = data.shape[1]
    nlags = self.nlags
    flags = self.flags
    npc = self.sfacs
    sfacs = npc
    if vnames == None:
        vnames = self.vnames
    auto_betta = {}
    auto_ee_dic = {}
    lxmatdic = xdata
    for varo in range(0,cols,1):
        Xmat = lxmatdic[vnames[varo]]
        Ymat = data[:,varo]
        row_min = min(Xmat.shape[0],Ymat.shape[0])
        Xmat_indi = Xmat.shape[0] - row_min
        Ymat_indi = Ymat.shape[0] - row_min
        Ymat = Ymat[Ymat_indi:]
        Xmat = Xmat[Xmat_indi:,:]
        if betta == None:
            if hweightm == None:
                XX = numpy.dot(Xmat.T,Xmat)
            else:
                XTW = numpy.dot(Xmat.T,hweightm)
                XX = numpy.dot(XTW,Xmat)
            XXi = numpy.linalg.inv(XX)
            if hweightm == None:
                Xy = numpy.dot(Xmat.T,Ymat)
            else:
                XTW = numpy.dot(Xmat.T,hweightm)
                Xy = numpy.dot(XTW,Ymat)
            auto_betta[vnames[varo]] = numpy.dot(XXi,Xy)
        if betta == None:
            auto_ee_dic[vnames[varo]] = Ymat - numpy.dot(Xmat,auto_betta[vnames[varo]])
        else:
            auto_ee_dic[vnames[varo]] = Ymat - numpy.dot(Xmat,betta[vnames[varo]])
            auto_betta = copy.deepcopy(betta)

    ee_mat = auto_ee_dic[vnames[0]]
    for varo in range(1,cols,1):
        ee_mat = numpy.vstack((ee_mat,auto_ee_dic[vnames[varo]]))
    ee_mat = ee_mat.T     
    if not func:
        self.auto_betta = copy.deepcopy(auto_betta)
        self.auto_ee_dic = copy.deepcopy(auto_ee_dic)
        self.ee_mat = copy.deepcopy(ee_mat)
        return self
    elif func:
        return auto_betta,auto_ee_dic,ee_mat


def genLX(self,vnames=None,ydata=None,nlags=None,mlags=None,func=False):
    '''
    Function to create a dictionary of the auto-regressors for each variable to eliminate serial corr
    so be careful, this will only create the regressor matrices involving OWN lags
    '''
    # Get all the needed variables from instance
    if nlags == None:
        nlags = copy.deepcopy(self.nlags)
    if mlags == None:
        mlags = nlags
    flags = copy.deepcopy(self.flags)
    if ydata == None:
        matd = copy.deepcopy(self.tdata)
    else:
        matd = ydata
    rows = matd.shape[0]
    cols = matd.shape[1]
    if vnames == None:
        vnames = copy.deepcopy(self.vnames)
    lxmatdic = {}
    for col in range(0,cols,1):
        nmatd = numpy.zeros([rows-(mlags),nlags])
        for lag in range(0,nlags,1):
            if lag == 0:
                nmatd[:,lag] = matd[lag:-((mlags)-lag),col]
            elif lag == nlags - 1:
                nmatd[:,lag] = matd[(mlags-1):-1,col]
            elif lag > 0:
                nmatd[:,lag] = matd[lag:-((mlags)-lag),col]
            # Flip matrix left to right in order to have p(1) lags first
        nmatd = numpy.fliplr(nmatd)
        lxmatdic[vnames[col]] = copy.deepcopy(nmatd)
    if not func:
        self.lxmatdic = copy.deepcopy(lxmatdic)
        return self
    elif func:
        return lxmatdic