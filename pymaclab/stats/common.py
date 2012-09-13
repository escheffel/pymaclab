import numpy
import copy

def genXX(self,matd=None,nlags=None,const='standard',func=False):
    # Get all the needed variables from instance
    if nlags == None: nlags = copy.deepcopy(self.nlags)
    if matd == None: matd = copy.deepcopy(self.data)
    if const == 'standard': const_term = copy.deepcopy(self.confdic['const_term'])
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
    if const_term == 'tt':
        constar = numpy.reshape(numpy.array([x for x in range(0,nmatd.shape[0],1)]),(nmatd.shape[0],1))
        nmatd = numpy.hstack((constar,nmatd))
    elif const_term == 'cc':
        constar = numpy.reshape(numpy.array([1,]*nmatd.shape[0]),(nmatd.shape[0],1))
        nmatd = numpy.hstack((constar,nmatd))
    elif const_term == 'ct':
        constar = numpy.reshape(numpy.array([1,]*nmatd.shape[0]),(nmatd.shape[0],1))
        consttrend = numpy.reshape(numpy.array([x for x in range(0,nmatd.shape[0],1)]),(nmatd.shape[0],1))
        nmatd = numpy.hstack((constar,consttrend,nmatd))
    elif const_term == 'ctt':
        constar = numpy.reshape(numpy.array([1,]*nmatd.shape[0]),(nmatd.shape[0],1))
        consttrend = numpy.reshape(numpy.array([x for x in range(0,nmatd.shape[0],1)]),(nmatd.shape[0],1))
        consttrend2 = consttrend**2
        nmatd = numpy.hstack((constar,consttrend,consttrend2,nmatd))
    if not func:
        self.nmatd = nmatd
        return self
    elif func:
        return nmatd