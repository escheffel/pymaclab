import numpy
import copy
from ..filters import hpfilter, bkfilter, cffilter


def stdX(self,data=None,standard=False,only_demean=False,func=False):
    '''
    Function used to standardize the original data, i.e. de-mean and divide by SD,
    also saves data's characteristics upon normalizing.
    
    :param data:             The data to be normalizes/transformed
    :type data:              arr2d
    :param standard:         Should data be normalized by standard deviation?
    :type standard:          bool
    :param only_demean:      Should data only be de-meaned but not divided by std?
    :type only_demean:       bool
    :param func:             If true then method is behaving like a normal functions which returns
    :type func:              bool
    
    :return self.tdata:      *(arr2d)* - The transformed data
    :return self.nstdata:    *(arr2d)* - The nonstandardized or raw data which was passed to the method
    :return self.trans_dic:  *(dic)* - A dictionary containing the mean and std of the data. Uses vnames.
    '''
    if data != None:
        data = copy.deepcopy(data)
    else:
        data = copy.deepcopy(self.data)
    vnames = self.vnames
    trans_dic = {}
    data_nonstd = copy.deepcopy(data)
    for colo in range(0,data.shape[1],1):
        trans_dic[vnames[colo]] = {}
        mean_i = numpy.average(data[:,colo])
        trans_dic[vnames[colo]]['mean'] = copy.deepcopy(mean_i)
        if standard:
            data[:,colo] = data[:,colo] - mean_i
        std_i = numpy.std(data[:,colo])
        trans_dic[vnames[colo]]['std'] = copy.deepcopy(std_i)
        data_nonstd[:,colo] = copy.deepcopy(data[:,colo])
        if standard and not only_demean: data[:,colo] = data[:,colo]/std_i
    if not func:
        self.tdata = copy.deepcopy(data)
        self.nstddata = copy.deepcopy(data_nonstd)
        self.trans_dic = copy.deepcopy(trans_dic)
        return self
    elif func:
        return data,trans_dic
    
def transX(self,data=None,func=False):
    '''
    Function to transform the raw data according to some classification
    1 = levels
    2 = first seasonal difference
    3 = second seasonal difference
    4 = log level
    5 = log first seasonal difference
    6 = log second seasonal difference
    7 = detrend log using hp filter monthly data
    8 = detrend log using hp filter quarterly data
    16 = log second seasonal difference
    17 = (1-L)(1-L^12)
    18 = log of (1-L)*annualizing factor (i.e. x4 for quarterly and x12 for monthly
    19 = linear detrend of level (raw) data
    20 = linear detrend of log data
    21 = quadratic detrend of level (raw) data
    22 = quadratic detrend of log data
    '''
    if data != None:
        data = copy.deepcopy(data)
    elif data == None:
        data = copy.deepcopy(self.data)
    vnames = self.vnames
    freq = self.confdic['freq']
    tcode = self.confdic['tcode']
    shift = 0
    freqshift = 1
    if freq == 'M':
        freqshift = 12
    elif freq == 'Q':
        freqshift = 4
    else:
        freqshift = 1
    # for level do nothing
    # first seasonal difference
    for i1,vname in enumerate(vnames):
        if tcode[vname] == 2:
            data[1*freqshift:,i1] = data[1*freqshift:,i1]-data[:-1*freqshift,i1]
            if shift < 1: shift = 1*freqshift
        elif tcode[vname] == 3:
            data[1*freqshift:,i1] = data[1*freqshift:,i1]-data[:-1*freqshift,i1]
            data[2*freqshift:,i1] = data[2*freqshift:,i1]-data[:-2*freqshift,i1]
            if shift < 2: shift = 2*freqshift
        elif tcode[vname] == 4:
            data[:,i1] = numpy.log(data[:,i1])
        elif tcode[vname] == 5:
            data[:,i1] = numpy.log(data[:,i1])
            data[1*freqshift:,i1] = data[1*freqshift:,i1]-data[:-1*freqshift,i1]
            if shift < 1: shift = 1*freqshift
        elif tcode[vname] == 6:
            data[:,i1] = numpy.log(data[:,i1])
            data[1*freqshift:,i1] = data[1*freqshift:,i1]-data[:-1*freqshift,i1]
            data[2*freqshift:,i1] = data[2*freqshift:,i1]-data[:-2*freqshift,i1]
            if shift < 2: shift = 2*freqshift
        elif tcode[vname] == 7:
            data[:,i1] = numpy.log(data[:,i1])
            data[:,i1] = hpfilter(data=data[:,i1],lam=129600)[0]
        elif tcode[vname] == 8:
            data[:,i1] = numpy.log(data[:,i1])
            data[:,i1] = hpfilter(data=data[:,i1],lam=1600)[0]
        elif tcode[vname] == 18:
            data[:,i1] = numpy.log(data[:,i1])
            if freq == 'Q':
                data[1:,i1] = (data[1:,i1]-data[:-1,i1])*4.0
            elif freq == 'M':
                data[1:,i1] = (data[1:,i1]-data[:-1,i1])*12.0
            if shift < 1: shift = 1
        elif tcode[vname] == 19:
            data[:,i1] = scipy.signal.filter(data[:,i1], type='linear')
        elif tcode[vname] == 20:
            data[:,i1] = scipy.signal.filter(numpy.log(data[:,i1]), type='linear')
        elif tcode[vname] == 21:
            data[:,i1] = scipy.signal.bsplines.quadratic(data[:,i1])
        elif tcode[vname] == 22:
            data[:,i1] = scipy.signal.bsplines.quadratic(numpy.log(data[:,i1]))
    if shift > 0: data = data[shift:,:]
    if not func:
        self.max_shift = shift
        self.data = copy.deepcopy(data)
        return self
    else:
        return data