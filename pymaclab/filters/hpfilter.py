from pymaclab.filters._hpfilter import hpfilter as hpf
import numpy as np

def hpfilter(data=None,lam=1600):
    if type(data) != type(np.array([1,2,3])):
        data = np.array(data)
    tlen = len(data)
    wdata = np.zeros((tlen,3))
    lam = 1600    
    return hpf(data,wdata,tlen,lam,0)