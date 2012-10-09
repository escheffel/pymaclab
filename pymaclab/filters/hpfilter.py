from pymaclab.filters._hpfilter import hpfilter as hpf
import numpy as np

def hpfilter(data=None,lam=1600):
    if type(data) == type(np.matlib.matrix([1,2,3])) and len(data.shape) == 2:
        if data.shape[0] < data.shape[1]: data = data.__array__().T
        else: data = data.__array__()
    elif type(data) != type(np.array([1,2,3])):
        data = np.array([x for x in data])
    elif type(data) == type(np.array([1,2,3])):
        if len(data.shape) == 2 and data.shape[0] < data.shape[1]:
            data = data.reshape(data.shape[1])
        elif len(data.shape) == 2 and data.shape[0] > data.shape[1]:
            data = data.reshape(data.shape[0])
    tlen = len(data)
    wdata = np.zeros((tlen,3))
    lam = 1600    
    return hpf(data,wdata,tlen,lam,0)