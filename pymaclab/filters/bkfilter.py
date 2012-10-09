from pymaclab.filters._bkfilter import bkfilter as bkf
import numpy as np

def bkfilter(data=None,up=6,dn=32,kkl=12):
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
    return bkf(data,up,dn,kkl,tlen)
