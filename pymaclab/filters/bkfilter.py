from pymaclab.filters._bkfilter import bkfilter as bkf
import numpy as np

def bkfilter(data=None,up=6,dn=32,kkl=12):
    if type(data) != type(np.array([1,2,3])):
        data = np.array(data)
    tlen = len(data)
    return bkf(data,up,dn,kkl,tlen)
