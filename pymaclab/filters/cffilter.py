from pymaclab.filters._cffilter import cffilter as cff
import numpy as np

def cffilter(data=None,low=6,high=32,drift=True):
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
    return cff(data,low=low,high=high,drift=drift) 
