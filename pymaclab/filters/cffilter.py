from pymaclab.filters import _cffilter as cff
import numpy as np

def cffilter(data=None,low=6,high=32,drift=True):  
    return cff(data=data,low=low,high=high,drift=drift) 
