from pymaclab.filters._cffilter import cffilter as cff
import numpy as np

def cffilter(data=None,low=6,high=32,drift=True):  
    return cff(data,low=low,high=high,drift=drift) 
