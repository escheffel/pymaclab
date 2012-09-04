import os
import tempfile
import pickle

from sympycore import *

def pickler(obj):
    fn = tempfile.mktemp()
    
    f = open(fn, 'wb')
    pickle.dump(obj, f)
    f.close()

    f = open(fn, 'rb')
    obj2 = pickle.load(f)
    f.close()
    os.remove(fn)

    return obj2

def test_pickle():

    x, y, z = map(Symbol,'xyz')

    obj = x+y
    assert obj==pickler(obj),`obj, pickler(obj)`

    obj = x+y/3
    assert obj==pickler(obj)

    obj = x+0.5*y
    assert obj==pickler(obj)

    obj = x+Sin(y)
    assert obj==pickler(obj)

    obj = x+Log(y)
    assert obj==pickler(obj)

