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

    obj = Matrix([[1,x],[y,2]])
    assert obj==pickler(obj)
