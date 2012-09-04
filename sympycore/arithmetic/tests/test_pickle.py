import os
import tempfile
import pickle

from sympycore.arithmetic import mpq, mpqc, Infinity, mpf, mpc

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

    obj = mpq((1,2))
    assert obj==pickler(obj)

    obj = mpqc(1,2)
    assert obj==pickler(obj)

    obj = Infinity(-1)
    assert obj==pickler(obj)

    obj = mpf('0.5')
    assert obj==pickler(obj)

    obj = mpc('0.5', '0.2')
    assert obj==pickler(obj)
