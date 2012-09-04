
import os
import tempfile
import pickle

from sympycore.core import Expr
from sympycore.utils import SYMBOL

class MyExpr(Expr):
    # see comment in sympycore/core.py for pickling pure Expr
    # instances.
    pass

def test_pickle():

    fn = tempfile.mktemp()
    obj = MyExpr(SYMBOL, 'x')
    
    f = open(fn, 'wb')
    pickle.dump(obj, f)
    f.close()

    f = open(fn, 'rb')
    obj2 = pickle.load(f)
    f.close()
    os.remove(fn)
    
    assert obj.head == obj2.head,`obj.head, obj2.head`
    assert obj.head is obj2.head,`obj.head, obj2.head, id(obj.head), id(obj2.head)`
    assert obj.data == obj2.data,`obj.data, obj2.data`
    assert obj.pair == obj2.pair,`obj.pair, obj2.pair`
    assert obj==obj2,`obj,obj2`

