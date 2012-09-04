
from sympycore import *

def test_number():
    true = Logic.true
    false = Logic.false
    assert str(true)=='True'
    assert str(false)=='False'
    assert not false
    assert true
    assert true==True
    assert false==False

    assert repr(true)=="Logic('True')"
    assert repr(false)=="Logic('False')"
    assert Logic(True)==true
    assert Logic(False)==false

    assert Logic('True')==true
    assert Logic('False')==false

    assert Logic('not True')==false
    assert Logic('not False')==true

    assert Logic('True or True')==true
    assert Logic('False or True')==true
    assert Logic('True or False')==true
    assert Logic('False or False')==false

    assert Logic('True and True')==true
    assert Logic('False and True')==false
    assert Logic('True and False')==false
    assert Logic('False and False')==false

def test_relational():
    x = Symbol('x')
    y = Symbol('y')

    assert (x<y) == Lt(x, y),`x<y, Lt(x,y)`
    assert (x<=y) == Le(x, y)
    assert (x>y) == Gt(x, y)
    assert (x>=y) == Ge(x, y)

    assert Logic('x<y') == Lt(x,y)
    assert Logic('x<=y') == Le(x,y)
    assert Logic('x>y') == Gt(x,y)
    assert Logic('x>=y') == Ge(x,y)
    assert Logic('x==y') == Eq(x,y)
    assert Logic('x!=y') == Ne(x,y)

    assert repr(Lt(x,y))=="Logic('x<y')"
    assert repr(Le(x,y))=="Logic('x<=y')"
    assert repr(Gt(x,y))=="Logic('x>y')"
    assert repr(Ge(x,y))=="Logic('x>=y')"
    assert repr(Eq(x,y))=="Logic('x==y')"
    assert repr(Ne(x,y))=="Logic('x!=y')"

    assert Not(Lt(x,y))==Ge(x,y)
    assert Not(Le(x,y))==Gt(x,y)
    assert Not(Gt(x,y))==Le(x,y)
    assert Not(Ge(x,y))==Lt(x,y)
    assert Not(Eq(x,y))==Ne(x,y)
    assert Not(Ne(x,y))==Eq(x,y)

def test_subs():
    true = Logic.true
    false = Logic.false
    assert Logic('x or y').subs('y', 'not x')==true
    assert Logic('x and y').subs('y', 'not x')==false
    assert Logic('not y').subs('y', 'not x')==Logic('x')

    assert Logic('x<=1 and x>=-2').subs('x', 2)==false
    assert Logic('x<=1 and x>=-2').subs('x', 1)==true
    assert Logic('x<=1 and x>=-2').subs('x', 0)==true
    assert Logic('x<=1 and x>=-2').subs('x', -2)==true
    assert Logic('x<=1 and x>=-2').subs('x', -3)==false

def test_in():
    assert Logic('x in y') == Logic.Element(classes.Verbatim('x'), classes.Set('y'))
    assert Logic('x in Integers') == Logic.Element(classes.Calculus('x'), classes.Set('Integers'))
    
