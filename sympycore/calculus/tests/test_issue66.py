
from sympycore import *
def test_reproduce_issue():
    e = Calculus('a*b')
    assert e.subs(dict (a=2,b=2))==4
    r = e.subs(dict (a=2,b=2.0))
    assert r==Calculus(4.0),`r.pair`
