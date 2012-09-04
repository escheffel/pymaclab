from sympycore import *
def test_david_bug():
    a = Calculus('-x*(x-1)*r')
    assert a.subs('x',3)==Calculus('-6*r')

    assert a.subs('x',3.4).subs('r',3.4)==Calculus('-3.4*(3.4-1)*3.4')
