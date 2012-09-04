
from sympycore import *

A = Calculus
x,y,z = map(A.Symbol,'xyz')
f = A.convert(1.2)

def test_ops():
    assert f+1==A.convert(2.2)
    assert f+1.2==A.convert(2.4)

    assert `f+x` in ["Calculus('1.2 + x')", "Calculus('x + 1.2')"],`f+x`
    assert `1.2+x` in ["Calculus('1.2 + x')", "Calculus('x + 1.2')"],`1.2+x`
    #assert `f+pi`=="Calculus('4.34159265359')"
    #assert `1.2+pi`=="Calculus('4.34159265359')"

