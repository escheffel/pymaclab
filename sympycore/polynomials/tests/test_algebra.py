
from sympycore import *

P = PolynomialRing

def test_default_ring():
    r = repr(P.convert(0))
    assert r=="PolynomialRing('0')", r
    assert P.zero==P.convert(0)
    assert str(0)=='0'
    assert `P.convert(2)`=="PolynomialRing('2')",repr(P.convert(2))

def test_X():
    X = PolynomialRing['x']
    C = X.convert
    assert `X.zero`=="PolynomialRing[('x',), Calculus]('0')",`X.zero`
    assert `X.one`=="PolynomialRing[('x',), Calculus]('1')", repr(X.one)
    assert str(X.zero)=='0',`str(X.zero)`
    assert str(X.one)=='1',`str(X.zero)`

    assert `C({2:3})`=="PolynomialRing[('x',), Calculus]('3*x**2')", `C({2:3})`
    assert `C({2:3,4:5})`=="PolynomialRing[('x',), Calculus]('5*x**4 + 3*x**2')",`C({2:3,4:5})`
    assert `C([(2,3)])`=="PolynomialRing[('x',), Calculus]('3*x**2')"
    assert `C([(2,3),(4,5)])`=="PolynomialRing[('x',), Calculus]('5*x**4 + 3*x**2')"

    assert str(C({2:3}))=="3*x**2"
    assert str(C([(2,3),(4,5)]))=="5*x**4 + 3*x**2"

    assert `C('x*3*x')`=="PolynomialRing[('x',), Calculus]('3*x**2')",`C('x*3*x')`
    assert `C('x**4+4*x**4-3*x*x+6*x**2')`=="PolynomialRing[('x',), Calculus]('5*x**4 + 3*x**2')"

    assert `C([1,2,3])`=="PolynomialRing[('x',), Calculus]('3*x**2 + 2*x + 1')",repr(C([1,2,3]))

def test_univariate():
    X = PolynomialRing['x']
    C = X.convert
    x = C([0,1])
    assert `x`=="PolynomialRing[('x',), Calculus]('x')",`x`
    assert C([1, 2, 3]) == 3*x**2 + 2*x + 1,`C([1, 2, 3]).pair, (3*x**2 + 2*x + 1).pair`
    assert C([1, 2, 3]).degree == 2
    assert C([1, 2, 3]).ldegree == 0
    assert C([0]).degree == 0
    assert C([0]).ldegree == 0
    assert (x**3 + x) + (3*x + x**4) == x**4 + x**3 + 4*x
    assert C([1, 2, 3])*2 == C([2, 4, 6])
    assert C([1, 2, 3])/2.0 == C([0.5, 1.0, 1.5])
    assert C([4,-3,5])(6) == 166
    assert C([1,2,3])(C([4,5,6])) == 108*x**4 + 180*x**3 + 231*x**2 + 130*x + 57
    assert (x-1)*(x+1) == x**2 - 1
    assert (x-1)*(x+1) / (x-1) == (x+1)
    assert (x**3 + x).diff() == 3*x**2 + 1
    assert str(C([3, 0, 1])) == 'x**2 + 3',str(C([3, 0, 1]))
    assert str(C([0])) == '0'
    assert str(C([1, 1])) == 'x + 1'
    assert str((5 + 3*x) / 5) == '3/5*x + 1', repr((5 + 3*x) / 5)
    assert str(C([4, 11, 6]) / C([6, 12])) == '1/2*x + 2/3'
    assert str(C([2, 3, 4]) % C([1, 2, 3])) == '1/3*x + 2/3'
