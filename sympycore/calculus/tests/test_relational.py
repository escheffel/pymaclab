
from sympycore import *
from sympycore.calculus.relational import *
x = Symbol('x')

A = Assumptions

def test_inequalities():

    assert A().check(2 < 3) is True
    assert A().check(Rational(2,3) < Rational(7,2)) is True
    assert A().check(pi**2 > 0) is True
    assert A().check(pi**2 < 0) is False
    assert A().check(3+4*pi > 0) is True
    assert A().check(3-pi > 0) is None
    assert A().check(-4*pi < 0) is True
    assert A().check(-4*pi > 0) is False
    assert A().check(-4*pi**Rational(1,2) < 0) is True
    assert A().check(-4*pi**Rational(1,2) > 0) is False
    assert A().check(-4*pi**2 - 1 < 0) is True
    assert A([x > 3]).check(x > 2) is True
    assert A([x > 3]).check(x > 3) is True
    assert A([x > 3]).check(x > 4) is None
    #assert A([x > 3]).check(x < 3) is False
    assert A([x > 3]).check(x > 2) is True

def XXXtest_double_inequalities():
    assert A([x>0, x<1]).check(x < 1) is True
    assert A([x>0, x<1]).check(x > 0) is True
    assert A([x>0, x<1]).check(x > 2) is False
    assert A([x>0, x<1]).check(x < -1) is False
