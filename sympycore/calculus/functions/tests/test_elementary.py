from __future__ import with_statement

from sympycore.calculus import Calculus as A
from sympycore.calculus.algebra import I
from sympycore.calculus.infinity import oo, zoo, undefined
from sympycore.calculus import Number, Symbol
from sympycore.calculus.functions.elementary import Sin, Cos, Tan, Cot, pi, E, Exp, Log
from sympycore.calculus.relational import Assumptions
from sympycore import FD, D

def test_exp_log():
    assert Exp(1) == E
    assert Exp(0) == 1
    assert Log(0) == -oo
    assert Log(1) == 0
    assert Log(oo) == oo
    assert Log(2,3) == Log(2)/Log(3)
    assert Log(0,3) == (-oo)/Log(3)
    assert Log(100,10) == 2
    assert Log(65536,2) == 16
    assert Log(101,10) != 2
    assert Log(100,0) == 0
    assert Log(5,pi) == Log(5)/Log(pi)
    assert Log(-3) == I*pi + Log(3)
    assert Log(-pi) == I*pi + Log(pi)
    assert Log(I) == I*pi/2
    assert Log(-I) == -I*pi/2
    assert Log(3*I) == Log(3) + I*pi/2
    assert Log(-3*I) == Log(3) - I*pi/2
    assert str(Log(1+I)) == 'Log(1 + I)'
    assert Log(5**Number(1,3)) == Number(1,3)*Log(5), str(Log(5**Number(1,3)))
    assert Log(5**(3+2*I)) != (3+2*I)*Log(5)

def test_log_assumptions():
    x = Symbol('x')
    assert Log(x**2) != 2*Log(x)
    with Assumptions([x > 0]):
        assert Log(x**2) == 2*Log(x),str(Log(x**2))

def test_trig_values():
    sqrt2 = A.convert('2**(1/2)')
    sqrt3 = A.convert('3**(1/2)')
    assert Sin(0) == 0
    assert Sin(pi) == 0
    assert Sin(4*pi) == 0
    assert Sin(3*pi/2) == -1, `Sin(3*pi/2)`
    assert Sin(5*pi/2) == 1
    assert Sin(pi/3) == sqrt3/2, `Sin(pi/3).pair, (sqrt3/2).pair`
    assert Sin(pi/2) == 1
    assert Cos(0) == 1
    assert Cos(pi) == -1
    assert Cos(8*pi) == 1
    assert Cos(-9*pi) == -1
    assert Cos(pi/2) == 0
    assert Cos(3*pi/2) == 0
    assert Cos(11*pi/2) == 0
    assert Cos(pi/12) == (1 + sqrt3) / (2 * sqrt2)
    assert Tan(7*pi/12) == Sin(7*pi/12)/Cos(7*pi/12)
    assert Tan(pi/2) == zoo,`Tan(pi/2)`
    assert Tan(pi) == 0
    assert Cot(pi/2) == 0
    assert Cot(pi) == zoo
    assert str(Sin(oo)) == 'Sin(oo)'
    assert Sin(undefined) == undefined

def test_trig_symmetry():
    x = Symbol('x')
    assert Sin(-x) == -Sin(x),`Sin(-x).pair, (-Sin(x)).pair`
    assert Cos(-x) == Cos(x)
    assert Tan(-x) == -Tan(x)
    assert Cot(-x) == -Cot(x)
    assert Sin(x+pi) == -Sin(x),`Sin(x+pi)`
    assert Sin(x+2*pi) == Sin(x)
    assert Sin(x+3*pi) == -Sin(x)
    assert Sin(x+4*pi) == Sin(x)
    assert Sin(x-5*pi) == -Sin(x)
    assert Cos(x+pi) == -Cos(x)
    assert Cos(x+2*pi) == Cos(x)
    assert Cos(x+3*pi) == -Cos(x)
    assert Cos(x+4*pi) == Cos(x)
    assert Cos(x-5*pi) == -Cos(x)
    assert Tan(x+pi) == Tan(x)
    assert Tan(x-3*pi) == Tan(x)
    assert Cot(x+pi) == Cot(x)
    assert Cot(x-3*pi) == Cot(x)
    assert Sin(pi/2-x) == Cos(x)
    assert Sin(3*pi/2-x) == -Cos(x)
    assert Sin(5*pi/2-x) == Cos(x)
    assert Cos(pi/2-x) == Sin(x)
    assert Cos(3*pi/2-x) == -Sin(x)
    assert Cos(5*pi/2-x) == Sin(x)
    assert Tan(pi/2-x) == Cot(x)
    assert Tan(3*pi/2-x) == Cot(x)
    assert Tan(5*pi/2-x) == Cot(x)
    assert Cot(pi/2-x) == Tan(x)
    assert Cot(3*pi/2-x) == Tan(x)
    assert Cot(5*pi/2-x) == Tan(x)
    assert Sin(pi/2+x) == Cos(x)
    assert Cos(pi/2+x) == -Sin(x)
    assert Tan(pi/2+x) == -Cot(x)
    assert Cot(pi/2+x) == -Tan(x)

def test_trig_diff():
    x = Symbol('x')
    assert Sin(x).diff(x) == Cos(x), `Sin(x).diff(x)`
    assert Cos(x).diff(x) == -Sin(x)
    assert Sin(2*x).diff(x) == 2*Cos(2*x)

    assert Tan(x).diff(x) == 1+Tan(x)**2
    assert Cot(x).diff(x) == -1-Cot(x)**2

    assert Log(x).diff(x) == 1/x
    assert Exp(x).diff(x) == Exp(x)


    assert (x*Sin(x)).diff(x) == x*Cos(x) + Sin(x)

def test_trig_diff_intorder():
    x = Symbol('x')
    assert Sin(x).diff(x,2) == -Sin(x)
    assert Sin(x).diff(x,3) == -Cos(x)
    assert Sin(x).diff(x,4) == Sin(x)
    assert Sin(x).diff(x,4000000) == Sin(x)

    assert Cos(x).diff(x,2) == -Cos(x)
    assert Cos(x).diff(x,3) == Sin(x)
    assert Cos(x).diff(x,4) == Cos(x)
    assert Cos(x).diff(x,4000000) == Cos(x)

    assert Tan(x).diff(x,2) == Tan(x).diff(x).diff(x)
    assert Cot(x).diff(x,2) == Cot(x).diff(x).diff(x)

def test_trig_diff_symorder():
    x = Symbol('x')
    n = Symbol('n')
    assert Sin(x).diff(x,n) == Sin(x + pi*n/2),`Sin(x).diff(x,n)`
    assert Cos(x).diff(x,n) == Cos(x + pi*n/2),`Cos(x).diff(x,n)`
    assert Tan(x).diff(x,n) == (FD[0]**n)(Tan)(x),`Tan(x).diff(x,n)`
    assert Cot(x).diff(x,n) == (FD[0]**n)(Cot)(x),`Cot(x).diff(x,n)`
