from sympycore.arithmetic.evalf import *
from sympycore.arithmetic.evalf import mpmath, compile_mpmath
from sympycore.calculus import Symbol, I, Number, Exp, Sin, Cos, E, pi
import math
import cmath

def test_evalf():
    expr1 = Number(1)/3
    expr2 = Sin(E)**2 + Cos(E)**2 - 1
    expr3 = Exp(I) - Cos(1) - I*Sin(1)
    assert abs(pi.evalf(15) - math.pi) < 1e-14, str( abs(pi.evalf(15) - math.pi))
    assert abs(expr1.evalf(30) - expr1) < 1e-29
    assert abs(expr2.evalf(30)) < 1e-29, `abs(expr2.evalf(30))`
    assert abs(expr2.evalf(100)) < 1e-99
    assert abs(expr2.evalf(300)) < 1e-99
    #assert abs(expr3.evalf(20)) < 1e-19

def test_compiled():
    x = Symbol('x')
    y = Symbol('y')
    f1 = compile_mpmath([], Exp(2))
    f2 = compile_mpmath('x', Exp(x))
    f3 = compile_mpmath(['x', 'y'], Cos(x)+Sin(y)*I)
    mpmath.mp.dps = 15
    assert abs(f1() - math.exp(2)) < 1e-14
    assert abs(f2(2) - math.exp(2)) < 1e-14
    assert abs(f3(3,4) - (cmath.cos(3)+cmath.sin(4)*1j)) < 1e-14
