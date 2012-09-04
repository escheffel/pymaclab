""" Provides calculus support.
"""
__docformat__ = "restructuredtext"

from .algebra import Calculus, I
from .infinity import  oo, undefined, moo, zoo
from .functions import (Exp, Log, Sqrt, Sin, Cos, Tan, Cot, pi, E, gamma,
    Sign, Mod, Ln, Factorial, ArcSin)
from .functions import CalculusFunctionRing, CalculusDifferentialRing, CalculusOperatorRing

def diff(expr, symbol, order=1):
    return expr.diff(symbol, order)

from .relational import Assumptions

Symbol = Calculus.Symbol

def Number(num, denom=None):
    n = Calculus.Number(Calculus.convert_coefficient(num))
    if denom is None:
        return n
    return n / denom

def Rational(num, denom):
    #XXX: use Div
    return Calculus.Number(Calculus.convert_coefficient(num)) / denom

Add = Calculus.Add
Mul = Calculus.Mul
Pow = Calculus.Pow

from ..arithmetic.number_theory import factorial as _factorial
def factorial(n):
    return Number(_factorial(n))
    
Polynom = Calculus.Polynom
