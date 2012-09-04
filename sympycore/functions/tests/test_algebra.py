
from sympycore import *

def test_add():
    x = CommutativeRing('x')
    f = FunctionRing('f')
    g = FunctionRing('g')
    assert str(f+g) in ['f + g', 'g + f'], str(f+g)
    assert str((f+g)(x)) in ['f(x) + g(x)', 'g(x) + f(x)'],  str((f+g)(x))
    assert str((f+g)(x, evaluate=False)) in ['(f + g)(x)',
                                             '(g + f)(x)'], str((f+g)(x, evaluate=False))

def test_defined():
    x = Calculus('x')
    f = CalculusFunctionRing(heads.CALLABLE, Sin)
    g = CalculusFunctionRing('g')
    assert str(f(x))=='Sin(x)'
    assert str(f+g) in ['Sin + g','g + Sin'], str(f+g)
    assert str((f+g)(x)) in ['Sin(x) + g(x)','g(x) + Sin(x)'],str((f+g)(x))
    assert str((f+g)(x, evaluate=False)) in ['(Sin + g)(x)','(g + Sin)(x)'],str((f+g)(x, evaluate=False))


