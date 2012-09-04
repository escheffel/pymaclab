from sympycore.calculus import *
A = Calculus.convert

x = Symbol('x')
y = Symbol('y')

def test_constructor():
    assert repr(A(1.2))=="Calculus('1.2')"
    assert repr(A('1.2'))=="Calculus('1.2')"

def test_linear():
    assert x + x == 2*x
    assert (2*x) / 2 == x
    assert (x/2)*2 == x
    assert 3*x + 4*x == 7*x
    assert 2*x/3 + x/3 == x
    assert I*x + 3*x == (3+I)*x

def test_powers():
    assert Number(3) ** Number(-2) == Number(1,9)
    assert Number(4) ** Number(1,2) == Number(2)
    assert 2**Number(3,2) * 2**Number(1,2) == 4, `2**Number(3,2) * 2**Number(1,2)`
    assert I**Number(3,2) * I**Number(1,2) == -1
    assert str(Number(2) ** Number(1,2)) == '2**(1/2)', str(Number(2) ** Number(1,2))
    assert A('4/9') ** A('1/2') == A('2/3')
    assert A('4/7') ** A('1/2') == A('7**(-1/2)')*2
    assert A('7/4') ** A('1/2') == A('7**(1/2)')/2
    assert str((I*x)**2)==str('-x**2')
    assert str((-I*x)**2)==str('-x**2')
    assert str((I*x)**4)==str('x**4')


def test_has_symbol():
    assert (1 + Cos(1+2**x)).has_symbol(x)
    assert (y + Cos(1+2**x)).symbols == set([x, y])

def test_subs():
    assert (oo*x + oo*y).subs(y,x) == oo*x
    assert (oo*x + oo*y).subs(y,-x) == undefined

def test_diff():
    x = Symbol('x')
    assert Number(2).diff(x) == 0
    assert x.diff(x) == 1
    assert (3*x).diff(x) == 3
    assert (3*x+1).diff(x) == 3
    assert (2*x**2 + x).diff(x) == 4*x + 1
    assert (x**3).diff(x) == 3*x**2
    assert ((1+x)*(2+x)*(3+x)**2).diff(x) in [((1+x)*(2+x)*(3+x)*2 + (1+x)*(3+x)**2 + (2+x)*(3+x)**2),
                                              (1 + x)*(2 + x)*(6 + 2*x) + (3 + x)**2*(1 + x) + (3 + x)**2*(2 + x)]

def test_subs():
    x = Symbol('x')
    y = Symbol('y')
    assert 2*x | (x, pi) == 2*pi
    assert x*y | {x:2, y:pi} == 2*pi
