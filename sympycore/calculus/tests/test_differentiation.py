from sympycore.calculus import *

x = Symbol('x')
n = Symbol('n')
M = 10**6

def test_diff_poly():
    assert diff(3**(pi+5), x) == 0
    assert diff(x, x) == 1
    assert diff(x, x, 2) == 0
    assert diff(3*x**2 + 5*x, x) == 6*x + 5
    assert diff(3*x**2 + 5*x, x, 2) == 6
    assert diff(x**M + 2*x + 3, x, M+1) == 0

def test_diff_powers():
    assert diff(Exp(x), x) == Exp(x)
    assert diff(Exp(x+1), x) == Exp(x+1)
    assert diff(Exp(x), x, M) == Exp(x)
    assert diff(Exp(x+1), x, M) == Exp(x+1)
    assert diff(Exp(3*x+1), x) == 3*Exp(3*x+1)
    assert diff(2**(3*x+1), x) == 3*2**(3*x+1)*Log(2)
    assert diff(2**(3*x+1), x, 1000) == 3**1000 * 2**(3*x+1) * Log(2)**1000
    assert diff((2*x+1)**(3*x), x) == (2*x+1)**(3*x) * (6*x/(2*x+1)+3*Log(2*x+1))

def test_diff_products():
    assert diff((1+x)*(2+x), x) == 3+2*x
    assert diff((1+x)**3*(2+x)**2, x) in [(2*(1+x)**3*(2+x) + 3*(1+x)**2 * (2+x)**2),
                                          (1 + x)**3*(4 + 2*x) + 3*(1 + x)**2*(2 + x)**2]
    assert diff(pi**2 * 2**Number(1,2) * x, x) == pi**2 * 2**Number(1,2)

def test_diff_Log():
    assert diff(Log(x), x) == 1/x
    assert diff(Log(3*x), x) == 1/x
    assert diff(Log(3*x+1), x) == 3/(3*x+1)
    assert diff(Log(3*x+1), x, 5) == 24 * 3**5 / (3*x+1)**5

def test_diff_trig():
    assert diff(Sin(x), x) == Cos(x)
    assert diff(Cos(x), x) == -Sin(x)
    assert diff(Sin(2*x), x) == 2*Cos(2*x)
    assert diff(Tan(x), x) == 1 + Tan(x)**2
    assert diff(Cot(x), x) == -1-Cot(x)**2
    assert diff(Sin(x), x, M) == Sin(x)
    assert diff(Sin(x), x, M+1) == Cos(x)
    assert diff(Cos(x), x, M) == Cos(x)
    assert diff(Cos(x), x, M+1) == -Sin(x)
    assert diff(Sin(3*x+1), x, 5) == 3**5 * Cos(3*x+1)
    assert diff(Cos(3*x+1), x, 5) == -3**5 * Sin(3*x+1)
    assert diff(x*Sin(x), x) == x*Cos(x) + Sin(x)
    assert diff(x*Sin(x)*Cos(x), x) == x*Cos(x)**2 + Cos(x)*Sin(x) - x*Sin(x)**2

def test_diff_symbolic_order():
    assert diff(Exp(x), x, n) == Exp(x)
    assert diff(Exp(-x), x, n) == (-1)**n * Exp(-x)
    assert diff(Cos(2*x+3), x, n) == Cos(2*x+3 + n*pi/2) * 2**n, `diff(Cos(2*x+3), x, n)`
    assert diff(Log(2*x+3), x, n) == (-1)**(n-1) * Factorial(n-1) * (3+2*x)**(-n) * 2**n
    assert diff(2**(3*x+4), x, n) == 2**(4+3*x) * 3**n * Log(2)**n
