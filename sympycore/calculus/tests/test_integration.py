from sympycore.calculus import *

x, y = map(Symbol, 'xy')

def test_integrate():
    assert Number(3).integrate(x) == 3*x
    assert x.integrate(x) == x**2 / 2
    assert (2*x).integrate(x) == x**2
    assert (x*y).integrate(x) == y * x**2 / 2
    assert (x*y**2).integrate(x) == y**2 * x**2 / 2
    r1 = (1 + 4*y*x + 3*x**2).integrate(x)
    r2 = x + 2*y*x**2 + x**3
    assert r1 == r2
    failures = [x**x, x**(1+y+5**x), 1/(2+x), x/(2+3*x+5), (3*x+2*y)*x][:-1]
    for f in failures:
        try:
            r = f.integrate(x)
            assert 0, ("integration of %s expected to fail but got %s" % (f, r))
        except NotImplementedError:
            pass

def test_integrate_defined():
    for g in [Number(4), 2*x, 6*x**2, 2*x+100*x**3]:
        f = g.integrate(x)
        f1 = f.subs(x,5)-f.subs(x,2)
        f2 = g.integrate((x,2,5))
        assert f2==f1,`str(g), str(f), str(f1), str(f2)`
