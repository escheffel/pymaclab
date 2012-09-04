from sympycore.polynomials import *
from sympycore.arithmetic import gcd, lcm

def test_polynomials():
    x = poly([0, 1])
    assert poly([1, 2, 3]) == 3*x**2 + 2*x + 1, `poly([1,2,3]), 3*x**2+2*x+1`
    assert poly([1, 2, 3]).degree == 2
    assert poly([1, 2, 3, 0]).degree == 2
    assert poly([0]).degree == -1
    assert (x**3 + x) + (3*x + x**4) == x**4 + x**3 + 4*x
    assert poly([1, 2, 3])*2 == poly([2, 4, 6])
    assert poly([1, 2, 3])/2.0 == poly([0.5, 1.0, 1.5])
    assert poly([4,-3,5])(6) == 166
    assert poly([1,2,3])(poly([4,5,6])) == 108*x**4 + 180*x**3 + 231*x**2 + 130*x + 57
    assert (x-1)*(x+1) == x**2 - 1
    assert (x-1)*(x+1) / (x-1) == (x+1)
    assert (x**3 + x).diff() == 3*x**2 + 1
    assert str(poly([3, 0, 1])) == '3 + x**2',str(poly([3, 0, 1]))
    assert str(poly([0])) == '0'
    assert str(poly([1, 1])) == '1 + x'
    assert str((5 + 3*x) / 5) == '1 + 3/5*x', str((5 + 3*x) / 5)
    assert str(poly([4, 11, 6]) / poly([6, 12])) == '2/3 + 1/2*x'
    assert str(poly([2, 3, 4]) % poly([1, 2, 3])) == '2/3 + 1/3*x'

def test_gcd():
    p1 = poly([1, 2, 1])
    p2 = poly([2, 3, 1])
    q, r = divmod(gcd(p1, p2), poly([1,1]))
    assert q.degree == 0 and not r
    q, r = divmod(lcm(p1, p2), p1*poly([2,1]))
    assert q.degree == 0 and not r

def test_poly_of_poly():
    p1 = poly([1,2,3],'x')
    p2 = poly([4,5,6],p1)
    assert str(p2)=="4 + 5*(1 + 2*x + 3*x**2) + 6*(1 + 2*x + 3*x**2)**2",str(p2)

def test_coeff_poly():
    p1 = poly([1,2,3],'x')
    p2 = poly([4,p1],'y')
    p2 += poly([0,3],'y')
    assert str(p2)=='4 + (4 + 2*x + 3*x**2)*y', str(p2)
