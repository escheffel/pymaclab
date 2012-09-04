
from sympycore import Ring, heads, CommutativeRing

def test_symbol():
    x,y = map(Ring, 'xy')
    assert str(x)=='x'
    assert repr(x)=="Ring('x')", repr(x)

def test_number():
    n = Ring(2)
    assert str(Ring(2))=='2'
    assert repr(Ring(2))=="Ring('2')", repr(Ring(2))
    assert str(Ring(-2))=='-2'
    assert repr(Ring(-2))=="Ring('-2')", repr(Ring(-2))

def test_neg():
    x,y,z = map(Ring, 'xyz')
    assert str(-x)=='-x', str(-x)
    assert str(-Ring(2))=='-2', str(-Ring(2))
    assert str(-(x+y)) in ['-x - y', '-y - x'], str(-(x+y))
    assert str(-x*y)=='-x*y', str(-x*y)

def test_add():
    x,y,z = map(Ring, 'xyz')
    assert str(x + y) in ['x + y', 'y + x'], str(x+y)
    assert repr(x + y) in ["Ring('x + y')", "Ring('y + x')"], repr(x+y)
    assert str(x + y + z) in ['x + y + z', 'y + x + z','x + z + y'], str(x + y + z)
    assert str(x + y + 1) in ['x + y + 1', 'y + x + 1','x + 1 + y'], str(x + y + 1)
    assert str(x + y + y) in ['x + 2*y', '2*y + x'], str(x + y + y)
    assert str(x + y + x) in ['2*x + y', 'y + 2*x'], str(x + y + x)
    assert str(x + 1 + 2)=='x + 3', str(x + 1 + 2)
    assert str(x + 1 + y) in ['x + 1 + y', 'y + x + 1'], str(x + 1 + y)

    assert str(Ring(1) + Ring(2))=='3', str(Ring(1) + Ring(2))
    assert str(2 + x) == 'x + 2', str(2+x)
    assert str(2 + x + y) in ['x + 2 + y', 'y + x + 2'], str(2+x + y)

def test_sub():
    x,y,z = map(Ring, 'xyz')
    assert str(x-y) in ['x - y','-y + x'], str(x-y)
    assert str(x-y-z) in ['x - y - z', '-y + x - z', 'x - z - y'], str(x-y-z)
    assert str(x-(y+z)) in ['x - y - z', '-y + x - z','x - z - y'], str(x-(y+z))
    assert str((x-y)-z) in ['x - y - z', '-y + x - z','x - z - y'], str((x-y)-z)

def test_ncmul():
    x,y,z = map(Ring, 'xyz')
    assert str(x*2)=='2*x', str(x*2)
    assert str(x*'2')=='2*x', str(x*'2')
    assert str('2'*x)=='2*x', str('2'*x)
    assert str(x*y)=='x*y', str(x*y)
    assert str(x*y*x)=='x*y*x', str(x*y*x)
    assert str(x*(y*x))=='x*y*x', str(x*(y*x))
    assert str(x*y*x*x)=='x*y*x**2', str(x*y*x*x)
    assert str(x*(y*x)*x)=='x*y*x**2', str(x*(y*x)*x)
    assert str(x*y*(x*x))=='x*y*x**2', str(x*y*(x*x))
    assert str((x*y)*(x*x))=='x*y*x**2', str((x*y)*(x*x))

    assert str(2*x*y)=='2*x*y', str(2*x*y)
    assert str(x*2*y)=='2*x*y', str(x*2*y)
    assert str(x*y*2)=='2*x*y', str(x*y*2)
    assert str(3*x*y*2)=='6*x*y', str(3*x*y*2)
    assert str(3*(x*y)*2)=='6*x*y', str(3*(x*y)*2)
    assert str(3*(x*x)*2)=='6*x**2', str(3*(x*x)*2)

def test_pow():
    x,y,z = map(Ring, 'xyz')
    assert str(x**2)=='x**2', str(x**2)
    assert str((2*x)**2)=='4*x**2', str((2*x)**2)
    assert str((2*x*y)**2)=='4*(x*y)**2', str((2*x*y)**2)
    assert str(x*x)=='x**2', str(x*x)
    assert str(x**y)=='x**y', str(x**y)
    assert str(2**x)=='2**x', str(2**x)

    assert str((2*x)**-2)=='1/4/x**2', str((2*x)**-2)
    assert str((2*x*y)**-2)=='1/4/(x*y)**2', str((2*x*y)**-2)

    assert str(((z*y)**2) * (1/(2*z*y)**2))=='1/4', str(((z*y)**2) * (1/(2*z*y)**2))

    assert str(x*y*x**-1 * x*y*x**-1 * x*y*x**-1)=='x*y**3/x',str(x*y*x**-1 * x*y*x**-1 * x*y*x**-1)
    assert str((x*y*x**-1)**3)=='x*y**3/x', str((x*y*x**-1)**3)

def test_div():
    x,y,z = map(Ring, 'xyz')
    assert str(1/y)=='1/y', str(1/y)
    assert str('1'/y)=='1/y', str('1'/y)
    assert str(y/'2')=='1/2*y', str(y/'2')
    assert str(1/y/y)=='1/y**2', str(1/y/y)
    assert str(x/y)=='x/y', str(x/y)
    assert str(x/y/y)=='x/y**2', str(x/y/y)
    assert str(x/y/y*y)=='x/y', str(x/y/y*y)
    assert str(x/y/y*x)=='x/y**2*x', str(x/y/y*x)
    assert str(1/(x+y)) in ['1/(x + y)', '1/(y + x)'], str(1/(x+y))
    assert str(1/(x*y))=='1/y/x', str(1/(x*y))
    assert str(y/(x*y))=='1/x', str(y/(x*y))
    assert str((x*y)*(1/(x*y)))=='1', str((x*y)*(1/(x*y)))
    assert str(1/(2*x*y))=='1/2/y/x', str(1/(2*x*y))
    assert str(1/(x*y/2))=='2/y/x', str(1/(x*y/2))

    assert str(x/0)=='zoo', str(x/0)
    assert str((x+y)/0)=='zoo', `str((x+y)/0), (x+y).pair`
    assert str((2*x)/(x-x))=='zoo', str((2*x)/(x-x))
