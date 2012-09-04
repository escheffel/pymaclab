
from sympycore import CommutativeRing, heads, Expr, core

def test_add():
    x,y,z = map(CommutativeRing, 'xyz')
    assert x+y == y+x
    assert x+2 == 2+x
    assert str(x+y) in ['x + y', 'y + x'],str(x+y)
    assert str(2+x) in ['2 + x', 'x + 2'], str(2+x)

def test_mul():
    x,y,z = map(CommutativeRing, 'xyz')
    assert x*y == y*x
    assert x*2 == 2*x
    assert str(x*y) in ['x*y', 'y*x'],str(x*y)
    assert str(2*x) == '2*x', str(2*x)
    assert str(x*2) == '2*x', str(x*2)

    assert y*(2*x)==(y*2)*x
    assert str(y*(2*x)) in ['2*x*y','2*y*x'], str(y*(2*x))

    assert x*(y*z)==(x*y)*z
    assert str(x*y*z) in ['x*y*z','y*x*z'], str(x*y*z)

    assert (2*x)*(y*z)==(x*y)*(2*z)
    assert str((2*x)*(y*z)) in ['2*x*y*z', '2*y*x*z'], str((2*x)*(y*z))
    
    assert (y*x)*(y*z)==x*y*z*y
    assert (y*y)*(x*z)==x*y**2*z
    assert str(x*y*z*y) in ['y**2*x*z'], str(x*y*z*y)

    assert (y*y)*y==y*(y*y)
    assert (y*y)*(y*y)==y*(y*y)*y
    assert str(y*y*y)=='y**3',str(y*y*y)

    assert (x*x)*(y*y)==x*(x*y)*y
    assert str((x*x)*(y*y))=='y**2*x**2',str((x*x)*(y*y))

def test_diff():
    x,y,z = map(CommutativeRing, 'xyz')
    assert x.diff(x, order=0)==x, str(x.diff(x, order=0))
    assert x.diff(x)==1, str(x.diff(x))
    assert x.diff('x')==1, str(x.diff('x'))

def test_to():
    x,y,z = map(CommutativeRing, 'xyz')
    assert x.to(heads.EXP_COEFF_DICT).data==(('x',), {core.IntegerList([1]):1})
    assert x.to(heads.EXP_COEFF_DICT, x, y).data==(('x','y'), {core.IntegerList([1,0]):1})
    assert (x+y).to(heads.EXP_COEFF_DICT).data==(('x', 'y'), {core.IntegerList([1,0]):1,core.IntegerList([0,1]):1})
    assert (x+y).to(heads.EXP_COEFF_DICT, x,y).data==(('x', 'y'), {core.IntegerList([1,0]):1,core.IntegerList([0,1]):1})
    assert (x+y).to(heads.EXP_COEFF_DICT, x).data==(('x', 'y'), {core.IntegerList([1,0]):1,core.IntegerList([0,1]):1})
    assert (x*2).to(heads.EXP_COEFF_DICT).data==(('x',), {core.IntegerList([1]):2}), (x*2).to(heads.EXP_COEFF_DICT).data

commutative_operations_results = '''\
(1)/(0):zoo
(x)/(0):zoo
(2*x)/(0):zoo

(0)+(2):2
(0)-(2):-2
(0)*(2):0
(0)/(2):0
(0)**(2):0
(0)+(x):x
(0)-(x):-x
(0)*(x):0
(0)/(x):0
(0)**(x):0**x
(0)+(2*x):2*x
(0)-(2*x):-2*x
(0)*(2*x):0
(0)/(2*x):0
(0)**(2*x):0**(2*x)
(0)+(y + x):y + x
(0)-(y + x):-y - x
(0)*(y + x):0
(0)/(y + x):0
(0)**(y + x):0**(y + x)
(0)+(x**2):x**2
(0)-(x**2):-x**2
(0)*(x**2):0
(0)/(x**2):0
(0)**(x**2):0**x**2
(0)+(y*x):y*x
(0)-(y*x):-y*x
(0)*(y*x):0
(0)/(y*x):0
(0)**(y*x):0**(y*x)
(0)+(Foo(x)):Foo(x)
(0)-(Foo(x)):-Foo(x)
(0)*(Foo(x)):0
(0)/(Foo(x)):0
(0)**(Foo(x)):0**Foo(x)
(0)+(Foo):Foo
(0)-(Foo):-Foo
(0)*(Foo):0
(0)/(Foo):0
(0)**(Foo):0**Foo

(1)+(2):3
(1)-(2):-1
(1)*(2):2
(1)/(2):1/2
(1)**(2):1
(1)+(x):x + 1
(1)-(x):-x + 1
(1)*(x):x
(1)/(x):1/x
(1)**(x):1
(1)+(2*x):2*x + 1
(1)-(2*x):-2*x + 1
(1)*(2*x):2*x
(1)/(2*x):1/2/x
(1)**(2*x):1
(1)+(y + x):y + x + 1
(1)-(y + x):-y - x + 1
(1)*(y + x):y + x
(1)/(y + x):1/(y + x)
(1)**(y + x):1
(1)+(x**2):1 + x**2
(1)-(x**2):1 - x**2
(1)*(x**2):x**2
(1)/(x**2):1/x**2
(1)**(x**2):1
(1)+(y*x):1 + y*x
(1)-(y*x):1 - y*x
(1)*(y*x):y*x
(1)**(y*x):1
+(2):2
-(2):-2
(2)+(1):3
(2)-(1):1
(2)-(0):2
(2)+(0):2
(2)*(1):2
(2)*(0):0
(2)/(1):2
(2)/(0):zoo
(2)**(1):2
(2)**(0):1
(2)+(2):4
(2)-(2):0
(2)*(2):4
(2)/(2):1
(2)**(2):4
(2)+(x):2 + x
(2)-(x):-x + 2
(2)*(x):2*x
(2)/(x):2/x
(2)**(x):2**x
(2)+(2*x):2*x + 2
(2)-(2*x):-2*x + 2
(2)*(2*x):4*x
(2)/(2*x):1/x
(2)**(2*x):2**(2*x)
(2)+(y + x):y + x + 2
(2)-(y + x):-y - x + 2
(2)*(y + x):2*y + 2*x
(2)/(y + x):2/(y + x)
(2)**(y + x):2**(y + x)
(2)+(x**2):2 + x**2
(2)-(x**2):2 - x**2
(2)*(x**2):2*x**2
(2)/(x**2):2/x**2
(2)**(x**2):2**x**2
(2)+(y*x):2 + y*x
(2)-(y*x):2 - y*x
(2)*(y*x):2*y*x
(2)**(y*x):2**(y*x)
+(x):x
-(x):-x
(x)+(1):x + 1
(x)-(1):x - 1
(x)*(1):x
(x)/(1):x
(x)**(1):x
(x)+(2):x + 2
(x)-(2):x - 2
(x)*(2):2*x
(x)/(2):1/2*x
(x)**(2):x**2
(x)+(x):2*x
(x)-(x):0
(x)*(x):x**2
(x)/(x):1
(x)**(x):x**x
(x)+(2*x):3*x
(x)-(2*x):-x
(x)*(2*x):2*x**2
(x)/(2*x):1/2
(x)**(2*x):x**(2*x)
(x)+(y + x):y + 2*x
(x)-(y + x):-y
(x)*(y + x):x*(y + x)
(x)/(y + x):x/(y + x)
(x)**(y + x):x**(y + x)
(x)+(x**2):x + x**2
(x)-(x**2):x - x**2
(x)*(x**2):x**3
(x)/(x**2):1/x
(x)**(x**2):x**x**2
(x)+(y*x):x + y*x
(x)-(y*x):x - y*x
(x)*(y*x):y*x**2
(x)/(y*x):1/y
(x)**(y*x):x**(y*x)
+(2*x):2*x
-(2*x):-2*x
(2*x)+(1):2*x + 1
(2*x)-(1):2*x - 1
(2*x)*(1):2*x
(2*x)/(1):2*x
(2*x)**(1):2*x
(2*x)+(2):2*x + 2
(2*x)-(2):2*x - 2
(2*x)*(2):4*x
(2*x)/(2):x
(2*x)**(2):4*x**2
(2*x)+(x):3*x
(2*x)-(x):x
(2*x)*(x):2*x**2
(2*x)/(x):2
(2*x)**(x):(2*x)**x
(2*x)+(2*x):4*x
(2*x)-(2*x):0
(2*x)*(2*x):4*x**2
(2*x)/(2*x):1
(2*x)**(2*x):(2*x)**(2*x)
(2*x)+(y + x):y + 3*x
(2*x)-(y + x):-y + x
(2*x)*(y + x):2*x*(y + x)
(2*x)/(y + x):2*x/(y + x)
(2*x)**(y + x):(2*x)**(y + x)
(2*x)+(x**2):2*x + x**2
(2*x)-(x**2):2*x - x**2
(2*x)*(x**2):2*x**3
(2*x)/(x**2):2/x
(2*x)**(x**2):(2*x)**x**2
(2*x)+(y*x):2*x + y*x
(2*x)-(y*x):2*x - y*x
(2*x)*(y*x):2*y*x**2
(2*x)/(y*x):2/y
(2*x)**(y*x):(2*x)**(y*x)
+(y + x):y + x
-(y + x):-y - x
(y + x)+(1):y + x + 1
(y + x)-(1):y + x - 1
(y + x)*(1):y + x
(y + x)/(1):y + x
(y + x)**(1):y + x
(y + x)+(2):y + x + 2
(y + x)-(2):y + x - 2
(y + x)*(2):2*y + 2*x
(y + x)/(2):1/2*y + 1/2*x
(y + x)**(2):(y + x)**2
(y + x)+(x):y + 2*x
(y + x)-(x):y
(y + x)*(x):x*(y + x)
(y + x)/(x):1/x*(y + x)
(y + x)**(x):(y + x)**x
(y + x)+(2*x):y + 3*x
(y + x)-(2*x):y - x
(y + x)*(2*x):2*x*(y + x)
(y + x)/(2*x):1/2/x*(y + x)
(y + x)**(2*x):(y + x)**(2*x)
(y + x)+(y + x):2*y + 2*x
(y + x)-(y + x):0
(y + x)*(y + x):(y + x)**2
(y + x)/(y + x):1
(y + x)**(y + x):(y + x)**(y + x)
(y + x)+(x**2):y + x + x**2
(y + x)-(x**2):y + x - x**2
(y + x)*(x**2):x**2*(y + x)
(y + x)/(x**2):1/x**2*(y + x)
(y + x)**(x**2):(y + x)**x**2
(y + x)+(y*x):y + x + y*x
(y + x)-(y*x):y + x - y*x
(y + x)*(y*x):y*x*(y + x)
(y + x)/(y*x):1/y/x*(y + x)
(y + x)**(y*x):(y + x)**(y*x)

+(x**2):x**2
-(x**2):-x**2
(x**2)+(1):1 + x**2
(x**2)-(1):-1 + x**2
(x**2)*(1):x**2
(x**2)/(1):x**2
(x**2)**(1):x**2
(x**2)+(2):2 + x**2
(x**2)-(2):-2 + x**2
(x**2)*(2):2*x**2
(x**2)/(2):1/2*x**2
(x**2)**(2):x**4
(x**2)+(x):x + x**2
(x**2)-(x):-x + x**2
(x**2)*(x):x**3
(x**2)/(x):x
(x**2)**(x):(x**2)**x
(x**2)+(2*x):2*x + x**2
(x**2)-(2*x):-2*x + x**2
(x**2)*(2*x):2*x**3
(x**2)/(2*x):1/2*x
(x**2)**(2*x):(x**2)**(2*x)
(x**2)+(y + x):y + x + x**2
(x**2)-(y + x):-y - x + x**2
(x**2)*(y + x):x**2*(y + x)
(x**2)/(y + x):x**2/(y + x)
(x**2)**(y + x):(x**2)**(y + x)
(x**2)+(x**2):2*x**2
(x**2)-(x**2):0
(x**2)*(x**2):x**4
(x**2)/(x**2):1
(x**2)**(x**2):(x**2)**x**2
(x**2)+(y*x):y*x + x**2
(x**2)-(y*x):-y*x + x**2
(x**2)*(y*x):y*x**3
(x**2)/(y*x):x/y;1/y*x
(x**2)**(y*x):(x**2)**(y*x)
+(y*x):y*x
-(y*x):-y*x
(y*x)+(1):1 + y*x
(y*x)-(1):-1 + y*x
(y*x)*(1):y*x
(y*x)/(1):y*x
(y*x)**(1):y*x
(y*x)+(2):2 + y*x
(y*x)-(2):-2 + y*x
(y*x)*(2):2*y*x
(y*x)/(2):1/2*y*x
(y*x)**(2):y**2*x**2
(y*x)+(x):x + y*x
(y*x)-(x):-x + y*x
(y*x)*(x):y*x**2
(y*x)/(x):y
(y*x)**(x):(y*x)**x
(y*x)+(2*x):2*x + y*x
(y*x)-(2*x):-2*x + y*x
(y*x)*(2*x):2*y*x**2
(y*x)**(2*x):(y*x)**(2*x)
(y*x)+(y + x):y + x + y*x
(y*x)-(y + x):-y - x + y*x
(y*x)*(y + x):y*x*(y + x)
(y*x)**(y + x):(y*x)**(y + x)
(y*x)+(x**2):y*x + x**2
(y*x)-(x**2):y*x - x**2
(y*x)*(x**2):y*x**3
(y*x)/(x**2):y/x
(y*x)**(x**2):(y*x)**x**2
(y*x)+(y*x):2*y*x
(y*x)-(y*x):0
(y*x)*(y*x):y**2*x**2
(y*x)/(y*x):1
(y*x)**(y*x):(y*x)**(y*x)
(1)/(y*x):1/y/x
(2)/(y*x):2/y/x
(y*x)/(2*x):1/2*y
(y*x)/(y + x):y*x/(y + x)
(1)+(Foo(x)):Foo(x) + 1
(1)-(Foo(x)):-Foo(x) + 1;1 - Foo(x)
(1)*(Foo(x)):Foo(x)
(1)/(Foo(x)):1/Foo(x)
(1)**(Foo(x)):1
(2)+(Foo(x)):2 + Foo(x)
(2)-(Foo(x)):-Foo(x) + 2
(2)*(Foo(x)):2*Foo(x)
(2)/(Foo(x)):2/Foo(x)
(2)**(Foo(x)):2**Foo(x)
(x)+(Foo(x)):x + Foo(x)
(x)-(Foo(x)):-Foo(x) + x;x - Foo(x)
(x)*(Foo(x)):x*Foo(x)
(x)/(Foo(x)):x/Foo(x)
(x)**(Foo(x)):x**Foo(x)
(2*x)+(Foo(x)):2*x + Foo(x)
(2*x)-(Foo(x)):2*x - Foo(x)
(2*x)*(Foo(x)):2*x*Foo(x)
(2*x)/(Foo(x)):2*x/Foo(x)
(2*x)**(Foo(x)):(2*x)**Foo(x)
(y + x)+(Foo(x)):y + x + Foo(x)
(y + x)-(Foo(x)):y + x - Foo(x)
(y + x)*(Foo(x)):Foo(x)*(y + x)
(y + x)/(Foo(x)):1/Foo(x)*(y + x)
(y + x)**(Foo(x)):(y + x)**Foo(x)
(x**2)+(Foo(x)):Foo(x) + x**2
(x**2)-(Foo(x)):-Foo(x) + x**2
(x**2)*(Foo(x)):Foo(x)*x**2
(x**2)/(Foo(x)):x**2/Foo(x)
(x**2)**(Foo(x)):(x**2)**Foo(x)
(y*x)+(Foo(x)):Foo(x) + y*x
(y*x)-(Foo(x)):-Foo(x) + y*x
(y*x)*(Foo(x)):y*x*Foo(x)
(y*x)/(Foo(x)):y*x/Foo(x)
(y*x)**(Foo(x)):(y*x)**Foo(x)
+(Foo(x)):Foo(x)
-(Foo(x)):-Foo(x)
(Foo(x))+(1):Foo(x) + 1
(Foo(x))-(1):Foo(x) - 1
(Foo(x))*(1):Foo(x)
(Foo(x))/(1):Foo(x)
(Foo(x))**(1):Foo(x)
(Foo(x))+(2):Foo(x) + 2
(Foo(x))-(2):Foo(x) - 2
(Foo(x))*(2):2*Foo(x)
(Foo(x))/(2):1/2*Foo(x)
(Foo(x))**(2):Foo(x)**2
(Foo(x))+(x):Foo(x) + x
(Foo(x))-(x):-x + Foo(x); Foo(x) - x
(Foo(x))*(x):Foo(x)*x
(Foo(x))/(x):Foo(x)/x
(Foo(x))**(x):Foo(x)**x
(Foo(x))+(2*x):2*x + Foo(x)
(Foo(x))-(2*x):-2*x + Foo(x)
(Foo(x))*(2*x):2*Foo(x)*x
(Foo(x))/(2*x):1/2*Foo(x)/x
(Foo(x))**(2*x):Foo(x)**(2*x)
(Foo(x))+(y + x):y + x + Foo(x)
(Foo(x))-(y + x):-y - x + Foo(x)
(Foo(x))*(y + x):Foo(x)*(y + x)
(Foo(x))/(y + x):Foo(x)/(y + x)
(Foo(x))**(y + x):Foo(x)**(y + x)
(Foo(x))+(x**2):Foo(x) + x**2
(Foo(x))-(x**2):Foo(x) - x**2
(Foo(x))*(x**2):Foo(x)*x**2
(Foo(x))/(x**2):Foo(x)/x**2
(Foo(x))**(x**2):Foo(x)**x**2
(Foo(x))+(y*x):Foo(x) + y*x
(Foo(x))-(y*x):Foo(x) - y*x
(Foo(x))*(y*x):y*x*Foo(x)
(Foo(x))/(y*x):1/y*Foo(x)/x
(Foo(x))**(y*x):Foo(x)**(y*x)
(Foo(x))+(Foo(x)):2*Foo(x)
(Foo(x))-(Foo(x)):0
(Foo(x))*(Foo(x)):Foo(x)**2
(Foo(x))/(Foo(x)):1
(Foo(x))**(Foo(x)):Foo(x)**Foo(x)

(1)+(Foo):1 + Foo
(1)-(Foo):1 - Foo
(1)*(Foo):Foo
(1)/(Foo):1/Foo
(1)**(Foo):1
(2)+(Foo):2 + Foo
(2)-(Foo):2 - Foo
(2)*(Foo):2*Foo
(2)/(Foo):2/Foo
(2)**(Foo):2**Foo
(x)+(Foo):x + Foo
(x)-(Foo):x - Foo
(x)*(Foo):x*Foo
(x)/(Foo):x/Foo
(x)**(Foo):x**Foo
(2*x)+(Foo):2*x + Foo
(2*x)-(Foo):2*x - Foo
(2*x)*(Foo):2*x*Foo
(2*x)/(Foo):2*x/Foo
(2*x)**(Foo):(2*x)**Foo
(y + x)+(Foo):y + x + Foo
(y + x)-(Foo):y + x - Foo
(y + x)*(Foo):(y + x)*Foo
(y + x)/(Foo):(y + x)/Foo
(y + x)**(Foo):(y + x)**Foo
(x**2)+(Foo):x**2 + Foo
(x**2)-(Foo):x**2 - Foo
(x**2)*(Foo):x**2*Foo
(x**2)/(Foo):x**2/Foo
(x**2)**(Foo):(x**2)**Foo
(y*x)+(Foo):y*x + Foo
(y*x)-(Foo):y*x - Foo
(y*x)*(Foo):y*x*Foo
(y*x)/(Foo):y*x/Foo
(y*x)**(Foo):(y*x)**Foo
(Foo(x))+(Foo):Foo(x) + Foo
(Foo(x))-(Foo):Foo(x) - Foo
(Foo(x))*(Foo):Foo(x)*Foo
(Foo(x))/(Foo):Foo(x)/Foo
(Foo(x))**(Foo):Foo(x)**Foo
+(Foo):Foo
-(Foo):-Foo
(Foo)+(1):1 + Foo
(Foo)-(1):-1 + Foo
(Foo)*(1):Foo
(Foo)/(1):Foo
(Foo)**(1):Foo
(Foo)+(2):2 + Foo
(Foo)-(2):-2 + Foo
(Foo)*(2):2*Foo
(Foo)/(2):1/2*Foo
(Foo)**(2):Foo**2
(Foo)+(x):x + Foo
(Foo)-(x):-x + Foo
(Foo)*(x):x*Foo
(Foo)/(x):1/x*Foo
(Foo)**(x):Foo**x
(Foo)+(2*x):2*x + Foo
(Foo)-(2*x):-2*x + Foo
(Foo)*(2*x):2*x*Foo
(Foo)/(2*x):1/2/x*Foo
(Foo)**(2*x):Foo**(2*x)
(Foo)+(y + x):y + x + Foo
(Foo)-(y + x):-y - x + Foo
(Foo)*(y + x):(y + x)*Foo
(Foo)/(y + x):1/(y + x)*Foo
(Foo)**(y + x):Foo**(y + x)
(Foo)+(x**2):x**2 + Foo
(Foo)-(x**2):-x**2 + Foo
(Foo)*(x**2):x**2*Foo
(Foo)/(x**2):1/x**2*Foo
(Foo)**(x**2):Foo**x**2
(Foo)+(y*x):y*x + Foo
(Foo)-(y*x):-y*x + Foo
(Foo)*(y*x):y*x*Foo
(Foo)/(y*x):1/y/x*Foo
(Foo)**(y*x):Foo**(y*x)
(Foo)+(Foo(x)):Foo(x) + Foo
(Foo)-(Foo(x)):-Foo(x) + Foo
(Foo)*(Foo(x)):Foo(x)*Foo
(Foo)/(Foo(x)):1/Foo(x)*Foo
(Foo)**(Foo(x)):Foo**Foo(x)
(Foo)+(Foo):2*Foo
(Foo)-(Foo):0
(Foo)*(Foo):Foo**2
(Foo)/(Foo):1
(Foo)**(Foo):Foo**Foo
(Foo(x))+(0):Foo(x)

(x)+(0):x
(x)-(0):x
(x)*(0):0
(x)**(0):1
(2*x)+(0):2*x
(2*x)-(0):2*x
(2*x)*(0):0
(2*x)**(0):1
(y + x)+(0):y + x
(y + x)-(0):y + x
(y + x)*(0):0
(y + x)/(0):zoo
(y + x)**(0):1
(x**2)+(0):x**2
(x**2)-(0):x**2
(x**2)*(0):0
(x**2)/(0):zoo
(x**2)**(0):1
(y*x)+(0):y*x
(y*x)-(0):y*x
(y*x)*(0):0
(y*x)/(0):zoo
(y*x)**(0):1

(Foo(x))-(0):Foo(x)
(Foo(x))*(0):0
(Foo(x))/(0):zoo
(Foo(x))**(0):1
(Foo)+(0):Foo
(Foo)-(0):Foo
(Foo)*(0):0
(Foo)/(0):zoo
(Foo)**(0):1

'''

def test_commutative_operations():
    Ring = CommutativeRing
    x,y,z = map(Ring, 'xyz')
    FRing = Ring.get_function_algebra()
    def Foo(x):
        return Ring(heads.APPLY, (foo, (Ring(x),)))
    foo = FRing(heads.CALLABLE, Foo)
    
    operands = [1,
                Ring(heads.NUMBER, 2),
                0,
                x,
                Ring(heads.TERM_COEFF, (x, 2)), # 2*x
                Ring(heads.TERM_COEFF_DICT, {x:1, y:1}), # x+y
                Ring(heads.POW, (x, 2)), # x**2
                Ring(heads.BASE_EXP_DICT, {x:1, y:1}), # x*y
                Foo(x), # Foo(x),
                foo,    # function foo
                ]
    unary_operations = ['+', '-']
    binary_operations = ['+', '-', '*', '/', '**']
    results = {}
    for line in commutative_operations_results.split('\n'):
        line = line.strip()
        if ':' not in line: continue
        expr, result = line.split(':')
        for e in expr.split(';'):
            e = e.strip()
            results[e] = [r.strip() for r in result.split(';')]

    for op1 in operands:
        if isinstance(op1, Expr):
            for op in unary_operations:
                expr = '%s(%s)' % (op, op1) 

                try:
                    result = str(eval('%s(op1)' % (op)))
                except Exception, msg:
                    print  expr,'failed with %s' % (msg)
                    raise
                
                if expr not in results:
                    print '%s:%s' % (expr, result)
                    continue
                assert result in results[expr], `results[expr], result`
        for op2 in operands:
            if not (isinstance(op1, Expr) or isinstance(op2, Expr)):
                continue
            for op in binary_operations:
                expr = '(%s)%s(%s)' % (op1, op, op2)

                try:
                    result = str(eval('(op1)%s(op2)' % op))
                except Exception, msg:
                    print  expr,'failed with %s' % (msg)
                    raise
                
                if expr not in results:
                    print '%s:%s' % (expr, result)
                    continue
                assert result in results[expr], `results[expr], result, op1, op2`
