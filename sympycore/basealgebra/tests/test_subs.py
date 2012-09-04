
from sympycore import CollectingField as Algebra

Symbol = Algebra.Symbol
Number = Algebra.Number
Add = Algebra.Add
Mul = Algebra.Mul
Pow = Algebra.Pow
Terms = Algebra.Terms
Factors = Algebra.Factors

x,y,z = map(Symbol, 'xyz')

def test_simple():
    assert x.subs(x,y)==y
    assert x.subs(z,y)==x
    assert (2*x).subs(x,y)==2*y
    assert (x+z).subs(x,y)==y+z
    assert (x**2).subs(x,y)==y**2
    assert (x*z).subs(x,y)==y*z
    assert (x**2*z).subs(x,y)==y**2*z
    assert (x + x*z).subs(x,y)==y+y*z
    
def test_repeated():
    assert x.subs([(x,y)])==y
    assert (x + x*z).subs([(x,y),(z,x)])==y+y*x,`(x + x*z).subs([(x,y),(z,x)])`
    assert (x+y).subs({x:1,y:2*z})==1+2*z,`((x+y).subs({x:1,y:2*z})).pair, (1+2*z).pair`
