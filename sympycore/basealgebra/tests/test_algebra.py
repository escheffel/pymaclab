
from sympycore import Algebra, Verbatim

def test_basic():
    a = Algebra('a')
    assert str(a)=='a'
    assert repr(a)=="Algebra('a')",`repr(a)`
    assert bool(a)==True
    assert a.as_algebra(Algebra)==a
    assert a.as_algebra(Verbatim)==Verbatim('a')
    assert Verbatim('a').as_algebra(Algebra)==a
