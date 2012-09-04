
from sympycore import Ring, core

def test_add():
    x,y,z = map(Ring, 'xyz')
    assert Ring.Add(x,y) == x+y

def test_div():
    x,y,z = map(Ring, 'xyz')
    assert Ring.Div(x,y) == x/y, str(Ring.Div(x,y))

def test_polynom():
    x,y,z = map(Ring, 'xyz')
    assert Ring.Polynom(x).data==(('x',), {core.IntegerList([1]):1})
    assert Ring.Polynom(x, variables=(x,y)).data==(('x','y'), {core.IntegerList([1,0]):1})
