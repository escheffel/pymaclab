
from sympycore import Ring, heads, core

P = Ring.Polynom

def test_zero():
    x,y,z = map(Ring, 'xyz')
    zero =  P(variables=(x,y))
    assert zero.pair==(heads.EXP_COEFF_DICT, ((x,y), {})), `zero.pair`
    assert zero==0
    zero = ((zero + x) - x)
    assert zero.pair==(heads.EXP_COEFF_DICT, ((x,y), {})), `zero.pair`
    assert zero==0,`zero`

    
def test_one():
    x,y,z = map(Ring, 'xyz')
    one = Ring.Polynom(1,variables=(x,y))
    assert one.pair == (heads.EXP_COEFF_DICT, ((x, y), {core.IntegerList([0,0]):1})), `one.pair`
    assert one==1
    one = Ring.Polynom(x,variables=(x,y))/x
    assert one.pair == (heads.EXP_COEFF_DICT, ((), {core.IntegerList([]):1})), `one.pair`
    assert one==1

def test_simple():
    x,y,z = map(Ring, 'xyz')
    
