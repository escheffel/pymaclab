
from sympycore.calculus import oo, moo, zoo, undefined, I, Calculus, pi

zero = Calculus.Number(0)
one = Calculus.Number(1)
mone = Calculus.Number(-1)
two = Calculus.Number(2)
mtwo = Calculus.Number(-2)
half = Calculus.Number(1)/2
mhalf = Calculus.Number(-1)/2
onehalf = Calculus.Number(3)/2
monehalf = Calculus.Number(-3)/2
fhalf = one/2.0

x,y,z = map(Calculus.Symbol, 'xyz')

def test_oo_product():
    assert str(oo * (x+y)) in ['oo*(x + y)', 'oo*(y + x)'], str(oo * (x+y))
    assert str(oo * (1+x+y)) in ['oo*(1 + x + y)', 'oo*(y + 1 + x)',
                                 'oo*(x + 1 + y)', 'oo*(y + x + 1)'],  str(oo * (1+x+y))
    assert str(oo * (2*x))=='oo*(x)'
    assert str(oo / (2/x))=='oo*(x)', `str(oo/(2/x))`
    assert str(oo / (2*x)) in ['oo*(x**(-1))','oo*(1/x)'],str(oo / (2*x))
    
def test_oo_sum():
    assert str(oo + (x+y)) in ['oo + (x + y)','oo + (y + x)'], str(oo + (x+y))
    assert str(oo + (2+x+y)) in ['oo + (2 + x + y)', 'oo + (y + 2 + x)',
                                 'oo + (x + 2 + y)', 'oo + (y + x + 2)'], str(oo + (2+x+y))
    assert str((2+x+y) + oo) in ['oo + (2 + x + y)', 'oo + (y + 2 + x)',
                                 'oo + (x + 2 + y)', 'oo + (y + x + 2)'], str((2+x+y) + oo)
    assert str(oo - (2+x+y)) in ['oo + (-2 - x - y)', 'oo + (-y - 2 - x)',
                                 'oo + (-x - 2 - y)', 'oo + (-y - x - 2)'], str(oo - (2+x+y))

    assert str((oo+y) + (x+y)) in ['oo + (y) + (x + y)', 'oo + (y) + (y + x)'], str((oo+y) + (x+y))
    assert str((x+y) + (oo+y)) in ['oo + (y) + (x + y)','oo + (y) + (y + x)'],  str((x+y) + (oo+y))
    assert str((oo+y) + (oo+x))=='oo + (y) + (oo + (x))', str((oo+y) + (oo+x))

    s = x+y
    s += oo
    assert str(s) in ['oo + (x + y)', 'oo + (y + x)'], str(s)

    assert str(Calculus.Add(oo, x))=='oo + (x)'
    assert str(Calculus.Add(x, oo))=='oo + (x)'
    assert str(Calculus.Add(y, oo, x))=='oo + (y) + (x)'
    assert str(Calculus.Add(y, oo, x, oo)) in ['oo + (y) + (oo + (x))', 'oo + (y) + (x) + (oo)'], str(Calculus.Add(y, oo, x, oo))

    s = x+y
    s -= -oo
    assert str(s) in ['oo + (x + y)', 'oo + (y + x)']

def test_oo_symbol():
    assert str(oo + x)=='oo + (x)'
    assert str(oo - x)=='oo + (-x)'
    assert str(x - oo)=='-oo + (x)'
    assert str(oo * x)=='oo*(x)'
    assert str(oo / x) in ['oo*(x**(-1))','oo*(1/x)'], str(oo / x)
    assert str(x/oo)=='0'
    assert str(oo ** x) in ['(oo)**x', 'oo**x'],str(oo**x)
    assert str(x**oo) in ['x**(oo)', 'x**oo'], str(x**oo)

    assert str(oo + pi)=='oo', str(oo+pi)
    assert str(oo - pi)=='oo'
    assert str(oo * pi)=='oo'
    assert str(oo / pi)=='oo'
    #assert str(oo ** pi)=='oo'
    assert str(pi - oo)=='-oo'
    assert str(pi/oo)=='0'
    #assert str(pi**oo)=='oo'

    assert oo + x == x + oo,`oo+x, x+oo`
    assert oo * x == x * oo
    assert oo + pi == pi + oo
    assert oo * pi == pi * oo

def test_moo_optable():
    assert moo + moo == moo
    assert moo + mtwo == moo
    assert moo + monehalf == moo
    assert moo + mone == moo
    assert moo + mhalf == moo
    assert moo + zero == moo
    assert moo + half == moo
    assert moo + fhalf == moo
    assert moo + one == moo
    assert moo + onehalf == moo
    assert moo + two == moo
    assert moo + oo == undefined
    assert moo + zoo == undefined
    assert moo + undefined == undefined

    assert moo - moo == undefined
    assert moo - mtwo == moo
    assert moo - monehalf == moo
    assert moo - mone == moo
    assert moo - mhalf == moo
    assert moo - zero == moo
    assert moo - half == moo
    assert moo - fhalf == moo
    assert moo - one == moo
    assert moo - onehalf == moo
    assert moo - two == moo
    assert moo - oo == moo
    assert moo - zoo == undefined
    assert moo - undefined == undefined

    assert moo * moo == oo
    assert moo * mtwo == oo
    assert moo * monehalf == oo
    assert moo * mone == oo
    assert moo * mhalf == oo
    assert moo * zero == undefined
    assert moo * half == moo
    assert moo * fhalf == moo
    assert moo * one == moo
    assert moo * onehalf == moo
    assert moo * two == moo
    assert moo * oo == moo
    assert moo * zoo == zoo
    assert moo * undefined == undefined

    assert moo / moo == undefined
    assert moo / mtwo == oo
    assert moo / monehalf == oo
    assert moo / mone == oo
    assert moo / mhalf == oo
    assert moo / zero == zoo
    assert moo / half == moo
    assert moo / fhalf == moo
    assert moo / one == moo
    assert moo / onehalf == moo
    assert moo / two == moo
    assert moo / oo == undefined
    assert moo / zoo == undefined
    assert moo / undefined == undefined

    assert moo ** moo == zero
    assert moo ** mtwo == zero
    assert moo ** monehalf == zero
    assert moo ** mone == zero
    assert moo ** mhalf == zero
    assert moo ** zero == one
    assert moo ** half == I*oo
    #assert moo ** fhalf == I*oo
    assert moo ** one == moo
    assert moo ** onehalf == -I*oo
    assert moo ** two == oo
    assert moo ** oo == zoo
    assert moo ** zoo == undefined
    assert moo ** undefined == undefined

def test_mtwo_optable():
    assert mtwo + moo == moo
    assert mtwo + oo == oo
    assert mtwo + undefined == undefined

    assert mtwo - moo == oo
    assert mtwo - oo == moo
    assert mtwo - undefined == undefined

    assert mtwo * moo == oo
    assert mtwo * oo == moo
    assert mtwo * zoo == zoo
    assert mtwo * undefined == undefined

    assert mtwo / moo == zero
    assert mhalf / zero == zoo
    assert mtwo / oo == zero
    assert mtwo / zoo == zero
    assert mtwo / undefined == undefined

    assert mtwo ** moo == zero
    assert mtwo ** oo == zoo
    assert mtwo ** zoo == undefined
    assert mtwo ** undefined == undefined

def test_monehalf_optable():
    assert monehalf + moo == moo
    assert monehalf + oo == oo
    assert monehalf + zoo == zoo
    assert monehalf + undefined == undefined

    assert monehalf - moo == oo
    assert monehalf - oo == moo
    assert monehalf - zoo == zoo
    assert monehalf - undefined == undefined

    assert monehalf * moo == oo
    assert monehalf * oo == moo
    assert monehalf * zoo == zoo
    assert monehalf * undefined == undefined

    assert monehalf / moo == zero
    assert mhalf / zero == zoo
    assert monehalf / oo == zero
    assert monehalf / zoo == zero
    assert monehalf / undefined == undefined

    assert monehalf ** moo == zero
    assert monehalf ** oo == zoo
    assert monehalf ** zoo == undefined
    assert monehalf ** undefined == undefined

def test_mone_optable():
    assert mone + moo == moo
    assert mone + oo == oo
    assert mone + undefined == undefined

    assert mone - moo == oo
    assert mone - oo == moo
    assert mone - undefined == undefined

    assert mone * moo == oo
    assert mone * oo == moo
    assert mone * undefined == undefined

    assert mone / moo == zero
    assert mhalf / zero == zoo
    assert mone / oo == zero
    assert mone / undefined == undefined

    assert mone ** moo == undefined
    assert mone ** oo == undefined
    assert mone ** undefined == undefined

def test_mhalf_optable():
    assert mhalf + moo == moo
    assert mhalf + oo == oo
    assert mhalf + zoo == zoo
    assert mhalf + undefined == undefined

    assert mhalf - moo == oo
    assert mhalf - oo == moo
    assert mhalf - zoo == zoo
    assert mhalf - undefined == undefined

    assert mhalf * moo == oo
    assert mhalf * oo == moo
    assert mhalf * zoo == zoo
    assert mhalf * undefined == undefined

    assert mhalf / moo == zero
    assert mhalf / zero == zoo
    assert mhalf / oo == zero
    assert mhalf / zoo == zero
    assert mhalf / undefined == undefined

    assert mhalf ** moo == zoo
    assert mhalf ** oo == zero
    assert mhalf ** zoo == undefined
    assert mhalf ** undefined == undefined

def test_zero_optable():
    assert zero + moo == moo
    assert zero + oo == oo
    assert zero + zoo == zoo
    assert zero + undefined == undefined

    assert zero - moo == oo
    assert zero - oo == moo
    assert zero - zoo == zoo
    assert zero - undefined == undefined

    assert zero * moo == undefined
    assert zero * oo == undefined
    assert zero * zoo == undefined
    assert zero * undefined == undefined

    assert zero / moo == zero
    assert zero / zero == undefined
    assert zero / oo == zero
    assert zero / zoo == zero
    assert zero / undefined == undefined

    assert zero ** moo == zoo
    assert zero ** oo == zero
    assert zero ** zoo == undefined
    assert zero ** undefined == undefined

def test_half_optable():
    assert half + moo == moo
    assert half + oo == oo
    assert half + zoo == zoo
    assert half + undefined == undefined

    assert half - moo == oo
    assert half - oo == moo
    assert half - zoo == zoo
    assert half - undefined == undefined

    assert half * moo == moo
    assert half * oo == oo
    assert half * zoo == zoo
    assert half * undefined == undefined

    assert half / moo == zero
    assert half / zero == zoo
    assert half / oo == zero
    assert half / zoo == zero
    assert half / undefined == undefined

    assert half ** moo == zoo
    assert half ** oo == zero
    assert half ** zoo == undefined
    assert half ** undefined == undefined

def test_fhalf_optable():
    assert fhalf + moo == moo
    assert fhalf + oo == oo
    assert fhalf + zoo == zoo
    assert fhalf + undefined == undefined

    assert fhalf - moo == oo
    assert fhalf - oo == moo
    assert fhalf - zoo == zoo
    assert fhalf - undefined == undefined

    assert fhalf * moo == moo
    assert fhalf * oo == oo
    assert fhalf * zoo == zoo
    assert fhalf * undefined == undefined

    assert fhalf / moo == zero
    assert fhalf / zero == zoo
    assert fhalf / oo == zero
    assert fhalf / zoo == zero
    assert fhalf / undefined == undefined

    assert fhalf ** moo == zoo
    assert fhalf ** oo == zero
    assert fhalf ** zoo == undefined
    assert fhalf ** undefined == undefined

def test_one_optable():
    assert one + moo == moo
    assert one + oo == oo
    assert one + zoo == zoo
    assert one + undefined == undefined

    assert one - moo == oo
    assert one - oo == moo
    assert one - zoo == zoo
    assert one - undefined == undefined

    assert one * moo == moo
    assert one * oo == oo
    assert one * zoo == zoo
    assert one * undefined == undefined

    assert one / moo == zero
    assert one / zero == zoo
    assert one / oo == zero
    assert one / zoo == zero
    assert one / undefined == undefined

    assert one ** moo == one
    assert one ** oo == one
    assert one ** zoo == one
    assert one ** undefined == one

def test_onehalf_optable():
    assert onehalf + moo == moo
    assert onehalf + oo == oo
    assert onehalf + zoo == zoo
    assert onehalf + undefined == undefined

    assert onehalf - moo == oo
    assert onehalf - oo == moo
    assert onehalf - zoo == zoo
    assert onehalf - undefined == undefined

    assert onehalf * moo == moo
    assert onehalf * oo == oo
    assert onehalf * zoo == zoo
    assert onehalf * undefined == undefined

    assert onehalf / moo == zero
    assert onehalf / zero == zoo
    assert onehalf / oo == zero
    assert onehalf / zoo == zero
    assert onehalf / undefined == undefined

    assert onehalf ** moo == zero
    assert onehalf ** oo == oo
    assert onehalf ** zoo == undefined
    assert onehalf ** undefined == undefined

def test_two_optable():
    assert two + moo == moo
    assert two + oo == oo
    assert two + zoo == zoo
    assert two + undefined == undefined

    assert two - moo == oo
    assert two - oo == moo
    assert two - zoo == zoo
    assert two - undefined == undefined

    assert two * moo == moo
    assert two * oo == oo
    assert two * zoo == zoo
    assert two * undefined == undefined

    assert two / moo == zero
    assert two / zero == zoo, `two/zero`
    assert two / oo == zero
    assert two / zoo == zero
    assert two / undefined == undefined

    assert two ** moo == zero
    assert two ** oo == oo
    assert two ** zoo == undefined
    assert two ** undefined == undefined

def test_oo_optable():
    assert oo + moo == undefined
    assert oo + mtwo == oo
    assert oo + monehalf == oo
    assert oo + mone == oo
    assert oo + mhalf == oo
    assert oo + zero == oo
    assert oo + half == oo
    assert oo + fhalf == oo
    assert oo + one == oo
    assert oo + onehalf == oo
    assert oo + two == oo
    assert oo + oo == oo
    assert oo + zoo == undefined
    assert oo + undefined == undefined

    assert oo - moo == oo
    assert oo - mtwo == oo
    assert oo - monehalf == oo
    assert oo - mone == oo
    assert oo - mhalf == oo
    assert oo - zero == oo
    assert oo - half == oo
    assert oo - fhalf == oo
    assert oo - one == oo
    assert oo - onehalf == oo
    assert oo - two == oo
    assert oo - oo == undefined
    assert oo - zoo == undefined
    assert oo - undefined == undefined

    assert oo * moo == moo
    assert oo * mtwo == moo
    assert oo * monehalf == moo
    assert oo * mone == moo
    assert oo * mhalf == moo
    assert oo * zero == undefined
    assert oo * half == oo
    assert oo * fhalf == oo
    assert oo * one == oo
    assert oo * onehalf == oo
    assert oo * two == oo
    assert oo * oo == oo
    assert oo * zoo == zoo
    assert oo * undefined == undefined

    assert oo / moo == undefined
    assert oo / mtwo == moo
    assert oo / monehalf == moo
    assert oo / mone == moo
    assert oo / mhalf == moo
    assert oo / zero == zoo
    assert oo / half == oo
    assert oo / fhalf == oo
    assert oo / one == oo
    assert oo / onehalf == oo
    assert oo / two == oo
    assert oo / oo == undefined
    assert oo / zoo == undefined
    assert oo / undefined == undefined

    assert oo ** moo == zero
    assert oo ** mtwo == zero
    assert oo ** monehalf == zero
    assert oo ** mone == zero
    assert oo ** mhalf == zero
    assert oo ** zero == one
    assert oo ** half == oo
    assert oo ** fhalf == oo
    assert oo ** one == oo
    assert oo ** onehalf == oo
    assert oo ** two == oo
    assert oo ** oo == oo
    assert oo ** zoo == undefined
    assert oo ** undefined == undefined

def test_zoo_optable():
    assert zoo + moo == undefined
    assert zoo + mtwo == zoo
    assert zoo + monehalf == zoo
    assert zoo + mone == zoo
    assert zoo + mhalf == zoo
    assert zoo + zero == zoo
    assert zoo + half == zoo
    assert zoo + fhalf == zoo
    assert zoo + one == zoo
    assert zoo + onehalf == zoo
    assert zoo + two == zoo
    assert zoo + oo == undefined
    assert zoo + zoo == undefined
    assert zoo + undefined == undefined

    assert zoo - moo == undefined
    assert zoo - mtwo == zoo
    assert zoo - monehalf == zoo
    assert zoo - mone == zoo
    assert zoo - mhalf == zoo
    assert zoo - zero == zoo
    assert zoo - half == zoo
    assert zoo - fhalf == zoo
    assert zoo - one == zoo
    assert zoo - onehalf == zoo
    assert zoo - two == zoo
    assert zoo - oo == undefined
    assert zoo - zoo == undefined
    assert zoo - undefined == undefined

    assert zoo * moo == zoo
    assert zoo * mtwo == zoo
    assert zoo * monehalf == zoo
    assert zoo * mone == zoo
    assert zoo * mhalf == zoo
    assert zoo * zero == undefined
    assert zoo * half == zoo
    assert zoo * fhalf == zoo
    assert zoo * one == zoo
    assert zoo * onehalf == zoo
    assert zoo * two == zoo
    assert zoo * oo == zoo
    assert zoo * zoo == zoo
    assert zoo * undefined == undefined

    assert zoo / moo == undefined
    assert zoo / mtwo == zoo
    assert zoo / monehalf == zoo
    assert zoo / mone == zoo
    assert zoo / mhalf == zoo
    assert zoo / zero == zoo
    assert zoo / half == zoo
    assert zoo / fhalf == zoo
    assert zoo / one == zoo
    assert zoo / onehalf == zoo
    assert zoo / two == zoo
    assert zoo / oo == undefined
    assert zoo / zoo == undefined
    assert zoo / undefined == undefined

    assert zoo ** moo == zero
    assert zoo ** mtwo == zero
    assert zoo ** monehalf == zero
    assert zoo ** mone == zero
    assert zoo ** mhalf == zero
    assert zoo ** zero == one
    assert zoo ** half == zoo
    assert zoo ** fhalf == zoo
    assert zoo ** one == zoo
    assert zoo ** onehalf == zoo
    assert zoo ** two == zoo
    assert zoo ** oo == zoo
    assert zoo ** zoo == undefined
    assert zoo ** undefined == undefined

def test_undefined_optable():
    assert undefined + moo == undefined
    assert undefined + mtwo == undefined
    assert undefined + monehalf == undefined
    assert undefined + mone == undefined
    assert undefined + mhalf == undefined
    assert undefined + zero == undefined
    assert undefined + half == undefined
    assert undefined + fhalf == undefined
    assert undefined + one == undefined
    assert undefined + onehalf == undefined
    assert undefined + two == undefined
    assert undefined + oo == undefined
    assert undefined + zoo == undefined
    assert undefined + undefined == undefined

    assert undefined - moo == undefined
    assert undefined - mtwo == undefined
    assert undefined - monehalf == undefined
    assert undefined - mone == undefined
    assert undefined - mhalf == undefined
    assert undefined - zero == undefined
    assert undefined - half == undefined
    assert undefined - fhalf == undefined
    assert undefined - one == undefined
    assert undefined - onehalf == undefined
    assert undefined - two == undefined
    assert undefined - oo == undefined
    assert undefined - zoo == undefined
    assert undefined - undefined == undefined

    assert undefined * moo == undefined
    assert undefined * mtwo == undefined
    assert undefined * monehalf == undefined
    assert undefined * mone == undefined
    assert undefined * mhalf == undefined
    assert undefined * zero == undefined
    assert undefined * half == undefined
    assert undefined * fhalf == undefined
    assert undefined * one == undefined
    assert undefined * onehalf == undefined
    assert undefined * two == undefined
    assert undefined * oo == undefined
    assert undefined * zoo == undefined
    assert undefined * undefined == undefined

    assert undefined / moo == undefined
    assert undefined / mtwo == undefined
    assert undefined / monehalf == undefined
    assert undefined / mone == undefined
    assert undefined / mhalf == undefined
    assert undefined / zero == undefined
    assert undefined / half == undefined
    assert undefined / fhalf == undefined
    assert undefined / one == undefined
    assert undefined / onehalf == undefined
    assert undefined / two == undefined
    assert undefined / oo == undefined
    assert undefined / zoo == undefined
    assert undefined / undefined == undefined

    assert undefined ** moo == undefined
    assert undefined ** mtwo == undefined
    assert undefined ** monehalf == undefined
    assert undefined ** mone == undefined
    assert undefined ** mhalf == undefined
    assert undefined ** zero == one
    assert undefined ** half == undefined
    assert undefined ** fhalf == undefined
    assert undefined ** one == undefined
    assert undefined ** onehalf == undefined
    assert undefined ** two == undefined
    assert undefined ** oo == undefined
    assert undefined ** zoo == undefined
    assert undefined ** undefined == undefined

def test_issue62():
    z = one / zoo
    assert z==0
    assert type(z) is Calculus, `type(z)`
