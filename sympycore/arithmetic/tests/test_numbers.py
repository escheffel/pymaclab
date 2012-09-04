from sympycore.arithmetic.numbers import *
from sympycore.arithmetic.infinity import *

mpq = normalized_fraction
mpc = mpqc

nan = undefined = Infinity(0)
oo = Infinity(1)
moo = Infinity(-1)
zoo = Infinity(undefined)

def test_Fraction():
    assert mpq(1) == 1
    assert mpq(1,1) == 1
    assert mpq(1,2) != 1
    assert mpq(1,2) == mpq(2,4)
    assert mpq(1,2) + mpq(5,6) == mpq(4,3)
    assert 2*mpq(1,2) == 1
    assert mpq(1,2)*2 == 1
    assert mpq(1,3) + mpq(2,3) == 1
    assert float(mpq(1,4)) == 0.25
    assert mpq(1,2) == mpq(-1,-2)
    assert mpq(-1,2) == mpq(1,-2)
    assert -mpq(2,3) == mpq(-2,3)
    assert +mpq(2,3) == mpq(2,3)
    assert mpq(1,2) - 1 == mpq(-1,2)
    assert 1 - mpq(1,2) == mpq(1,2)
    assert mpq(2,3)**0 == 1
    assert mpq(2,3)**1 == mpq(2,3)
    assert mpq(2,3)**2 == mpq(4,9)
    assert mpq(-2,3)**2 == mpq(4,9)
    assert mpq(-2,3)**-2 == mpq(9,4)
    assert mpq(-2,3)**3 == mpq(-8,27)
    assert mpq(-2,3)**-3 == mpq(-27,8)
    assert div(1,2) == mpq(1,2)
    assert div(3,mpq(1,2)) == 6
    assert div(mpq(1,2),mpq(3,2)) == mpq(1,3)
    assert mpq(1234,15) < 83
    assert mpq(1234,15) > 82
    assert mpq(2,3) < mpq(3, 4)
    assert mpq(28,3) // 4 == 2
    assert mpq(-28,3) // 4 == -3
    assert 4 // mpq(7,3) == 1
    assert 4 // mpq(-7,3) == -2
    assert mpq(5,4) % 1 == mpq(1,4)
    assert mpq(-5,4) % 1 == mpq(3,4)
    assert mpq(17,4) % mpq(6,7) == mpq(23,28)
    assert mpq(-17,4) % mpq(6,7) == mpq(1,28)

def test_mpf():
    # XXX
    setdps(15)
    assert mpf(2) != 3
    assert mpf(2) == 2
    assert mpf(1.1) == mpf('1.1')
    assert mpf(3) * mpf(4) == 12
    assert mpf(3) + mpf(4) == mpf(7)
    assert mpf(3) - mpf(4) == -1
    assert -mpf(3) == -3
    assert 4*mpf(3) == 12
    assert mpf(3)*4 == 12
    assert 2+mpf(5) == 7
    assert 0.5 + mpf(1) == 1.5
    assert mpf(1) + 0.5 == 1.5
    assert hash(mpf(3)) == hash(3)
    assert hash(mpf(1.5)) == hash(1.5)
    assert mpf(3) < mpf(4)
    assert mpf(3) > mpf(2)
    assert mpf(3) < 4

def test_mpc():
    assert mpc(2,3) == mpc(2,3)
    assert mpc(2,3) != mpc(2,4)
    assert mpc(3,0) == 3
    assert mpc(3,1) != 3
    assert mpc(2,3)*2 == mpc(4,6)
    assert 2*mpc(2,3) == mpc(4,6)
    assert mpc(2,3)+2 == mpc(4,3)
    assert 2+mpc(2,3) == mpc(4,3)
    assert mpc(2,3)-2 == mpc(0,3)
    assert 2-mpc(2,3) == mpc(0,-3)
    assert mpc(mpq(2,3), 1) + mpc(mpq(1,3), 1) == mpc(1, 2)
    assert -mpc(2,3) == mpc(-2,-3)
    assert +mpc(2,3) == mpc(2,3)
    assert mpc(0,1)**2 == -1
    assert mpc(0,1)**3 == mpc(0,-1)
    assert mpc(0,1)**4 == 1
    assert mpc(0,1)**(10**9) == 1
    assert mpc(0,1)**(10**9+1) == mpc(0,1)
    assert mpc(0,1)**(10**9+2) == -1
    assert mpc(0,1)**(10**9+3) == mpc(0,-1)
    assert mpc(3,4)**10 == mpc(-9653287,1476984)
    assert mpc(mpq(1,2), mpq(-5,6))**3 == mpc(mpq(-11,12), mpq(-5,108))
    assert str(mpc(4,-3)) == '4 - 3*I'
    assert str(mpc(0,mpq(1,2))) == '1/2*I'
    assert mpc(0,1)**(-2) == -1
    assert abs(mpc(0,2)) == 2
    assert abs(mpc(mpf(3.0),mpf(4.0))) == 5.0

def test_extended_numbers():
    assert 1*oo == oo
    assert -3*oo == -oo
    assert oo*oo == oo
    assert oo+oo == oo
    assert mpq(1,3)*oo == oo
    assert -oo == (-1)*oo
    assert -(-oo) == oo
    assert oo - oo == nan
    assert -oo + -oo == -oo
    assert -oo - (-oo) == nan
    assert oo + nan == nan
    assert (mpc(0,1)*oo)*mpc(0,1) == -oo
    assert zoo * oo == zoo
    assert zoo * -oo == zoo
    assert zoo * nan == nan
    assert zoo + nan == nan
    assert 0*oo == nan
    assert 0*zoo == nan
    assert 0*nan == nan
    assert abs(oo) == abs(-oo) == abs(zoo) == oo
    assert abs(nan) == nan

def test_extended_cmp():
    assert -oo < oo
    assert -oo < 3
    assert -oo < mpq(-5,4)
    assert -oo <= 3
    #assert -oo <= -oo
    assert not (-oo > 3)
    assert not (-oo > -oo)
    assert not (-oo < -oo)
    assert oo > 3
    assert oo > -oo
    #assert oo <= oo
    assert not (oo < 3)
    assert not (oo < oo)
    assert not (oo > oo)
    assert not (nan < nan)
    assert not (nan > nan)
    assert max(2, 3, 1, nan, 2) == 3
    assert min(2, 3, 1, nan, 2) == 1
    assert max(2, -oo, oo, 3) == oo
    assert min(2, -oo, oo, 3) == -oo

def test_int_roots():
    assert int_root(1,1) == (1, True)
    assert int_root(0,1) == (0, True)
    assert int_root(0,3) == (0, True)
    assert int_root(10000, 1) == (10000, True)
    assert int_root(4,2) == (2, True)
    assert int_root(16,2) == (4, True)
    assert int_root(26,2) == (5, False)
    assert int_root(1234567**7, 7) == (1234567, True)
    assert int_root(1234567**7+1, 7) == (1234567, False)
    assert int_root(1234567**7-1, 7) == (1234566, False)
    b = 25**1000
    assert int_root(b, 1000) == (25, True)
    assert int_root(b+1, 1000) == (25, False)
    assert int_root(b-1, 1000) == (24, False)
    c = 10**400
    c2 = c**2
    assert int_root(c2, 2) == (c, True)
    assert int_root(c2+1, 2) == (c, False)
    assert int_root(c2-1, 2) == (c-1, False)
    assert int_root(2,10**10) == (1, False)

def test_powers():
    assert try_power(3, 2) == (9, [])
    assert try_power(3, -2) == (mpq(1, 9), [])
    assert try_power(0, -1) == (zoo, [])
    assert try_power(0, 0) == (1, [])
    assert try_power(mpq(1,2), 0) == (1, [])
    assert try_power(mpc(0, 1), 2) == (-1, [])
    assert try_power(mpc(-1,2), -2) == (mpc(mpq(-3,25), mpq(4,25)), [])
    assert try_power(2, mpq(1, 2)) == (1, [(2, mpq(1, 2))])
    assert try_power(4, mpq(1, 2)) == (2, [])
    assert try_power(4, mpq(3, 2)) == (8, [])
    assert try_power(4, mpq(-3, 2)) == (mpq(1,8), [])
    assert try_power(-1, mpq(1, 2)) == (mpc(0, 1), [])
    assert try_power(-729, mpq(1, 6)) == (3, [(-1, mpq(1, 6))])
    assert try_power(mpq(4,9), mpq(1,2)) == (mpq(2, 3), [])
    assert try_power(mpq(-4,9), mpq(1,2)) == (mpc(0, mpq(2, 3)), [])
    assert try_power(mpq(-4,9), mpq(-1,2)) == (mpc(0, mpq(-3, 2)), [])
    assert try_power(mpq(4,7), mpq(1,2)) == (2, [(7, mpq(-1, 2))])
    assert try_power(mpq(7,4), mpq(1,2)) == (mpq(1, 2), [(7, mpq(1, 2))])
    assert try_power(oo, 1) == (oo, [])
    assert try_power(oo, 3) == (oo, [])
    assert try_power(oo, 0) == (1, [])
    assert try_power(oo, -1) == (0, [])
    assert try_power(oo, -3) == (0, [])
    assert try_power(oo, mpq(-2,3)) == (0, [])
    assert try_power(-oo, 1) == (-oo, [])
    assert try_power(-oo, 2) == (oo, [])
    assert try_power(-oo, 3) == (-oo, [])
    assert try_power(-oo, mpq(5,2)) == (oo*mpc(0,1), [])
    #assert try_power(oo*mpc(0,1), mpq(5,2)) == (oo, [(mpc(0, 1), mpq(5, 2))])
    assert try_power(zoo, 2) == (zoo, [])
    assert try_power(-zoo, 3) == (zoo, [])
    assert try_power(zoo, -1) == (0, [])
