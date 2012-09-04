from sympycore.arithmetic.numbers import normalized_fraction as mpq
from sympycore.arithmetic.numbers import mpf
from sympycore.arithmetic.number_theory import *
from sympycore.arithmetic.number_theory import binomial_coefficients, binomial_coefficients_list

def test_factorial():
    assert factorial(0) == 1
    assert factorial(1) == 1
    assert factorial(2) == 2
    assert factorial(10) == 3628800
    # Verify that memoization didn't cause trouble
    assert factorial(10) == 3628800
    assert factorial(130) % 1009 == 926

def test_gcd():
    assert gcd() == 0
    assert gcd(2) == 2
    assert gcd(2,3) == 1
    assert gcd(2*3*5, 2*3*7) == 6
    assert gcd(2*3*5, 2*5, 3*5) == 5
    assert gcd(mpq(2,3),mpq(5,8)) == mpq(1,24)

def test_lcm():
    assert lcm(0) == 0
    assert lcm(1) == 1
    assert lcm(12,18) == 36
    assert lcm(7,11,7) == 77
    assert lcm(12,16,15,9) == 720
    assert lcm(mpq(1,2),mpq(5,6)) == mpq(5,2)

def test_integer_digits():
    assert integer_digits(0) == [0]
    assert integer_digits(1234, 10) == [1,2,3,4]
    assert integer_digits(37, 2) == [1,0,0,1,0,1]

def test_real_digits():
    assert real_digits(0) == ([0], 1)
    assert real_digits(28) == ([2,8], 2)
    assert real_digits(mpq(3,2)) == ([1,5], 1)
    assert real_digits(mpq(1,1000)) == ([1], -2)
    assert real_digits(mpq(355,113)) == ([3,1,4,1,5,9,2,9,2,0], 1)
    assert real_digits(mpq(1,10),2) == ([1,1,0,0,1,1,0,0,1,1], -3)
    # 1234.1234...
    p = mpq(1234123456789,1000000000)
    assert real_digits(p, 10, 8) == ([1,2,3,4,1,2,3,4],4)
    assert real_digits(p, 10, 12) == ([1,2,3,4,1,2,3,4,5,6,7,8], 4)

def test_binomial_coefficients_list():
    assert binomial_coefficients_list(0) == [1]
    assert binomial_coefficients_list(1) == [1,1]
    assert binomial_coefficients_list(2) == [1,2,1]
    assert binomial_coefficients_list(3) == [1,3,3,1]
    assert binomial_coefficients_list(4) == [1,4,6,4,1]
    assert binomial_coefficients_list(5) == [1,5,10,10,5,1]
    assert binomial_coefficients_list(6) == [1,6,15,20,15,6,1]

def test_binomial_coefficients():
    for n in range(15):
        c = binomial_coefficients(n)
        l = [c[k] for k in sorted(c)]
        assert l==binomial_coefficients_list(n)

def test_f2q():
    assert f2q(mpf(20)/521)==mpq(20, 521)
    assert f2q(mpf(1)/3)==mpq(1, 3), `f2q(mpf(1)/3)`
    assert f2q(mpf(1232)/5224)==mpq(154, 653), `f2q(mpf(1232)/5224)`
