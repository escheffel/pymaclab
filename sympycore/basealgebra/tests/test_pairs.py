
from sympycore import CollectingField as Algebra

from sympycore.heads import *

Symbol = Algebra.Symbol
Number = Algebra.Number
Add = Algebra.Add
Sub = Algebra.Sub
Mul = Algebra.Mul
Pow = Algebra.Pow
Terms = Algebra.Terms
Factors = Algebra.Factors
Apply = Algebra.Apply

def test_symbol():
    a = Symbol('a')
    assert str(a)=='a'

def test_number():
    n = Number(2)
    assert str(n)=='2'

def test_add():
    a = Symbol('a')
    n = Number(2)
    s = Add(n, a, Number(3), 0, Number(0), -3)
    assert str(s) in ['2 + a','a + 2'], str(s)
    s += n
    assert str(s) in ['4 + a', 'a + 4'], `str(s),repr(s)`
    s += 3
    assert str(s) in ['7 + a', 'a + 7']
    s += 0
    assert str(s) in ['7 + a', 'a + 7']
    s += Number(0)
    assert str(s) in ['7 + a', 'a + 7']
    s += -7
    assert str(s)=='a'
    assert s.head is SYMBOL
    s += 2
    assert str(s) in ['2 + a', 'a + 2']
    s += Number(-2)
    assert str(s)=='a'
    assert s.head is SYMBOL
        
    assert str(Add('x+y','-x'))=='y',str(Add('x+y','-x').pair)

    s = Sub(a,n, Number(4), 0, Number(0), -3)
    assert str(s) in ['a - 3', '-3 + a'], str(s)
    s -= n
    assert str(s) in ['a - 5', '-5 + a'], str(s)
    s -= -7
    assert str(s) in ['2 + a', 'a + 2']
    s -= 0
    assert str(s) in ['2 + a', 'a + 2']
    s -= Number(0)
    assert str(s) in ['2 + a', 'a + 2']
    s -= 2
    assert str(s)=='a'
    assert s.head is SYMBOL
    s -= -6
    assert str(s) in ['6 + a', 'a + 6']
    s -= Number(6)
    assert str(s)=='a'
    
def test_mul():
    a = Symbol('a')
    n = Number(2)
    s = Mul(n,a)
    assert str(s)=='2*a'

    assert str((2/a) / (2*a)) in ['a**(-2)', 'a**-2', '1/a**2'],str((2/a) / (2*a))

def test_pow():
    a = Symbol('a')
    n = Number(2)
    s = Pow(a,3)
    assert str(s)=='a**3'

def test_terms():
    a = Symbol('a')
    n = Number(2)
    s = Terms((a,2))
    assert str(s)=='2*a'

def test_factors():
    a = Symbol('a')
    n = Number(2)
    s = Factors((a,2))
    assert str(s)=='a**2'

def test_new():
    a = Algebra.convert('a')
    assert str(a)=='a'
    assert isinstance(a, Algebra)==True
    assert (a is Algebra.convert(a))==True

def _test_copy():
    a = Symbol('a')
    assert (a.copy() is a)==True

def _test_func_args():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    s = a + n
    s1 = 2*a
    m = a ** 2
    m2 = a*b
    assert a.func(*a.args)==a
    assert n.func(*n.args)==n
    assert s.func(*s.args)==s
    assert s1.func(*s1.args)==s1
    assert m.func(*m.args)==m
    assert m2.func(*m2.args)==m2

def _test_apply_args():
    a = Symbol('a')
    b = Symbol('b')
    f = Symbol('f')
    fab = Apply(f, a,b)
    assert fab.args==(a,b)
    assert fab.func(*fab.args)==fab
    
def _test_Add_args():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    s = a + n
    s1 = 2*a
    m = a ** 2
    m2 = a*b

    assert Add(*a.as_Add_args())==a
    assert Add(*n.as_Add_args())==n
    assert Add(*s.as_Add_args())==s
    assert Add(*s1.as_Add_args())==s1
    assert Add(*m.as_Add_args())==m
    assert Add(*m2.as_Add_args())==m2

def _test_Mul_args():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    s = a + n
    s1 = 2*a
    m = a ** 2
    m2 = a*b

    assert Mul(*a.as_Mul_args())==a
    assert Mul(*n.as_Mul_args())==n
    assert Mul(*s.as_Mul_args())==s
    assert Mul(*s1.as_Mul_args())==s1
    assert Mul(*m.as_Mul_args())==m
    assert Mul(*m2.as_Mul_args())==m2

def _test_Pow_args():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    s = a + n
    s1 = 2*a
    m = a ** 2
    m2 = a*b

    assert Pow(*a.as_Pow_args())==a
    assert Pow(*n.as_Pow_args())==n
    assert Pow(*s.as_Pow_args())==s
    assert Pow(*s1.as_Pow_args())==s1
    assert Pow(*m.as_Pow_args())==m
    assert Pow(*m2.as_Pow_args())==m2

def _test_Terms_args():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    s = a + n
    s1 = 2*a
    m = a ** 2
    m2 = a*b

    assert Terms(*a.as_Terms_args())==a
    assert Terms(*n.as_Terms_args())==n
    assert Terms(*s.as_Terms_args())==s
    assert Terms(*s1.as_Terms_args())==s1
    assert Terms(*m.as_Terms_args())==m
    assert Terms(*m2.as_Terms_args())==m2

def _test_Factors_args():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    s = a + n
    s1 = 2*a
    m = a ** 2
    m2 = a*b

    assert Factors(*a.as_Factors_args())==a
    assert Factors(*n.as_Factors_args())==n
    assert Factors(*s.as_Factors_args())==s
    assert Factors(*s1.as_Factors_args())==s1
    assert Factors(*m.as_Factors_args())==m
    assert Factors(*m2.as_Factors_args())==m2

def test_as_verbatim():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    t = Number(-3)
    s = a + n
    s1 = 2*a
    m = a ** 2
    m2 = a*b

    assert str(a.as_verbatim())=='a'
    assert str(n.as_verbatim())=='2'
    assert str(t.as_verbatim())=='-3', `str(t.as_verbatim())`
    assert str(s.as_verbatim()) in ['2 + a','a + 2', 'a*1 + 1*2', '1*2 + a*1'],`str(s.as_verbatim())`
    assert str(s1.as_verbatim()) in ['2*a','a*2'], str(s1.as_verbatim())
    assert str(m.as_verbatim())=='a**2',` str(m.as_verbatim())`
    assert str(m2.as_verbatim()) in ['a*b','b*a'], str(m2.as_verbatim())

def test_Mul():
    a = Symbol('a')
    b = Symbol('b')
    n = Number(2)
    assert Mul(a,a**-1)==1
    assert Mul(n,n**-1)==1

def test_neg():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert str(-n)==str('-3')
    assert str(-a)==str('-a'),str(-a)
    assert str(-s) in [str('-2 - a'), '-a - 2'], str(-s)
    assert str(-s1)==str('-2*a')
    assert str(-m) in ['-a*b','-(a*b)','-(b*a)'],`str(-m)`

def test_number_add():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert n+2==5
    assert 2+n==5
    assert n-1==2
    assert 1-n==-2
    assert n+n==6
    assert n-n==0
    assert str(n+a) in [str('3 + a'), 'a + 3']
    assert str(n+s) in [str('5 + a'), 'a + 5']
    assert str((-2)+s)==str('a')
    assert str(s+(-2))==str('a')
    assert str(n+s1) in [str('3 + 2*a'), '2*a + 3']
    assert str(n+m) in ['3 + a*b', '3 + b*a'], str(n+m)
    assert str(n+m1) in [str('3 + a**2'), 'a**2 + 3']
    assert str(n-a) in [str('3 - a'), '-a + 3']
    assert str(n-s) in [str('1 - a'), '-a + 1']
    assert str(n-s1) in [str('3 - 2*a'), '-2*a + 3'], str(n-s1)
    assert str(n-m) in ['3 - a*b', '3 - (a*b)', '3 - (b*a)'],`str(n-m)`
    assert str(n-m1)==str('3 - a**2')
    assert str(s - 0)==str(s)
    assert str(s - Number(0))==str(s)
    assert str(s - 2)==str('a')
    assert str(s - Number(2))==str('a')
    assert str(2 - (2 - a))==str('a')
    assert str(Number(2) - (2 - a))==str('a')

def test_number_mul():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert n*2==6
    assert n*n==9
    assert str(n*a)==str('3*a')
    assert str(n*s) in [str('6 + 3*a'), '3*a + 6']
    assert str(n*s1)==str('6*a')
    assert str(n*m) in ['3*a*b','3*b*a'],str(n*m)
    assert str(n*m1)==str('3*a**2')
    assert str(m1*0)==str('0')
    assert str(m1*Number(0))==str('0')

def test_number_pow():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert n**1==3
    assert n**0==1
    assert n**2==9
    assert n**n==27
    assert str(n**a)==str('3**a')
    assert str(n**s) in [str('3**(2 + a)'),'3**(a + 2)'], str(n**s)
    assert str(n**s1)==str('3**(2*a)')
    assert str(n**m) in ['3**(a*b)','3**(b*a)'], str(n**m)
    assert str(n**m1) in ['3**(a**2)',
                          '3**a**2'], str(n**m1)

def test_symbol_add():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert str(a+2) in [str('2 + a'), 'a + 2']
    assert str(a+a)==str('2*a')
    assert str(a+n) in [str('3 + a'), 'a + 3']
    assert str(a+s) in [str('2 + 2*a'), '2*a + 2']
    assert str(a+s1)==str('3*a')
    assert str(b+s1) in [str('b + 2*a'),'2*a + b'], str(b+s1)
    assert str(a+m) in [str('a + a*b'), 'a*b + a','a + b*a'], str(a+m)
    assert str(a+m1) in [str('a + a**2'), 'a**2 + a'], str(a+m1)

def test_symbol_mul():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert str(a*2)=='2*a'
    assert str(a*a)=='a**2'
    assert str(a*n)==str('3*a')
    assert str(a*s) in [str('a*(2 + a)'), 'a*(a + 2)']
    assert str(a*s1)==str('2*a**2')
    assert str(b*s1) in ['2*a*b','2*b*a'],str(b*s1)
    assert str(a*m) in [str('b*a**2'),'a**2*b']
    assert str(a*m1)==str('a**3')
    assert str(b*m1) in [str('b*a**2'),'a**2*b']

def test_symbol_pow():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert a**1==a
    assert a**0==1
    assert str(a**2)=='a**2'
    assert str(a**n)=='a**3'
    assert str(a**a)==str('a**a')
    assert str(a**s) in [str('a**(2 + a)'), 'a**(a + 2)']
    assert str(a**s1)==str('a**(2*a)')
    assert str(a**m) in ['a**(a*b)','a**(b*a)'],str(a**m)
    assert str(a**m1) in [str('a**(a**2)'), 'a**a**2']

def test_add_add():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    s2 = a/2
    s3 = 2/a
    m = a*b
    m1 = a**2
    assert str(s+2) in [str('4 + a'), 'a + 4']
    assert str(s+a) in [str('2 + 2*a'), '2*a + 2']
    assert str(s+n) in [str('5 + a'), 'a + 5']
    assert str(s+s) in [str('4 + 2*a'), '2*a + 4']
    assert str(s+s1) in [str('2 + 3*a'), '3*a + 2']
    assert str(s+m) in [str('2 + a + a*b'), '2 + a*b + a', 'a + 2 + a*b','2 + a + b*a'], str(s+m)
    assert str(s+m1) in [str('2 + a + a**2'), '2 + a**2 + a', 'a + 2 + a**2'], str(s+m1)

    assert str(s1+2) in [str('2 + 2*a'), '2*a + 2']
    assert str(s1+a)==str('3*a')
    assert str(s1+n) in [str('3 + 2*a'), '2*a + 3']
    assert str(s1+s) in [str('2 + 3*a'), '3*a + 2']
    assert str(s1+s1)==str('4*a')
    assert str(s1+m) in [str('2*a + a*b'), 'a*b + 2*a','2*a + b*a'],str(s1+m)
    assert str(s1+m1) in ['2*a + a**2','a**2 + 2*a']

    assert str(s-2)==str('a')
    assert str(s-a)==str('2')
    assert str(s-n) in [str('a - 1'),'-1 + a']
    assert str(s-s)==str('0')
    assert str(s-s1) in [str('2 - a'), '-a + 2']
    assert str(s-m) in [str('2 + a - a*b'), '2 - a*b + a', 'a + 2 - a*b',
                        'a + 2 - (a*b)', '2 + a - (b*a)'], str(s-m)
    assert str(s-m1) in [str('2 + a - a**2'), '2 - a**2 + a', 'a + 2 - a**2'], str(s-m1)

    assert str(s1-2) in [str('2*a - 2'),'-2 + 2*a'], str(s1-2)
    assert str(s1-a)==str('a')
    assert str(s1-n) in [str('2*a - 3'), '-3 + 2*a'], str(s1-n)
    assert str(s1-s) in [str('a - 2'), '-2 + a'], str(s1-s)
    assert str(s1-s1)==str('0')
    assert str(s1-m) in [str('2*a - a*b'), '-a*b + 2*a', '2*a - (a*b)',
                         '2*a - (b*a)'], str(s1-m)
    assert str(s1-m1) in [str('2*a - a**2'), '-a**2 + 2*a'], str(s1-m1)

    assert str((2*a)*(a/2))==str('a**2')
    assert str((2*a)*(2/a))==str('4')
    assert str((a+b) * (2/(a+b))) == str('2')

def test_add_mul():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert str(s*2) in ['4 + 2*a','2*a + 4']
    assert str(s*a) in ['a*(2 + a)', 'a*(a + 2)']
    assert str(s*n) in [str('6 + 3*a'), '3*a + 6']
    assert str(s*s) in [str('(2 + a)**2'), '(a + 2)**2'], str(s*s)
    assert str(s*s1) in [str('2*a*(2 + a)'), '2*a*(a + 2)']
    assert str(s*m) in [str('a*b*(2 + a)'), 'a*(2 + a)*b', 'a*(a + 2)*b',
                        'b*a*(2 + a)', 'a*b*(a + 2)'], str(s*m)
    assert str(s*m1) in ['(2 + a)*a**2','a**2*(2 + a)', 'a**2*(a + 2)'], str(s*m1)

    assert str(s1*2)=='4*a'
    assert str(s1*a)=='2*a**2'
    assert str(s1*n)==str('6*a')
    assert str(s1*s) in [str('2*a*(2 + a)'), '2*a*(a + 2)'], str(s1*s)
    assert str(s1*s1)==str('4*a**2')
    assert str(s1*m) in [str('2*b*a**2'), '2*a**2*b'],str(s1*m)
    assert str(s1*m1)==str('2*a**3')

    assert str(s*0)==str('0')
    assert str(s*Number(0))==str('0')

def test_add_pow():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    
    assert s**1==s
    assert s**0==1
    assert str(s**2) in ['(2 + a)**2', '(a + 2)**2']
    assert str(s**-2) in ['(2 + a)**(-2)', '(2 + a)**-2', '(a + 2)**(-2)',
                          '(a + 2)**-2', '1/(a + 2)**2', '1/(2 + a)**2'], str(s**-2)
    assert str(s**n) in ['(2 + a)**3','(a + 2)**3']
    assert str(s**a) in [str('(2 + a)**a'), '(a + 2)**a']
    assert str(s**s) in [str('(2 + a)**(2 + a)'), '(a + 2)**(a + 2)']
    assert str(s**s1) in [str('(2 + a)**(2*a)'), '(a + 2)**(2*a)'], str(s**s1)
    assert str(s**m) in [str('(2 + a)**(a*b)'), '(a + 2)**(a*b)',
                         '(2 + a)**(b*a)'], str(s**m)
    assert str(s**m1) in ['(2 + a)**(a**2)','(2 + a)**a**2', '(a + 2)**a**2',
                          '(a + 2)**(a**2)'],str(s**m1)

    assert s1**1==s1
    assert s1**0==1
    assert str(s1**2)=='4*a**2'
    assert str(s1**-2) in ['1/4*a**(-2)', '1/4*a**-2','1/4/a**2'],str(s1**-2)
    assert str(s1**n)=='8*a**3'
    assert str(s1**a)==str('(2*a)**a')
    assert str(s1**s) in [str('(2*a)**(2 + a)'),'(2*a)**(a + 2)'], str(s1**s)
    assert str(s1**s1)==str('(2*a)**(2*a)')
    assert str(s1**m) in ['(2*a)**(a*b)','(2*a)**(b*a)'],str(s1**m)
    assert str(s1**m1) in [str('(2*a)**(a**2)'),'(2*a)**a**2'], str(s1**m1)

def test_mul_add():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert str(m+2) in ['2 + a*b','2 + b*a', 'a*b + 2'],str(m+2)
    assert str(m+a) in [str('a + a*b'), 'a*b + a','a + b*a'], str(m+a)
    assert str(m+n) in ['3 + a*b','3 + b*a', 'a*b + 3'],  str(m+n)
    assert str(m+s) in [str('2 + a + a*b'),'2 + a*b + a', 'a + 2 + a*b',
                        '2 + a + b*a', 'a*b + a + 2'],str(m+s)
    assert str(m+s1) in [str('2*a + a*b'),'a*b + 2*a', '2*a + b*a'],str(m+s1)
    assert str(m+m) in ['2*a*b','2*b*a'],str(m+m)
    assert str(m+m1) in ['a*b + a**2','a**2 + a*b','b*a + a**2'], str(m+m1)

    assert str(m1+2) in [str('2 + a**2'), 'a**2 + 2'], str(m1+2)
    assert str(m1+a) in [str('a + a**2'), 'a**2 + a'],str(m1+a)
    assert str(m1+n) in [str('3 + a**2'), 'a**2 + 3'],str(m1+n)
    assert str(m1+s) in [str('2 + a + a**2'), '2 + a**2 + a', 'a + 2 + a**2', 'a**2 + a + 2'],str(m1+s)
    assert str(m1+s1) in ['2*a + a**2','a**2 + 2*a']
    assert str(m1+m) in ['a*b + a**2','a**2 + a*b', 'b*a + a**2'],str(m1+m)
    assert str(m1+m1)==str('2*a**2')

def test_mul_mul():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    c = (-2*a*b)**(1/Number(2))
    d = (2*a*b)**(1/Number(2))
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2
    assert str(m*2) in ['2*a*b', '2*b*a'],str(m*2)
    assert str(m*a) in ['b*a**2','a**2*b'], str(m*a)
    assert str(m*n) in ['3*a*b','3*b*a'], str(m*n)
    assert str(m*s) in [str('a*b*(2 + a)'),'a*(2 + a)*b', 'a*(a + 2)*b',
                        'b*a*(2 + a)', 'a*b*(a + 2)'],str(m*s)
    assert str(m*s1) in [str('2*b*a**2'),'2*a**2*b'],str(m*s1)
    assert str(m*m) in ['a**2*b**2','b**2*a**2'],str(m*m)
    assert str(m*m1) in [str('b*a**3'),'a**3*b'],str(m*m1)
    assert str(c*c) in ['-2*a*b', '-2*b*a'],str(c*c)
    assert str(d*d) in ['2*a*b','2*b*a'],str(d*d)

    assert str(m1*2)=='2*a**2'
    assert str(m1*a)=='a**3'
    assert str(m1*n)==str('3*a**2')
    assert str(m1*s) in ['(2 + a)*a**2','a**2*(2 + a)', 'a**2*(a + 2)'], str(m1*s)
    assert str(m1*s1)==str('2*a**3')
    assert str(m1*m) in [str('b*a**3'),'a**3*b'],str(m1*m)
    assert str(m1*m1)==str('a**4')

    assert str((a**2)*(a**-2))==str('1')
    assert str(((a*b)**2)*(a**-2))==str('b**2')
    assert str((a**b)*(a**-b))==str('1')
    assert str((a**b*b**b)*(a**-b))==str('b**b')
    assert str(((2*a)**b)*((2*a)**-b))==str('1')
    assert str((a)**(1/n) * (a)**(2/n))==str('a')
    assert str((2*a)**(1/n) * (2*a)**(2/n))==str('2*a')
    assert str(((a)**(1/n)*2**(2/n)) * ((a)**(2/n) * 2**(1/n))) in [str('2*a'),'a*2'], str(((a)**(1/n)*2**(2/n)) * ((a)**(2/n) * 2**(1/n)))
    assert str(((2*a)**(1/n)*2**(2/n)) * ((2*a)**(2/n) * 2**(1/n)))==str('4*a'), str(((2*a)**(1/n)*2**(2/n)) * ((2*a)**(2/n) * 2**(1/n)))

def test_mul_pow():
    n = Number(3)
    a = Symbol('a')
    b = Symbol('b')
    s = 2 + a
    s1 = 2*a
    m = a*b
    m1 = a**2

    assert m**1==m
    assert m**0==1
    assert str(m**2) in ['a**2*b**2', 'b**2*a**2'],str(m**2)
    assert str(m**n) in ['a**3*b**3','b**3*a**3'],  str(m**n)
    assert str(m**a) in ['(a*b)**a', '(b*a)**a'], str(m**a)
    assert str(m**s) in [str('(a*b)**(2 + a)'), '(a*b)**(a + 2)',
                         '(b*a)**(2 + a)'],str(m**s)
    assert str(m**s1) in ['(a*b)**(2*a)','(b*a)**(2*a)'], str(m**s1)
    assert str(m**m) in ['(a*b)**(a*b)','(b*a)**(b*a)'], str(m**m)
    assert str(m**m1) in [str('(a*b)**(a**2)'),'(a*b)**a**2',
                          '(b*a)**a**2'],str(m**m1)

    assert str(m1**2)=='a**4'
    assert str(m1**n)=='a**6'
    assert str(m1**a)==str('(a**2)**a'), `str(m1**a), (m1**a)`
    assert str(m1**s) in [str('(a**2)**(2 + a)'), '(a**2)**(a + 2)'], str(m1**s)
    assert str(m1**s1)==str('(a**2)**(2*a)')
    assert str(m1**m) in ['(a**2)**(a*b)', '(a**2)**(b*a)'],str(m1**m)
    assert str(m1**m1) in [str('(a**2)**(a**2)'), '(a**2)**a**2'], str(m1**m1)
    assert str((a**(1/n))**n)==str('a')

def test_expand():
    x,y,z = map(Symbol, 'xyz')
    assert ((x+y)**2).expand()==x**2+y**2+2*x*y
    assert str(((x+y)**2).expand()) in ['2*x*y + x**2 + y**2',
                                        'x**2 + y**2 + 2*x*y',
                                        'x**2 + y**2 + 2*y*x',
                                        'x**2 + 2*x*y + y**2',
                                        '2*y*x + x**2 + y**2'], str(((x+y)**2).expand())
    assert ((x-y)**2).expand()==x**2+y**2-2*x*y, `((x-y)**2).expand().pair, (x**2+y**2-2*x*y).pair`
    assert ((x+y)**3).expand()==x**3+y**3+3*x**2*y+3*x*y**2, ((x+y)**3).expand()
    assert ((x-y)**3).expand()==x**3-y**3-3*x**2*y+3*x*y**2
    assert ((x+y)**3).expand()==-((-x-y)**3).expand()
    assert str((x*(x+y)).expand()) in ['x*y + x**2','x**2 + x*y', 'x**2 + y*x',
                                       'y*x + x**2'],str((x*(x+y)).expand())
    assert str(((x+y)*x).expand()) in ['x*y + x**2','x**2 + x*y', 'x**2 + y*x',
                                       'y*x + x**2'],str(((x+y)*x).expand())
    
    assert str(((x+y+z)**2).expand()) in ['2*x*y + 2*x*z + 2*y*z + x**2 + y**2 + z**2',
                                          'x**2 + y**2 + z**2 + 2*x*y + 2*x*z + 2*y*z',
                                          'x**2 + y**2 + z**2 + 2*y*x + 2*x*z + 2*y*z',
                                          '2*y*z + z**2 + x**2 + 2*y*x + y**2 + 2*x*z',
                                          'z**2 + 2*x*z + 2*x*y + 2*z*y + y**2 + x**2',
                                          '2*y*z + z**2 + y**2 + 2*y*x + x**2 + 2*x*z',
                                          '2*y*z + 2*x*z + y**2 + z**2 + 2*y*x + x**2'],  str(((x+y+z)**2).expand())
    r = ' + '.join(sorted(str(((x+y+z)**3).expand()).split(' + ')))
    assert r in ['3*x*y**2 + 3*x*z**2 + 3*y*x**2 + 3*y*z**2 + 3*z*x**2 + 3*z*y**2 + 6*x*y*z + x**3 + y**3 + z**3',
                 '3*x**2*z + 3*x*z**2 + 3*y**2*x + 3*y**2*z + 3*y*x**2 + 3*y*z**2 + 6*y*x*z + x**3 + y**3 + z**3',
                 '3*x**2*y + 3*x**2*z + 3*x*y**2 + 3*x*z**2 + 3*z**2*y + 3*z*y**2 + 6*x*z*y + x**3 + y**3 + z**3'],\
                 r

    assert str(((2*x+y)**2).expand()) in ['4*x**2 + 4*x*y + y**2',
                                          'y**2 + 4*x**2 + 4*x*y',
                                          '4*x**2 + y**2 + 4*y*x',
                                          '4*y*x + 4*x**2 + y**2'], str(((2*x+y)**2).expand())
    assert str(((2*x-y)**2).expand()) in ['4*x**2 + y**2 - 4*x*y',
                                          'y**2 + 4*x**2 - 4*x*y',
                                          '4*x**2 + y**2 - 4*y*x',
                                          '4*x**2 - 4*x*y + y**2',
                                          '-4*y*x + 4*x**2 + y**2'],\
                                          str(((2*x-y)**2).expand())
    assert str(((x+y)**2-x**2-y**2).expand()) in ['2*x*y','2*y*x'],  str(((x+y)**2-x**2-y**2).expand())
    assert str((((x+y)**2-x**2-y**2)*(x*y)).expand()) in ['2*x**2*y**2','2*y**2*x**2'], str((((x+y)**2-x**2-y**2)*(x*y)).expand())
    assert str(((x*y)*((x+y)**2-x**2-y**2)).expand()) in ['2*x**2*y**2','2*y**2*x**2'], str(((x*y)*((x+y)**2-x**2-y**2)).expand())
    assert str(((3*x*y)*((x+y)**2-x**2-y**2)).expand()) in ['6*x**2*y**2','6*y**2*x**2'], str(((3*x*y)*((x+y)**2-x**2-y**2)).expand())
    assert str(((1/x)*((x+y)**2-x**2-y**2)).expand())=='2*y', str(((1/x)*((x+y)**2-x**2-y**2)).expand())
    assert str(((3/x)*((x+y)**2-x**2-y**2)).expand())=='6*y', str(((3/x)*((x+y)**2-x**2-y**2)).expand())
    assert str(((x**2)*((x+y)**2-x**2-y**2)).expand()) in ['2*y*x**3','2*x**3*y'],  str(((x**2)*((x+y)**2-x**2-y**2)).expand())
    assert str((((x+y)**2-x**2-y**2)*(x**2)).expand()) in ['2*y*x**3','2*x**3*y'], str((((x+y)**2-x**2-y**2)*(x**2)).expand())
    two = ((x+y)**2-x**2-y**2)/x/y
    assert str((two).expand())=='2', str(two.expand())
    assert str((two*x).expand())=='2*x'
    assert str((two*(2*x)).expand())=='4*x'
    assert str((two*(x+y)).expand()) in ['2*x + 2*y', '2*y + 2*x'], str((two*(x+y)).expand())
    assert str((x*two).expand())=='2*x'
    assert str(((x-y)*two).expand()) in ['2*x - 2*y', '-2*y + 2*x'],str(((x-y)*two).expand())

    two_y = ((x+y)**2-x**2-y**2)/x
    assert str((two_y).expand())=='2*y'
    assert str((two_y*x).expand()) in ['2*x*y','2*y*x'],str((two_y*x).expand())
    assert str((two_y*(2*x)).expand()) in ['4*x*y','4*y*x'], str((two_y*(2*x)).expand())
    assert str((two_y*(x+y)).expand()) in ['2*y**2 + 2*x*y','2*x*y + 2*y**2',
                                           '2*y**2 + 2*y*x','2*y*x + 2*y**2'],str((two_y*(x+y)).expand())
    assert str((x*two_y).expand()) in ['2*x*y','2*y*x'], str((x*two_y).expand())
    assert str(((x-y)*two_y).expand()) in ['2*x*y - 2*y**2','-2*y**2 + 2*y*x',
                                           '2*y*x - 2*y**2'],str(((x-y)*two_y).expand())

    x2 = ((x+y)**2-x**2-y**2)/y*x/2
    assert str((x2).expand())=='x**2'
    assert str((2*x2).expand())=='2*x**2'
    assert str((x*x2).expand())=='x**3'
    assert str((x**2*x2).expand())=='x**4'
    assert str((y*x2).expand()) in ['y*x**2','x**2*y'],str((y*x2).expand())

    assert str(((x+y)*(x+y+z)).expand()) in ['x*z + y*z + 2*x*y + x**2 + y**2',
                                             'x**2 + y**2 + 2*x*y + x*z + y*z',
                                             'x**2 + y**2 + 2*y*x + x*z + y*z',
                                             'x**2 + y**2 + y*z + 2*y*x + x*z',
                                             'x**2 + x*z + 2*x*y + z*y + y**2',
                                             'y*z + 2*y*x + x**2 + y**2 + x*z'],str(((x+y)*(x+y+z)).expand())
    assert str(((x+y+z)*(x+y)).expand()) in ['x*z + y*z + 2*x*y + x**2 + y**2',
                                             'x**2 + y**2 + 2*x*y + x*z + y*z',
                                             'x**2 + y**2 + 2*y*x + x*z + y*z',
                                             'x**2 + y**2 + y*z + 2*y*x + x*z',
                                             'x**2 + x*z + 2*x*y + z*y + y**2',
                                             'y*z + 2*y*x + x**2 + y**2 + x*z'], str(((x+y+z)*(x+y)).expand())
    assert str(((1/x+x)*x).expand()) in ['1 + x**2', 'x**2 + 1'], str(((1/x+x)*x).expand())
    assert str((x**2*(1/x+x)**2).expand()) in ['1 + 2*x**2 + x**4',
                                               '1 + x**4 + 2*x**2',
                                               'x**4 + 1 + 2*x**2',
                                               'x**4 + 2*x**2 + 1'],str((x**2*(1/x+x)**2).expand())

    r = ' + '.join(sorted(str(((1+x+y)**3).expand()).split(' + ')))

    assert r in ['1 + x**3 + y**3 + 3*x + 3*x**2 + 3*x*y**2 + 3*y + 3*y**2 + 3*y*x**2 + 6*x*y',
                 '1 + 3*x**2 + x**3 + 3*y**2*x + 3*y + 3*x + 3*y**2 + 3*y*x**2 + 6*y*x + y**3',
                 '1 + 3*x + 3*x**2 + 3*y + 3*y**2 + 3*y**2*x + 3*y*x**2 + 6*y*x + x**3 + y**3',
                 '1 + 3*x + 3*x**2 + 3*x**2*y + 3*x*y**2 + 3*y + 3*y**2 + 6*x*y + x**3 + y**3',
                 #'1 + 2*y*x + 3*x + 3*x**2 + 3*y + 3*y**2 + 3*y**2*x + 3*y*x**2 + 4*y*x + x**3 + y**3'
                 ], r

def _test_has_symbol():
    x, y = map(Symbol, 'xy')
    assert (3+x+2**y)._get_symbols_data() == set(['x', 'y'])
    assert (3+x+2**y).symbols == set([x, y])
    assert x.has_symbol(x)
    assert y.has_symbol(y)
    assert not x.has_symbol(y)
    assert (1+x).has_symbol(x)
    assert not (1+x).has_symbol(y)
    assert (3*x).has_symbol(x)
    assert not (3*x).has_symbol(y)
    assert (3*x*y).has_symbol(x)
    assert (3*x*y).has_symbol(y)
    assert (x**2).has_symbol(x)
    assert (2**x).has_symbol(x)
    assert (2**(1+3*x)).has_symbol(x)
    assert (y + 2**(1+3*x)).has_symbol(x)

def _test_to_polynomial_data():
    x, y = map(Symbol,'xy')
    assert x.to_polynomial_data([x])[0] == {1:1}
    assert x.to_polynomial_data([x], True)[0] == {1:1}
    assert x.to_polynomial_data([x,y])[0] == {(1,0):1}
    assert x.to_polynomial_data([x,y], True)[0] == {(1,0):1}
    assert x.to_polynomial_data([y])[0] == {(0,1):1}
    assert x.to_polynomial_data([y], True)[0] == {0:x}

    assert (2*x).to_polynomial_data([x])[0] == {1:2}
    assert (2*x).to_polynomial_data([x,y])[0] == {(1,0):2}
    assert (2*x).to_polynomial_data([x,y],True)[0] == {(1,0):2}
    assert (2*x).to_polynomial_data([y])[0] == {(0,1):2}
    assert (2*x).to_polynomial_data([y], True)[0] == {0:2*x}

    assert (x*y**2).to_polynomial_data([x])[0] == {(1,2):1}
    assert (x*y**2).to_polynomial_data([x], True)[0] == {1:y**2}
    assert (x*y**2).to_polynomial_data([x,y])[0] == {(1,2):1}
    assert (x*y**2).to_polynomial_data([x,y], True)[0] == {(1,2):1}
    assert (x*y**2).to_polynomial_data([y])[0] == {(2,1):1}
    assert (x*y**2).to_polynomial_data([y], True)[0] == {2:x}

    assert (3*x*y**2).to_polynomial_data([x])[0] == {(1,2):3}
    assert (3*x*y**2).to_polynomial_data([x],True)[0] == {1:3*y**2}
    assert (3*x*y**2).to_polynomial_data([x,y])[0] == {(1,2):3}
    assert (3*x*y**2).to_polynomial_data([x,y],True)[0] == {(1,2):3}
    assert (3*x*y**2).to_polynomial_data([y])[0] == {(2,1):3}
    assert (3*x*y**2).to_polynomial_data([y],True)[0] == {2:3*x}

    assert (x+2*y).to_polynomial_data([x])[0] == {(1,0):1,(0,1):2}
    assert (x+2*y).to_polynomial_data([x], True)[0] == {1:1,0:2*y}
    assert (x+2*y).to_polynomial_data([x,y])[0] == {(1,0):1,(0,1):2}
    assert (x+2*y).to_polynomial_data([x,y], True)[0] == {(1,0):1,(0,1):2}
    assert (x+2*y).to_polynomial_data([y])[0] == {(1,0):2,(0,1):1}
    assert (x+2*y).to_polynomial_data([y], True)[0] == {1:2,0:x}

    assert (x**(Number(3)/2)).to_polynomial_data()[0] == {3:1}

def test_bug_issue61():
    a = Symbol('a')
    assert (a+1-1)==a
