"""Provides algorithms from number theory.
"""

from ..core import init_module

init_module.import_lowlevel_operations()

from .numbers import mpq, div
from . import mpmath

__all__ = ['gcd', 'lcm', 'factorial',
           'integer_digits', 'real_digits',
           'multinomial_coefficients', 'f2q']

__docformat__ = "restructuredtext en"

def factorial(n, memo=[1, 1]):
    """Return n factorial (for integers n >= 0 only)."""
    if n < 0:
        raise ValueError("expected non-negative integer, got %r" % (n)) #pragma NO COVER
    k = len(memo)
    if n < k:
        return memo[n]
    p = memo[-1]
    while k <= n:
        p *= k
        k += 1
        if k < 100:
            memo.append(p)
    return p

def gcd(*args):
    """Calculate the greatest common divisor (GCD) of the arguments."""
    L = len(args)
    if L == 0: return 0
    if L == 1: return args[0]
    if L == 2:
        a, b = args
        while b:
            a, b = b, a % b
        return a
    return gcd(gcd(args[0], args[1]), *args[2:])

def lcm(*args):
    """Calculate the least common multiple (LCM) of the arguments."""
    L = len(args)
    if L == 0: return 0
    if L == 1: return args[0]
    if L == 2: return div(args[0], gcd(*args))*args[1]
    return lcm(lcm(args[0], args[1]), *args[2:])

# TODO: this could use the faster implementation in mpmath
def integer_digits(n, base=10):
    """Return a list of the digits of abs(n) in the given base."""
    assert base > 1
    assert isinstance(n, (int, long))
    n = abs(n)
    if not n:
        return [0]
    L = []
    while n:
        n, digit = divmod(n, base)
        L.append(int(digit))
    return L[::-1]

# TODO: this could (also?) be implemented as an endless generator
def real_digits(x, base=10, truncation=10):
    """Return ``(L, d)`` where L is a list of digits of ``abs(x)`` in
    the given base and ``d`` is the (signed) distance from the leading
    digit to the radix point.

    For example, 1234.56 becomes ``([1, 2, 3, 4, 5, 6], 4)`` and 0.001
    becomes ``([1], -2)``. If, during the generation of fractional
    digits, the length reaches `truncation` digits, the iteration is
    stopped."""
    assert base > 1
    assert isinstance(x, (int, long, mpq))
    if x == 0:
        return ([0], 1)
    x = abs(x)
    exponent = 0
    while x < 1:
        x *= base
        exponent -= 1
    integer, fraction = divmod(x, 1)
    L = integer_digits(integer, base)
    exponent += len(L)
    if fraction:
        p, q = fraction
        for i in xrange(truncation - len(L)):
            p = (p % q) * base
            if not p:
                break
            L.append(int(p//q))
    return L, exponent

def binomial_coefficients(n):
    """Return a dictionary containing pairs {(k1,k2) : C_kn} where
    C_kn are binomial coefficients and n=k1+k2.

    INPUT:
        n -- an integer

    OUTPUT:
        dict
 
    EXAMPLES:
        >>> sorted(binomial_coefficients(3).items())
        [((0, 3), 1), ((1, 2), 3), ((2, 1), 3), ((3, 0), 1)]
 
    Notice the coefficients above are the same as below:
        x**3 + 3*x**2*y + 3*x*y**2 + y^3

    """
    d = {(0, n):1, (n, 0):1}
    a = 1
    for k in xrange(1, n//2+1):
        a = (a * (n-k+1))//k
        d[k, n-k] = d[n-k, k] = a
    return d

def binomial_coefficients_list(n):
    """Return a list of binomial coefficients as rows of the Pascal's
    triangle.
    """
    d = [1] * (n+1)
    a = 1
    for k in xrange(1, n//2+1):
        a = (a * (n-k+1))//k
        d[k] = d[n-k] = a
    return d

def multinomial_coefficients(m, n, _tuple=tuple, _zip=zip):
    """Return a dictionary containing pairs ``{(k1,k2,..,km) : C_kn}``
    where ``C_kn`` are multinomial coefficients such that
    ``n=k1+k2+..+km``.

    INPUT:
        m -- integer
        n -- integer
        _tuple, _zip -- hacks for speed; don't set these as a user.

    OUTPUT:
        dict

    EXAMPLES:
        >>> sorted(multinomial_coefficients(2,5).items())
        [((0, 5), 1), ((1, 4), 5), ((2, 3), 10), ((3, 2), 10), ((4, 1), 5), ((5, 0), 1)]

    Notice that these are the coefficients of ``(x+y)**5``:
         x**5 + 5*x**4*y + 10*x**3*y**2 + 10*x**2*y**3 + 5*x*y**4 + y**5
 
        >>> sorted(multinomial_coefficients(3,2).items())
        [((0, 0, 2), 1), ((0, 1, 1), 2), ((0, 2, 0), 1), ((1, 0, 1), 2), ((1, 1, 0), 2), ((2, 0, 0), 1)]


    ALGORITHM:
    The algorithm we implement for computing the multinomial coefficients
    is based on the following result:
    
       Consider a polynomial and its ``n``-th exponent::
       
         P(x) = sum_{i=0}^m p_i x^k
         P(x)^n = sum_{k=0}^{m n} a(n,k) x^k

       We compute the coefficients ``a(n,k)`` using the
       J.C.P. Miller Pure Recurrence [see D.E.Knuth, Seminumerical
       Algorithms, The art of Computer Programming v.2, Addison
       Wesley, Reading, 1981;]::
       
         a(n,k) = 1/(k p_0) sum_{i=1}^m p_i ((n+1)i-k) a(n,k-i),

       where ``a(n,0) = p_0^n``.
    """
    if m==0:
        return {}
    if m==2:
        return binomial_coefficients(n)
    symbols = [(0,)*i + (1,) + (0,)*(m-i-1) for i in range(m)]
    s0 = symbols[0]
    p0 = [_tuple([aa-bb for aa,bb in _zip(s,s0)]) for s in symbols]
    r = {_tuple([aa*n for aa in s0]):1}
    r_get = r.get
    r_update = r.update
    l = [0] * (n*(m-1)+1)
    l[0] = r.items()
    for k in xrange(1, n*(m-1)+1):
        d = {}
        d_get = d.get
        for i in xrange(1, min(m,k+1)):
            nn = (n+1)*i-k
            if not nn:
                continue
            t = p0[i]
            for t2, c2 in l[k-i]:
                tt = _tuple([aa+bb for aa,bb in _zip(t2,t)])
                cc = nn * c2
                dict_add_item(None, d, tt, cc)
        r1 = [(t, c//k) for (t, c) in d.iteritems()]
        l[k] = r1
        r_update(r1)
    return r

def reldiff(x,y):
    """abs(x-y)/((abs(x)+abs(y))/2)"""
    return abs(x-y)/((abs(x)+abs(y))/2)

def f2q(f, minerr=None):
    '''
    Find (almost) the best rational representation of a float that has
    finite precision.

    The algorithm is based on Stern-Brocot representation of
    irrational/rational numbers and its relation to continuants [1].

    References:
    [1] "Concrete Mathematics" by Graham, Knuth, and Patashnik. Addison-Wesley. 1989
    '''
    if mpmath.libmp.BACKEND=='gmpy':
        tmpz = type(mpmath.libmp.MPZ_ZERO)
        if minerr is None:
            r = mpmath.libmp.gmpy.f2q(f)
        else:
            r = mpmath.libmp.gmpy.f2q(f, minerr)
        if isinstance (r, tmpz):
            return int(r)
        n,d = int(r.numer ()), int(r.denom ())
        if d==1:
            return n
        return mpq((n,d))
    mpf = mpmath.mpf
    f = mpf(f)
    if f<0: return -f2q(-f,minerr)
    al = f                          # Initialize Eq.(6.142) [1]
    if not minerr:
        minerr = 1/mpf(2**al.context.prec) # min number for the given precision
    elif minerr<0:
        minerr = 1/mpf(2**-minerr)
    a = mpmath.floor(al)   # a = floor(al) # Initialize Eq.(6.142) [1]
    r1=[0,0,1]                           # r1[1:] is 1st row of Eq.(6.138) [1] 
    r2=[0,1,a]                           # r2[1:] is 2nd row of Eq.(6.138) [1]
    err = reldiff(f,a)
    while err>minerr:
        al = 1/(al-a)                    # Eq.(6.142) [1]
        a = mpmath.floor(al) # a = floor(al) # Eq.(6.142) [1]
        r1[0] = r1[1]                    # roll
        r1[1] = r1[2]                    # roll
        r1[2] = r1[1]*a+r1[0]            # use Eq.(6.127) [1]
        r2[0] = r2[1]                    # roll
        r2[1] = r2[2]                    # roll
        r2[2] = r2[1]*a+r2[0]            # use Eq.(6.127) [1]
        newerr = reldiff(f,r2[2]/r1[2])
        if err<=newerr:
            n, d = (int(r2[1]),int(r1[1]))
            if d==1: return n
            return mpq((n, d))
        err = newerr
    n, d = (int(r2[2]),int(r1[2]))
    if d==1: return n
    return mpq((n, d))


#from .combinatorics import multinomial_coefficients
