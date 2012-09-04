"""
This module implements efficient "low-level" number types for purely
numerical tasks.

We use the following types:

int/long -- (Python built-ins) exact integers
mpf -- arbitrary-precision floats
mpc -- arbitrary-precision complex floats
mpq -- nonintegral fractions
mpqc -- nonreal complex rationals

In the interest of speed, there are some quirks:

  * mpq derives from tuple. mpq.__init__ does *not* validate its input.
    mpqs *should only* be created from fully reduced (p, q) pairs,
    and p/q should *not* be integer-valued.

  * To safely create a fraction from two integers, use the function
    normalized_fraction()

  * Arithmetic operations on mpqs automatically normalize back to Python
    ints when the results are integer-valued. Likewise, mpqc with zero
    real part normalize to their real parts.

  * Note that ``mpq(2)/mpq(3)`` does *not* work as expected; it does the
    same thing as ``2/3`` in Python. To perform safe division, use the
    provided ``div`` function.

  * mpqs can be compared with Python ints, but cannot be (correctly)
    compared to Python floats. (Use the mpf class instead.)

Powers, except when the exponent is a positive integer, should be
computed with the ``try_power()`` function which detects when the
result is not exact.
"""
#
# Author: Fredrik Johansson
# Created: January 2008

import math

from . import mpmath

__docformat__ = "restructuredtext"

from ..core import init_module

init_module.import_heads()

@init_module
def _init(module):
    from ..basealgebra.verbatim import Verbatim
    module.Verbatim = Verbatim

from ..utils import str_SUM, str_PRODUCT, str_POWER, str_APPLY, str_SYMBOL, str_NUMBER, NUMBER, SYMBOL

inttypes = (int, long)

from .mpmath import mpf, mpc, mp
from .mpmath.libmp import from_rational, round_nearest



def mpf_to_str_data(self, sort=True):
    if self < 0:
        return str_SUM, str(self)
    return str_NUMBER, str(self)

def mpc_to_str_data(self, sort=True):
    return str_NUMBER, str(self)

mpf.to_str_data = mpf_to_str_data
mpc.to_str_data = mpc_to_str_data

rounding = round_nearest

def getdps():
    return mp.dps

def setdps(n):
    p = mp.dps
    mp.dps = int(n)
    return p

def getprec():
    return mp.prec

#----------------------------------------------------------------------------
# Fractions
#

def normalized_fraction(p, q=1):
    """ Return a normalized fraction.
    """
    x, y = p, q
    while y:
        x, y = y, x % y
    if x != 1:
        p //= x
        q //= x
    if q == 1:
        return p
    return mpq((p, q))

class mpq(tuple):
    """Represents a fraction."""

    # These methods are inherited directly from tuple for speed. This works
    # as long as all mpqs are normalized:
    # __new__/__init__
    # __nonzero__
    # __eq__
    # __hash__

    # These methods are generated and defined in methods.py:
    # __add__, __sub__, __rsub__, __mul__, __div__, __rdiv__, __pow__
    # __lt__, __le__, __gt__, __ge__

    __slots__ = []

    @property
    def _mpf_(self):
        p, q = self
        return from_rational(p, q, getprec(), rounding)

    def as_verbatim(self):
        p, q = self
        if p<0:
            return -(Verbatim(NUMBER, -p) / Verbatim(NUMBER, q))
        return Verbatim(NUMBER, p) / Verbatim(NUMBER, q)

    def to_str_data(self,sort=True):
        if self[0]<0:
            return str_SUM, str(self)
        return str_PRODUCT, str(self)

    def __str__(self):
        return "%i/%i" % self

    def __repr__(self):
        p, q = self
        return "%s((%s, %s))" % (type(self).__name__, p, q)

    def __float__(self):
        p, q = self
        return float(p) / q

    def __int__(self):
        p, q = self
        return p // q

    def __neg__(self):
        p, q = self
        return mpq((-p, q))

    def __pos__(self):
        return self

    def __abs__(self):
        p, q = self
        if p < 0:
            return mpq((-p, q))
        return self

    def __floordiv__(a, b):
        return int(a / b)

    def __rfloordiv__(a, b):
        return int(b / a)

    def __mod__(a, b):
        return a - (a//b)*b

    def __rmod__(a, b):
        return b - (b//a)*a

    def __divmod__(a, b):
        return (a-a%b)/b, a%b

    def __rdivmod__(a, b):
        return (b-b%a)/a, b%a

    def __rpow__(a, b):
        z, sym = try_power(b, a)
        if not sym:
            return z
        return NotImplemented


#----------------------------------------------------------------------------
# Complex numbers
#

def innerstr(x):
    if isinstance(x, mpq):
        return "%s/%s" % x
    return str(x)


class mpqc(object):
    """Represents an exact complex number.

    The integer power of a complex number is computed using
    `exponentiation by squaring`__ method.

    __ http://en.wikipedia.org/wiki/Exponentiation_by_squaring
    """

    __slots__ = ['real', 'imag']

    def __new__(cls, real, imag=0):
        if not imag:
            return real
        if isinstance(real, (float, mpf)) or isinstance(imag, (float, mpf)):
            return mpc(real, imag)
        self = object.__new__(mpqc)
        self.real = real
        self.imag = imag
        return self

    def __repr__(self):
        return "%s(%r, %r)" % (self.__class__.__name__, self.real, self.imag)

    def __hash__(self):
        return hash((self.real, self.imag))

    def __mpcval__(self):
        return mpc(self.real, self.imag)

    def __getstate__(self):
        return (self.real, self.imag)

    def __setstate__(self, (real, imag)):
        self.real = real
        self.imag = imag

    def to_str_data(self,sort=True):
        re, im = self.real, self.imag
        if re==0:
            if im == 1: return str_SYMBOL,"I"
            if im == -1: return str_SUM, "-I"
            return str_PRODUCT, str(self.imag) + "*I"
        restr = innerstr(self.real)
        if im == 1: return str_SUM, "%s + I" % restr
        if im == -1: return str_SUM, "%s - I" % restr
        if im > 0: return str_SUM, "%s + %s*I" % (restr, innerstr(self.imag))
        if im < 0: return str_SUM, "%s - %s*I" % (restr, innerstr(-self.imag))
        raise NotImplementedError(`self`)

    def __str__(self):
        return self.to_str_data()[1]

    def as_verbatim(self):
        re, im = self.real, self.imag
        if re < 0: 
            r = -Verbatim(NUMBER, -self.real)
        else:
            r = Verbatim(NUMBER, self.real)
        if im < 0:
            i = -Verbatim(NUMBER, -self.imag)
            ni = Verbatim(NUMBER, -self.imag)
        else:
            i = Verbatim(NUMBER, self.imag)
        I = Verbatim(NUMBER, mpqc(0,1))
        if re==0:
            if im == 1: return I
            if im == -1: return -I
            return i * I
        if im == 1: return r + I
        if im == -1: return r - I
        if im < 0:
            return r - ni * I
        else:
            return r + i * I

    def as_numer_denom(self):
        re, im = self.real, self.imag
        if type(re) is mpq:
            r1, r2 = re
        else:
            r1, r2 = re, 1
        if type(im) is mpq:
            i1, i2 = im
        else:
            i1, i2 = im, 1
        r = r1*i2
        i = i1*r2
        d = r2*i2
        c = gcd(r,i,d)
        return mpqc(r//c,i//c), d//c
            
    def __eq__(self, other):
        if hasattr(other, "imag"):
            return self.real == other.real and self.imag == other.imag
        return False

    def __pos__(self): return self
    def __neg__(self): return mpqc(-self.real, -self.imag)

    def __abs__(self):
        re, im = self.real, self.imag
        if not re:
            return abs(im)
        m2 = re**2 + im**2
        if isinstance(m2, inttypes):
            m, flag = int_root(m2,2)
            if flag:
                return m
        if isinstance(m2, mpq):
            p,fp = int_root(m2[0],2)
            q,fq = int_root(m2[1],2)
            if fp and fq:
                return mpq((p,q))
        raise NotImplementedError('abs(%r)' % (self))

#----------------------------------------------------------------------------
# Interface functions
#

inttypes = (int, long)
rationaltypes = (mpq,)
numbertypes = (int, long, float, complex, mpq, mpqc, mpf, mpc)
realtypes = (int, long, float, mpq, mpf)
complextypes = (complex, mpqc, mpc)

inttypes_set = frozenset(inttypes)
rationaltypes_set = frozenset(rationaltypes)
realtypes_set = frozenset(realtypes)
complextypes_set = frozenset(complextypes)
numbertypes_set = frozenset(numbertypes)

def number_div(Algebra, a, b):
    if type(b) in inttypes_set:
        if not b:
            if a:
                return Infinity(Infinity(0))
            return Infinity(0)
        if b == 1:
            return a
        if type(a) in inttypes_set:
            return normalized_fraction(a, b)
    return a / b

def div(a, b):
    """Safely compute a/b.

    If a or b is an integer, this function makes sure to convert it to
    a rational.
    """
    if type(b) in inttypes_set:
        if not b:
            return Infinity(a)
            raise ZeroDivisionError('%r / %r' % (a, b))
        if b == 1:
            return a
        if type(a) in inttypes_set:
            return normalized_fraction(a, b)
    return a / b

def int_root(y, n):
    """ Return a pair ``(floor(y**(1/n)), x**n == y)``.

    Given integers y and n, return a tuple containing ``x =
    floor(y**(1/n))`` and a boolean indicating whether ``x**n == y``
    exactly.
    """
    if y < 0: raise ValueError, "y must not be negative"
    if n < 1: raise ValueError, "n must be positive"
    if y in (0, 1): return y, True
    if n == 1: return y, True
    if n > y: return 1, False
    # Get initial estimate for Newton's method. Care must be taken to
    # avoid overflow
    try:
        guess = int(y ** (1./n)+0.5)
    except OverflowError:
        try:
            guess = int(math.exp(math.log(y)/n)+0.5)
        except OverflowError:
            guess = 1 << int(math.log(y, 2)/n)
    # Newton iteration
    xprev, x = -1, guess
    while 1:
        t = x**(n-1)
        xprev, x = x, x - (t*x-y)//(n*t)
        if abs(x - xprev) < 1:
            break
    # Compensate
    t = x**n
    while t > y:
        x -= 1
        t = x**n
    return x, t == y

def try_power(x, y):
    """\
    Attempt to compute ``x**y`` where ``x`` and ``y`` must be of the
    types int, long, mpq, mpf, or mpqc. The function
    returns::

      z, symbolic

    where ``z`` is a number (i.e. a complex rational) and ``symbolic``
    is a list of ``(b, e)`` pairs representing symbolic factors
    ``b**e``.

    Examples::
    
      try_power(3, 2) --> (9, [])
      try_power(2, 1/2) --> (1, [(2, 1/2)])
      try_power(45, 1/2) --> (3, [(5, 1/2)])
    """
    if not y or x == 1:
        # (anything)**0 -> 1
        # 1**(anything) -> 1
        return 1, []
    if isinstance(x, Infinity) or isinstance(y, Infinity):
        return x**y, []
    if isinstance(y, inttypes):
        if y >= 0:
            return x**y, []
        elif not x:
            return Infinity.get_zoo(), []
        elif isinstance(x, inttypes):
            return mpq((1, x**(-y))), []
        return x**y, []
    elif isinstance(x, inttypes) and isinstance(y, mpq):
        if x < 0:
            if x==-1:
                p, q = y
                if q==2:
                    return mpqc(0, 1)**p, []
                return 1, [(x,y)]
            else:
                z, sym = try_power(-x, y)
                z1, sym1 = try_power(-1, y)
                return z * z1, sym+sym1
        else:
            p, q = y
            r, exact = int_root(x, q)
            if exact:
                if r==1 or not p:
                    return 1, []
                g = r**p if p>0 else mpq((1, r**(-p)))
                return g, []
    elif isinstance(x, mpq) and isinstance(y, mpq):
        a, b = x
        r, rsym = try_power(a, y)
        s, ssym = try_power(b, y)
        ssym = [(b, -e) for b, e in ssym]
        return (div(r,s), rsym + ssym)
    elif isinstance(x, mpf) or isinstance(y, mpf):
        return x ** y, []
    return 1, [(x, y)]

from .number_theory import gcd
from .evalf import evalf
from .infinity import Infinity

from .methods import (\
    fraction_add, fraction_sub, fraction_rsub, fraction_mul,
    fraction_div, fraction_rdiv, fraction_pow,
    fraction_lt, fraction_le, fraction_gt, fraction_ge,
    complex_add, complex_sub, complex_rsub, complex_mul,
    complex_div, complex_rdiv, complex_pow)

mpq.__add__ = mpq.__radd__ = fraction_add
mpq.__sub__ = fraction_sub
mpq.__rsub__ = fraction_rsub
mpq.__mul__ = mpq.__rmul__ = fraction_mul
mpq.__div__ = fraction_div
mpq.__rdiv__ = fraction_rdiv
mpq.__pow__ = fraction_pow
mpq.__lt__ = fraction_lt
mpq.__le__ = fraction_le
mpq.__gt__ = fraction_gt
mpq.__ge__ = fraction_ge

mpqc.__add__ = mpqc.__radd__ = complex_add
mpqc.__sub__ = complex_sub
mpqc.__rsub__ = complex_rsub
mpqc.__mul__ = mpqc.__rmul__ = complex_mul
mpqc.__div__ = complex_div
mpqc.__rdiv__ = complex_rdiv
mpqc.__pow__ = complex_pow
