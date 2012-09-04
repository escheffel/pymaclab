#
# Author: Fredrik Johansson
# Created: January 2008
""" Provides UnivariatePolynomial class.
"""
__docformat__ = "restructuredtext"
__all__ = ['UnivariatePolynomial', 'poly']

from ..utils import ADD, POW, MUL, NUMBER, SYMBOL, DENSE_POLY
from ..basealgebra.verbatim import Verbatim
from ..basealgebra import Algebra
from ..arithmetic.numbers import div

def poly(coefs=[], symbol='x', coerce_input=True):
    coefs = coefs or [0]
    if coerce_input:
        coefs = map(UnivariatePolynomial.convert_coef, coefs)
    return UnivariatePolynomial(DENSE_POLY, (symbol, normalized(coefs)))

def normalized(coefs):
    while coefs and not coefs[-1]:
        coefs = coefs[:-1]
    return tuple(coefs)

class UnivariatePolynomial(Algebra):
    """ Represents (dense) univariate polynomial.
    """

    @classmethod
    def convert(cls, data, typeerror=True):
        if isinstance(data, (list, tuple)):
            coefs = data or [0]
            return cls(DENSE_POLY, ('x', normalized(coefs)))
        return super(Algebra, cls).convert(data, typeerror=True)

    def __nonzero__(self):
        return not not self.data[1]

    @property
    def degree(self):
        return len(self.data[1])-1

    @classmethod
    def convert_coef(cls, x):
        # permit anything by default
        return x

    def as_verbatim(self):
        if self.degree == -1:
            return Verbatim(NUMBER, 0)
        t = []
        head, (symbol, data) = self.pair
        x = Verbatim(SYMBOL, symbol)
        for exp, coef in enumerate(data):
            if coef == 0:
                continue
            scoef = Verbatim(NUMBER, coef)
            if exp == 0:
                monomial = scoef
            else:
                monomial = x
                if exp != 1: monomial = monomial ** Verbatim(NUMBER, exp)
                if coef != 1: monomial = scoef * monomial
            t.append(monomial)
        if t:
            return Verbatim(ADD, tuple(t))
        return Verbatim(NUMBER, 0)

    @property
    def tree(self):
        return self.as_verbatim().tree

    def __call__(self, x):
        head, (symbol, data) = self.pair
        cls = type(self)
        if isinstance(x, UnivariatePolynomial):
            p = cls(x.head, (symbol,()))
            xp = cls(x.head, (symbol,(1,)))
            for i, c in enumerate(data):
                p += c * xp
                xp *= x
            return p
        else:
            y = cls.convert_coef(0)
            for c in reversed(data):
                y = y*x+c
            return y

    def __neg__(self):
        head, (symbol, data) = self.pair
        return type(self)(head, (symbol, tuple([-c for c in data])))

    def __add__(self, other):
        head, (symbol, data) = self.pair
        cls = type(self)
        if not isinstance(other, UnivariatePolynomial):
            other = cls(head, (symbol,(other,)))
        ohead, (osymbol, odata) = other.pair
        if symbol != osymbol:
            return NotImplemented
        coefs = [a+b for (a, b) in zip(data, odata)]
        if self.degree >= other.degree:
            coefs += data[other.degree+1:]
        else:
            coefs += odata[self.degree+1:]
        return cls(head, (symbol, normalized(coefs)))

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        head, (symbol, data) = self.pair
        cls = type(self)
        if not isinstance(other, UnivariatePolynomial):
            try:
                other = self.convert_coef(other)
                coefs = [other*c for c in data]
                return cls(head, (symbol, normalized(coefs)))
            except:
                pass
            return NotImplemented
        ohead, (osymbol, odata) = other.pair
        if osymbol != osymbol:
            return NotImplemented
        p = [self.convert_coef(0)] * (len(data)+len(odata))
        for i, c in enumerate(data):
            for j, d in enumerate(odata):
                p[i+j] += c*d
        return cls(head, (symbol, normalized(p)))

    __rmul__ = __mul__

    def __pow__(self, n):
        # TODO: Miller's algorithm
        assert isinstance(n, int) and n >= 0
        head, (symbol, data) = self.pair
        if n == 0: return type(self)(head, (symbol, (1,)))
        if n == 1: return self
        if n & 1 == 0: return (self*self)**(n//2)
        if n & 1 == 1: return self*((self*self)**(n//2))
        return NotImplemented

    def __divmod__(self, other):
        cls = type(self)
        head, (symbol, data) = self.pair
        if isinstance(other, UnivariatePolynomial):
            ohead, (osymbol, odata) = other.pair
            if symbol != osymbol:
                return NotImplemented
            if other.degree < 0:
                raise ZeroDivisionError, "polynomial division"
            n = self.degree
            nv = other.degree
            u = data
            v = odata
            r = list(u)
            q = [0] * (len(r)+1)
            for k in range(n-nv, -1, -1):
                q[k] = div(r[nv+k], v[nv])
                for j in range(nv+k-1, k-1, -1):
                    r[j] -= q[k]*v[j-k]
            for j in range(nv, n+1, 1):
                r[j] = 0
            return cls(head, (symbol, normalized(q))), cls(head, (symbol, normalized(r)))
        else:
            if other == 0:
                raise ZeroDivisionError, "polynomial division"
            rec = div(1, other)
            q = cls(head, (symbol, tuple([c*rec for c in data])))
            r = cls(head, (symbol, ()))
            return q, r

    def __div__(self, other):
        return divmod(self, other)[0]

    __truediv__ = __div__

    def __mod__(self, other):
        return divmod(self, other)[1]

    def diff(self):
        head, (symbol, data) = self.pair
        return type(self)(head, (symbol, tuple([c*(k+1) for k, c in enumerate(data[1:])])))

