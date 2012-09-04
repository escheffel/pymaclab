""" Provides CalculusInfinity class.
"""
__docformat__ = "restructuredtext"
__all__ = ['CalculusInfinity', 'oo', 'moo', 'undefined', 'zoo']

from ..utils import str_PRODUCT, NUMBER, TERMS, FACTORS, SYMBOL
from ..heads import SPECIAL
from ..arithmetic import Infinity
from .algebra import Calculus

CalculusInfinity = Infinity

class _CalculusInfinity(Infinity):

    one = Calculus.one
    zero = Calculus.zero

    def __new__(cls, direction):
        if isinstance(direction, Calculus):
            head, data = direction.head, direction.data
            if head is NUMBER:
                return cls(data)
            elif head is TERMS:
                if len(data)==1:
                    t, c = data.items()[0]
                    if c > 0:
                        return cls(t)
            r = direction.get_direction()
            if r is not NotImplemented:
                return cls(r)
        return Infinity.__new__(cls, direction)

    def to_str_data(self, sort=True):
        return str_PRODUCT, str(self)

    def _get_symbols_data(self):
        return set()

    def subs(self, *args):
        d = self.data
        if isinstance(d, Calculus):
            return type(self)(d.subs(*args))
        return self

    @classmethod
    def IsUnbounded(cls, x):
        if isinstance(x, Calculus):
            head = x.head
            if head is NUMBER:
                return cls.IsUnbounded(x.data)
            f = x.is_bounded
            if f or f is not None:
                return Calculus.convert(not f)
        if isinstance(x, cls):
            return Calculus.one
        r = Infinity.IsUnbounded(x)
        if r is not NotImplemented:
            return Calculus.convert(r)
        x = Calculus.convert(x)
        return Calculus.Apply(cls.IsUnbounded, x)

    @classmethod
    def EqualArg(cls, x, y):
        if isinstance(x, Calculus):
            d = x.get_direction()
            if d is not NotImplemented:
                return cls.EqualArg(d, y)
        if isinstance(y, Calculus):
            d = y.get_direction()
            if d is not NotImplemented:
                return cls.EqualArg(x, d)
            if isinstance(x, Calculus):
                xy = x / y
                if xy.head is NUMBER:
                    return Calculus.convert(xy.data > 0)

        r = Infinity.EqualArg(x, y)
        if r is not NotImplemented:
            return Calculus.convert(r)
        x = Calculus.convert(x)
        y = Calculus.convert(y)
        return Calculus.Apply(cls.EqualArg, x,y)

    @classmethod
    def IsPositive(cls, x):
        if isinstance(x, Calculus):
            if x.head is NUMBER:
                return cls.IsPositive(x.data)
        r = Infinity.IsPositive(x)
        if r is not NotImplemented:
            return Calculus.convert(r)
        x = Calculus.convert(x)
        return Calculus.apply(cls.IsPositive, x)

    def __pow__(self, other):
        if isinstance(other, Calculus):
            head, data = other.head, other.data
            if head is NUMBER:
                other = data
        r = Infinity.__pow__(self, other)
        if r is not NotImplemented:
            if isinstance(r, type(self)):
                return r
            return Calculus.convert(r)
        x = Calculus.convert(other)
        return Calculus(FACTORS, {Calculus(SPECIAL, self): x})

    def __rpow__(self, other):
        if isinstance(other, Calculus):
            if other.head is SYMBOL:
                other = other.evalf(2)
            head, data = other.head, other.data
            if head is NUMBER:
                other = data
        r = Infinity.__rpow__(self, other)
        if r is not NotImplemented:
            if isinstance(r, type(self)):
                return r
            return Calculus.convert(r)
        x = Calculus.convert(other)
        return Calculus(FACTORS, {x: self})

oo = CalculusInfinity(1)
moo = CalculusInfinity(-1)
undefined = CalculusInfinity(0)
zoo = CalculusInfinity(undefined)

Calculus.oo = oo
Calculus.zoo = zoo
Calculus.undefined = undefined
Calculus.Infinity = CalculusInfinity
