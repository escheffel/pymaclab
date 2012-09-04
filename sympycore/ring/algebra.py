
__all__ = ['Ring', 'CommutativeRing']

from ..basealgebra import Algebra
from .interface import RingInterface

from ..core import init_module, classes
init_module.import_heads()
init_module.import_numbers()

@init_module
def _init(m):
    from ..arithmetic import mpq
    Ring.coefftypes = (int, long, mpq)

class Ring(Algebra, RingInterface):
    """
    Ring represents algebraic ring (R, +, *) where (R, +) is abelian
    group, (R,*) is monoid, with distributivity.
    """

    @classmethod
    def get_function_algebra(cls):
        return classes.FunctionRing

    def __str__(self):
        h, d = self.pair
        return h.data_to_str_and_precedence(type(self), d)[0]

    def __pos__(self):
        return self

    def __neg__(self):
        return self.head.neg(type(self), self)

    def __add__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.add_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.add(cls, self, other)

    __radd__ = __add__

    def __iadd__(self, other):
        cls = type(self)
        if type(other) is not cls:
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.inplace_add(cls, self, other)

    def __sub__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.sub_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.sub(cls, self, other)

    def __rsub__(self, other):
        return other + (-self)

    def __isub__(self, other):
        cls = type(self)
        if type(other) is not cls:
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.inplace_add(cls, self, -other)

    def __mul__(self, other):
        cls = type(self)
        tother = type(other)
        if cls is not tother:
            if tother in numbertypes_set:
                return self.head.non_commutative_mul_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.non_commutative_mul(cls, self, other)

    def __rmul__(self, other):
        cls = type(self)
        tother = type(other)
        if cls is not tother:
            if tother in numbertypes_set:
                return self.head.non_commutative_rmul_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return other.head.non_commutative_mul(cls, other, self)

    def __pow__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.pow_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.pow(cls, self, other)

    def __rpow__(self, other):
        cls = type(self)
        tother = type(other)
        if cls is not tother:
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return other.head.pow(cls, other, self)

    def __div__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.non_commutative_div_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self * other**-1

    def __rdiv__(self, other):
        cls = type(self)
        tother = type(other)
        if cls is not tother:
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return other * self**-1

    __truediv__ = __div__

    def expand(self):
        return self.head.expand(type(self), self)

    def evalf(self, n=None):
        return self.head.evalf(type(self), self, n)

class CommutativeRing(Ring):

    def __mul__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.commutative_mul_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented:
                return NotImplemented
        return self.head.commutative_mul(cls, self, other)

    __rmul__ = __mul__

    def __imul__(self, other):
        cls = type(self)
        if type(other) is not cls:
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented:
                return NotImplemented
        return self.head.inplace_commutative_mul(cls, self, other)

    def __div__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.commutative_div_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented:
                return NotImplemented
        return self.head.commutative_div(cls, self, other)

    def __rdiv__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.commutative_rdiv_number(cls, self, other)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented:
                return NotImplemented
        return other * self**-1

    def to(self, target, *args):
        """ Convert expression to target representation.

        The following targets are recognized:

          EXP_COEFF_DICT - convert expression to exponents-coefficient
            representation, additional arguments are variables. When
            no arguments are specified, variables will be all symbols
            and non-power expressions used in the expression.

          TERM_COEFF_DICT - convert expression to term-coefficient
            representation. Note that the returned result may have
            actual head NUMBER, SYMBOL, TERM_COEFF, POW, or
            BASE_EXP_DICT instead of TERM_COEFF_DICT.
        """
        head, data = self.pair
        if target is head:
            return self
        if target is EXP_COEFF_DICT:
            return head.to_EXP_COEFF_DICT(type(self), data, self, args or None)
        if target is TERM_COEFF_DICT:
            return head.to_TERM_COEFF_DICT(type(self), data, self)
        raise NotImplementedError('%s.to(target=%r)' % (type(self), target))

    def diff(self, symbol, order=1):
        if order==0:
            return self
        cls = type(self)
        if type(symbol) is cls:
            assert symbol.head is SYMBOL,`symbol.pair`
            symbol = symbol.data
        elif isinstance(symbol, str):
            pass
        else:
            raise TypeError('diff(symbol, order) first argument must be str or %s instance but got %s instance' % (cls.__name__, type(symbol).__name__))
        try:
            cache = {}
            result = self.head.diff(cls, self.data, self, symbol, order, cache=cache)
        finally:
            cache.clear()
        return result

    def integrate(self, x):
        cls = type(self)

        t = type(x)
        if t is tuple:
            x, a, b = x
            t = type(x)
        else:
            a, b = None, None

        if t is cls:
            assert x.head is SYMBOL,`x.pair`
            x = x.data
        elif t is str:
            pass
        else:
            raise TypeError('integrate(x,..), x must be str or %s instance but got %s instance' % (cls.__name__, type(symbol).__name__))
    
        if a is None:
            return self.head.integrate_indefinite(cls, self.data, self, x)
        if type(a) is not cls:
            a = cls(a)
        if type(b) is not cls:
            b = cls(b)
        return self.head.integrate_definite(cls, self.data, self, x, a, b)

classes.Ring = Ring
classes.CommutativeRing = CommutativeRing
