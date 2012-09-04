#
# Created January 2008 by Pearu Peterson
#
""" Provides Calculus class.
"""
__docformat__ = "restructuredtext"
__all__ = ['Calculus', 'I']

from ..core import classes, defined_functions
from ..utils import TERMS, str_PRODUCT, SYMBOL, NUMBER
from ..heads import APPLY, CALLABLE, TERM_COEFF, TERM_COEFF_DICT
from ..heads import BASE_EXP_DICT as FACTORS
from ..basealgebra import Algebra, Verbatim
from ..ring import CommutativeRing

from ..arithmetic.numbers import normalized_fraction, mpq, mpf, mpc, mpqc, try_power

from ..arithmetic import mpmath, setdps
from ..arithmetic.evalf import evalf

algebra_numbers = (int, long, mpq, mpqc, mpf, mpc)
convertible_numbers = algebra_numbers + (float, complex)

float_one = mpf('1.0')

class Calculus(CommutativeRing):
    """ Represents an element of a symbolic algebra.

    The set of a symbolic algebra is a set of expressions.
    """

    _hash = None

    coefftypes = algebra_numbers
    coefftypes_set = frozenset(coefftypes)
    exptypes = algebra_numbers
    exptypes_set = frozenset(exptypes)

    @classmethod
    def get_function_algebra(cls):
        return classes.CalculusFunctionRing

    def __abs__(self):
        head, data = self.pair
        if head is NUMBER:
            return type(self)(NUMBER, abs(data))
        raise NotImplementedError('abs(%r)' % (self))

    def __call__(self, *args, **kws):
        head, data = self.pair
        if head is CALLABLE:
            return data(*args, **kws)
        raise TypeError(`self` + ' is not callable')

    def __lt__(self, other):
        if type(other) is type(self):
            return self.data < other.data
        return self.data < other
    def __le__(self, other):
        if type(other) is type(self):
            return self.data <= other.data
        return self.data <= other
    def __gt__(self, other):
        if type(other) is type(self):
            return self.data > other.data
        return self.data > other
    def __ge__(self, other):
        if type(other) is type(self):
            return self.data >= other.data
        return self.data >= other
    def __ne__(self, other):
        return not (self == other)
    
    def _sympy_(self):
        import sympy as m
        head, data = self.pair
        if head is SYMBOL:
            return m.Symbol(data)
        if head is NUMBER:
            return m.Number(data)
        if head is TERMS:
            return m.Add(*self.as_Add_args())
        if head is FACTORS:
            return m.Mul(*self.as_Mul_args())
        raise NotImplemented(`self`)

    def as_algebra(self, cls, typeerror=True):
        """ Convert algebra to another algebra.
        """
        if cls is Verbatim:
            return self.as_verbatim()
        if cls is classes.Unit:
            return cls(NUMBER, self)
        if issubclass(cls, PolynomialRing):
            return self.as_polynom(cls)
        return self.as_verbatim().as_algebra(cls)

    @classmethod
    def get_predefined_symbols(cls, name):
        if name=='I':
            return I
        if name=='D':
            return classes.DFactory()
        return getattr(defined_functions, name, None)

    @classmethod
    def get_defined_function(cls, name):
        func = getattr(defined_functions, name, None)
        if func is None:
            func = getattr(defined_functions, name.title(), None)
        return func
    
    @classmethod
    def convert_coefficient(cls, obj, typeerror=True):
        """ Convert obj to coefficient algebra.
        """
        if isinstance(obj, float):
            return mpf(obj)
        if isinstance(obj, complex):
            return mpc(obj.real, obj.imag)
        if isinstance(obj, algebra_numbers):
            return obj
        if isinstance(obj, cls):
            head, data = obj.pair
            if head is NUMBER:
                return data
        if typeerror:
            raise TypeError('%s.convert_coefficient: failed to convert %s instance'\
                            ' to coefficient algebra, expected int|long object'\
                            % (cls.__name__, obj.__class__.__name__))
        else:
            return NotImplemented

    @classmethod
    def convert_exponent(cls, obj, typeerror=True):
        """ Convert obj to exponent algebra.
        """
        if isinstance(obj, cls):
            return obj
        if isinstance(obj, algebra_numbers):
            return obj
        if isinstance(obj, float):
            return mpf(obj)
        if isinstance(obj, complex):
            return mpc(obj.real, obj.imag)

        # parse algebra expression from string:
        if isinstance(obj, (str, unicode, Verbatim)):
            return Verbatim.convert(obj).as_algebra(cls)

        # convert from another algebra:
        if isinstance(obj, Algebra):
            return obj.as_algebra(cls)

        if typeerror:
            raise TypeError('%s.convert_exponent: failed to convert %s instance'\
                            ' to exponent algebra, expected int|long object'\
                            % (cls.__name__, obj.__class__.__name__))
        else:
            return NotImplemented
    
    @classmethod
    def Number(cls, num, denom=None):
        if denom is None:
            return cls(NUMBER, num)
        return cls(NUMBER, normalized_fraction(num, denom))

    @classmethod
    def Log(cls, arg, base=None):
        log = defined_functions.Log
        if base is None:
            return log(arg)
        return log(arg)/log(base)

    @classmethod
    def Exp(cls, arg):
        return defined_functions.Exp(arg)

    @classmethod
    def Mod(cls, x, y):
        return defined_functions.Mod(x, y)

    def evalf(self, n=None):
        if n is not None:
            setdps(n)
        head, data = self.pair
        def evalf(cls, head, data, target):
            if head is NUMBER:
                data = data * float_one
            if head is SYMBOL:
                if hasattr(data, 'evalf'):
                    return cls(data.evalf(n))
            if head is APPLY:
                func, args = data
                name = str(func).lower()
                if hasattr(mpmath, name):
                    f = getattr(mpmath, name)
                    l = []
                    for a in args:
                        h, d = a.pair
                        if h is NUMBER:
                            l.append(d)
                        else:
                            l = None
                            break
                    if l is not None:
                        return cls(f(*l))
            if head is TERM_COEFF:
                term, coeff = data
                coeff = coeff * float_one
                return cls(head, (term, coeff))
            if head is TERM_COEFF_DICT:
                new_data = {}
                for term, coeff in data.iteritems():
                    new_data[term] = coeff * float_one
                return cls(head, new_data)
            return target
        return head.walk(evalf, type(self), data, self)
        
        if n:
            setdps(n)
        head, data = self.pair
        if head is NUMBER:
            return self.Number(data * float_one)
        if head is SYMBOL:
            r = getattr(data, 'evalf', lambda p: NotImplemented)(n)
            if r is not NotImplemented:
                return self.Number(r)
            return self
        if head is APPLY:
            func = data[0]
            h, func_data = func.pair
            assert h is CALLABLE, `self, func` # todo: support symbolic functions
            args = data[1]
            assert len (args)==1,`self, args` # todo: support multivariate functions
            v = args[0].evalf(n)
            h, d = v.pair
            if h is NUMBER:
                return self.Number(getattr(mpmath, func_data.__name__.lower())(d))
            else:
                return type(self)(APPLY, (func, (v,)))
        convert = self.convert
        return self.func(*[convert(a).evalf(n) for a in self.args])

    def get_direction(self):
        head, data = self.pair
        if head is NUMBER:
            if isinstance(data, (int, long)):
                return data
            return getattr(data, 'get_direction', lambda : NotImplemented)()
        if head is TERMS:
            if len(data)==1:
                t, c = data.items()[0]
                r = t.get_direction()
                if r is not NotImplemented:
                    return r * c
        if head is FACTORS:
            direction = 1
            cls = type(self)
            for t,c in self.data.iteritems():
                d = t.get_direction()
                if d is NotImplemented:
                    return d
                if not isinstance(c, (int, long)):
                    return NotImplemented
                d = self.Pow(cls.convert(d), c).get_direction()
                if d is NotImplemented:
                    return d
                direction *= d
            return direction
        return getattr(data, 'get_direction', lambda : NotImplemented)()

    @property
    def is_bounded(self):
        head, data = self.pair
        if head is NUMBER:
            if isinstance(data, (int, long)):
                return True
            return getattr(data, 'is_bounded', None)
        if head is SYMBOL:
            return getattr(data, 'is_bounded', None)
        if head is TERMS:
            for t, c in data.iteritems():
                b = t.is_bounded
                if not b:
                    return b
                if isinstance(c, (int, long)):
                    continue
                b = getattr(c, 'is_bounded', None)
                if not b:
                    return b
            return True
        return

    def as_polynom(self, ring_cls=None):
        """ Convert expression to an element of polynomial ring.
        
        If the polynomial ring is not given then it will be created.

        For example,

          >>> x,y,z = map(Symbol,'xyz')
          >>> (x+2*y+3*z).as_polynom() + y*z
          PolynomialRing[(x, y, z), Calculus]('2*x + y*z + y + 3*z')
          >>> P = PolynomialRing[x,y]
          >>> (x+2*y+3*z).as_polynom(P) + y*z
          PolynomialRing[(x, y), Calculus]('x + (2 + z)*y + 3*z')

        """
        if ring_cls is None:
            data, variables = self.to_polynomial_data()
            return PolynomialRing[tuple(variables)].convert(data)
        else:
            data, variables = self.to_polynomial_data(ring_cls.variables, True)
            return ring_cls.convert(data)

    def __divmod__(self, other):
        if isinstance(other, Calculus):
            lhs = self.as_polynom()
            rhs = other.as_polynom(type(lhs))
            return divmod(lhs, rhs)
        return NotImplemented

    def __float__(self):
        head, data = self.pair
        if head is NUMBER:
            return float(data)
        head, data = self.evalf().pair
        if head is NUMBER:
            return float(data)
        raise TypeError('Cannot convert %r to float' % (self))

    def __or__(self, other):
        if type(other) is tuple:
            return self.subs(*other)
        return self.subs(other)

classes.Calculus = Calculus

# as a result, x<y will return Logic(LT, (x, y)):
Calculus.enable_symbolic_comparison('inequality')

one = Calculus.Number(1)
zero = Calculus.Number(0)
Calculus.one = one
Calculus.zero = zero

I = Calculus.Number(mpqc(0,1))

from ..polynomials.algebra import PolynomialRing, AdditiveTuple
from .infinity import CalculusInfinity
