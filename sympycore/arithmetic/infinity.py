"""Provides base class for extended numbers: Infinity.
"""
__docformat__ = "restructuredtext"
__all__ = ['Infinity']

from ..core import init_module, Expr
init_module.import_heads()

from .numbers import realtypes, complextypes, mpqc, numbertypes, div

def is_undefined(obj):
    if isinstance(obj, Infinity):
        return not obj.data
    return False

def is_zoo(obj):
    if isinstance(obj, Infinity):
        return is_undefined(obj.data)
    return False

class Infinity(object):
    """ Base class for directional infinities.

    The formal definition of directional infinity is

    ::

      Infinity(direction) == oo * direction

    where ``oo`` represents a symbol of a formal limiting process
    ``limit(r, r->oo)``.

    Direction can be

      * 0                - defines ``undefined``
      * complex number   - defines direction in complex plain
      * undefined        - defines infinity with undefined direction (``zoo``)

    Derived classes may redefine classmethods::

      IsUnbounded(x)        - returns 1 if x is unbounded expression
      IsZero(x)             - returns 1 if x is zero
      EqualArg(x,y)         - returns 1 if Arg(x)==Arg(y).
      
    to add algebra support (including symbolic functions) to direction
    expression.

    The following notation is used::

      +oo = Infinity(1)
      -oo = Infinity(-1)
      undefined = Infinity(0)
      zoo = Infinity(undefined)

    """

    one = 1
    zero = 0

    _add_ops = None

    @classmethod
    def EqualArg(cls, x, y):
        if isinstance(x, numbertypes) and isinstance(y, numbertypes):
            if x and y:
                if isinstance(x, complextypes):
                    if isinstance(y, complextypes):
                        return x.real * y.imag == x.imag * y.real
                    else:
                        return False
                if isinstance(y, complextypes):
                    return False
                return cmp(x,0)==cmp(y,0)
            return False
        if is_zoo(x) or is_zoo(y) or is_undefined(x) or is_undefined(y):
            return False
        if x==y:
            return True
        if not x or not y:
            return False
        if isinstance(x, Expr): algebra = type(x)
        elif isinstance(y, Expr): algebra = type(y)
        else: raise NotImplementedError(`cls, x, y`)
        return algebra.Apply(cls.EqualArg, x, y)

    @classmethod
    def IsUnbounded(cls, x):
        if isinstance(x, Expr):
            head, data = x.pair
            if head is NUMBER:
                return cls.IsUnbounded(data)
            is_bounded = getattr(x, 'is_bounded', None)
            if is_bounded is not None:
                return not is_bounded
            y = x.evalf(2)
            if y.head is NUMBER:
                return cls.IsUnbounded(y.data)
            return APPLY.new(type(x), (cls.IsUnbounded, (x,)))
        if isinstance(x, Infinity):
            return True
        if isinstance(x, numbertypes):
            # XXX: handle mpf, nan, inf
            return False
        if isinstance(x, Expr):
            algebra = type(x)
        else:
            return NotImplemented
        return algebra.Apply(cls.IsUnbounded, x)
    
    @classmethod
    def IsPositive(cls, x):
        if isinstance(x, Expr):
            head, data = x.pair
            if head is NUMBER:
                return cls.IsPositive(data)
            y = x.evalf(2)
            if y.head is NUMBER:
                return cls.IsPositive(y.data)
            return APPLY.new(type(x), (cls.IsPositive, (x,)))
        if isinstance(x, numbertypes):
            return cmp(x, 0)==1
        if isinstance(x, cls):
            return cls.IsPositive(x.data)
        if isinstance(x, Expr): algebra = type(x)
        else: return NotImplemented
        return algebra.Apply(cls.IsPositive, x)

    def __new__(cls, direction):
        obj = object.__new__(cls)
        if isinstance(direction, Expr):
            head, data = direction.pair
            if head is NUMBER:
                return cls(data)
            elif head is TERM_COEFF:
                t, c = data
                if c>0:
                    return cls(t)
            r = direction.get_direction()
            if r is not NotImplemented:
                return cls(r)
        if is_undefined(direction):
            obj.data = direction
        elif isinstance(direction, realtypes):
            obj.data = cmp(direction, 0)
        elif isinstance(direction, complextypes):
            re, im = direction.real, direction.imag
            if re:
                if abs(im) > abs(re):
                    obj.data = mpqc(re, im)/im
                else:
                    obj.data = mpqc(re, im)/re
            elif im:
                obj.data = mpqc(re, im)/im
            else:
                obj.data = 0
        else:
            obj.data = direction
        return obj

    def __repr__(self):
        return '%s(%r)' % (type(self).__name__, self.data)

    def __str__(self):
        d = self.data
        if is_zoo(self):
            return 'zoo'
        if is_undefined(self):
            return 'undefined'
        if d==1:
            return 'oo'
        if d==-1:
            return '-oo'
        a = self._add_ops
        if a is not None:
            lhs, rhs = a
            return '%s + (%s)' % (lhs, rhs)
        return 'oo*(%s)' % (d,)

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.data == other.data
        return False

    def __abs__(self):
        """ ``abs(oo*x) -> oo*abs(x)``
        """
        if is_undefined(self):
            return self
        return type(self)(1) # +oo

    def __neg__(self):
        """ negation of extended number
        
        ::

          -(oo*U1(x)) -> oo*(-U1(x))
          -zoo        -> zoo
        """
        if is_zoo(self) or is_undefined(self):
            return self
        return type(self)(-self.data)

    def __add__(self, other):
        """ addition of extended numbers

        ::

          oo*x + oo*y     -> oo*(EqualArg(x, y) * x)
          oo*x + y        -> oo*((1+IsUnbounded(y)*(EqualArg(x,y)-1))*x)
          y + oo*x        -> oo*x + y
        """
        cls = type(self)
        if is_undefined(self) or is_undefined(other) or is_zoo(other):
            return cls(0)
        if is_zoo(self) and isinstance(other, cls):
            return cls(0)
        x = self.data
        if isinstance(other, cls):
            y = other.data
            r = cls(cls.EqualArg(x, y) * x)
            r._add_ops = (self, other)
        else:
            y = other
            y_unbounded = cls.IsUnbounded(y)
            if not y_unbounded:
                return self
            z = (1+y_unbounded*(cls.EqualArg(x, y)-1))
            r = cls(z * x)
            r._add_ops = (self, other)
        return r

    __radd__ = __add__

    def __sub__(self, other):
        """ substraction of extended numbers

        ::

          oo*x - y -> oo*x + (-y)
        """
        return self + (-other)

    def __rsub__(self, other):
        """ substraction of extended numbers

        ::

          y - oo*x -> oo*(-x) + y
        """
        return (-self) + other

    def __mul__(self, other):
        """ multiplication of extended numbers

        ::
        
          (oo*x) * (oo*y) -> oo*(x*y)
          (oo*x) * y          -> oo*(x*y)
        """
        cls = type(self)
        if is_undefined(self) or is_undefined(other) or not other:
            return cls(0)
        if isinstance(other, Infinity):
            return cls(self.data * other.data)
        return cls(self.data * other)

    __rmul__ = __mul__

    def __div__(self, other):
        """ division of extended numbers

        ::

          (oo*x) / (oo*y)     -> undefined
          (oo*x) / 0          -> zoo
          (oo*x) / y          -> oo*(x/y)
        """
        cls = type(self)
        if is_undefined(self) or is_undefined(other) or isinstance(other, cls):
            return cls(0)
        if not other:
            return cls(cls(0))
        if isinstance(other, Expr):
            head, data = other.pair
            if head is NUMBER:
                other = data
            else:
                head, data = other.evalf(2).pair
                if head is NUMBER:
                    other = data
        if isinstance(other, realtypes):
            if other < 0:
                return cls(-self.data)
            return self
        return cls(self.data / other)

    def __rdiv__(self, other):
        """ division of extended numbers

        ::

          y / oo*x -> 0*(1/x)
        """
        if is_undefined(self):
            return self
        return type(other)(0)

    def __pow__(self, other):
        """ exponentiation of extended numbers

        ::

          (oo*x)**0  -> 1
          (oo*x)**(oo*y) -> 0 if y<0; oo*(IsPositive(x)) if y > 0
          (oo*x)**y      -> 0 if y<0; oo*(x**y) if y > 0
        """
        if not other:
            return self.one
        cls = type(self)
        if is_undefined(self) or is_undefined(other) or is_zoo(other):
            return cls(0)
        x = self.data
        zero = self.zero
        algebra = None
        if isinstance(other, Expr):
            algebra = type(other)
            head, data = other.pair
            if head is NUMBER:
                other = data
        if isinstance(other, realtypes):
            if other < 0:
                return zero
            if isinstance(x, realtypes):
                if x > 0:
                    return cls(1)
                return cls((-1)**other)
            return cls(x ** other)
        if isinstance(other, cls):
            y = other.data
            if isinstance(y, realtypes):
                if y < 0:
                    return zero
                if isinstance(x, realtypes):
                    if x > 0:
                        return cls(1)
                    if y==1:
                        return cls(cls(0))
                if is_zoo(self):
                    return self
                return cls(cls.IsPositive(x))
        if algebra is not None:
            exp = algebra(other)
            return algebra(POW, (algebra(SPECIAL, self), exp))
        return NotImplemented

    def __rpow__(self, other):
        """ exponentiation of extended numbers

        ::

          1 ** (oo*x) -> 1
          y ** (oo*x) -> (y**x)**(oo)
        """
        cls = type(self)
        x = self.data
        zero = self.zero
        one = self.one
        algebra = None
        if isinstance(other, Expr):
            algebra = type(other)
            head, data = other.pair
            if head is NUMBER:
                other = data
        if isinstance(other, realtypes):
            if other==1:
                return one
            elif other>1:
                if x==1:
                    return self
                if x==-1:
                    return zero
            elif other>-1:
                if x==1:
                    return zero
                if x==-1:
                    return cls(cls(0))
            elif other==-1:
                return cls(0)
            else:
                if x==1:
                    return cls(cls(0))
                if x==-1:
                    return zero
            if is_zoo(self):
                return cls(0)
            if is_undefined(self):
                return self
        if algebra is not None:
            base = algebra(other)
            return algebra(POW, (base, algebra(SPECIAL, self)))
        return NotImplemented

    def __lt__(self, other):
        if is_undefined(self) or is_undefined(other) or is_zoo(self) or is_zoo(other):
            return False
        x = self.data
        if isinstance(other, type(self)):
            y = other.data
            if x in [-1,1] and y in [-1,1]:
                return x < y
        if isinstance(other, realtypes):
            if x in [-1, 1]:
                return 2*x < cmp(other, 0)
        return NotImplemented

    __le__ = __lt__

    def __gt__(self, other):
        if is_undefined(self) or is_undefined(other) or is_zoo(self) or is_zoo(other):
            return False
        x = self.data
        if isinstance(other, type(self)):
            y = other.data
            if x in [-1,1] and y in [-1,1]:
                return x > y
        if isinstance(other, realtypes):
            if x in [-1, 1]:
                return 2*x > cmp(other, 0)
        return NotImplemented

    __ge__ = __gt__

    @classmethod
    def get_oo(cls):
        return cls(1)

    @classmethod
    def get_moo(cls):
        return cls(-1)

    @classmethod
    def get_zoo(cls):
        return cls(cls(0))

    @classmethod
    def get_undefined(cls):
        return cls(0)

    def is_undefined(self):
        return self.data==0
    
