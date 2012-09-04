""" Provides CommutativeRing class.
"""

__docformat__ = "restructuredtext"
__all__ = ['CommutativeRing']

from .algebra import Algebra

class CommutativeRing(Algebra):
    """ Base class to commutative rings.

    Derived classes may redefine the following methods::
    
      Symbol(cls, obj)
      Number(cls, obj)
      Add(cls, *seq)
      Mul(cls, *seq),
      Pow(cls, base, exponent)
      Terms(cls, *seq)
      Factors(cls, *seq)
      as_Add_args(self)
      as_Mul_args(self)
      as_Pow_args(self),
      as_Terms_args(self)
      as_Factors_args(self)
    """
    __slots__ = ['_symbols']
    _symbols = None


    @classmethod
    def Pos(cls, arg): return +arg

    @classmethod
    def Neg(cls, arg): return -arg

    @classmethod
    def Add(cls, *seq):
        """ Compute sum over seq containing algebra elements.
        """
        raise NotImplementedError('%s must define classmethod Add'    #pragma NO COVER
                                  % (cls.__name__))                   #pragma NO COVER

    @classmethod
    def Sub(cls, *seq):
        """ Compute ``seq[0] - Add(*seq[0])``.
        """
        if seq:
            return seq[0] - cls.Add(*seq[1:])
        return cls.zero
    
    @classmethod
    def Mul(cls, *seq):
        """ Compute product over seq containing algebra elements.
        """
        raise NotImplementedError('%s must define classmethod Mul'    #pragma NO COVER
                                  % (cls.__name__))                   #pragma NO COVER

    @classmethod
    def Div(cls, *seq):
        """ Compute ``seq[0] / Mul(*seq[0])``.
        """
        if seq:
            return seq[0] / cls.Mul(*seq[1:])
        return cls.one

    @classmethod
    def Pow(cls, base, exponent):
        """ Compute power from base and exponent.
        
        Argument base must be an algebra element and exponent must be
        an element of exponent algebra.
        """
        raise NotImplementedError('%s must define classmethod Pow'     #pragma NO COVER
                                  % (cls.__name__))                    #pragma NO COVER

    @classmethod
    def Log(cls, arg, base=None):
        """ Compute logarithm of arg in base.

        Argument arg must be an element of exponent algebra and base
        is an element of an algebra.
        """
        raise NotImplementedError('%s must define classmethod Pow'     #pragma NO COVER
                                  % (cls.__name__))                    #pragma NO COVER

    @classmethod
    def Terms(cls, *seq):
        """ Compute sum over seq containing pairs ``(element, coefficient)``.

        ``element``-s must belong to algebra, ``coefficient``-s must
        belong to the coefficient algebra.
        """
        raise NotImplementedError('%s must define classmethod Terms'   #pragma NO COVER
                                  % (cls.__name__))                    #pragma NO COVER

    @classmethod
    def Factors(cls, *seq):
        """ Compute product over seq containing pairs ``(element, exponent)``.

        ``element``-s must belong to algebra, ``exponent``-s must belong
        to the exponent algebra.
        """
        raise NotImplementedError('%s must define classmethod Factors' #pragma NO COVER
                                  % (self.__class__.__name__))         #pragma NO COVER

    def as_Add_args(self):
        """ Return a sequence such that ``Add(*self.as_Add_args()) == self``
        """
        raise NotImplementedError('%s must define method as_Add_args'  #pragma NO COVER
                                  % (self.__class__.__name__))         #pragma NO COVER
    
    def as_Mul_args(self):
        """ Return a sequence such that ``Mul(*self.as_Mul_args()) == self``
        """
        raise NotImplementedError('%s must define method as_Mul_args'  #pragma NO COVER
                                  % (self.__class__.__name__))         #pragma NO COVER
    
    def as_Pow_args(self):
        """ Return a 2-tuple such that ``Pow(*self.as_Pow_args()) == self``
        """
        raise NotImplementedError('%s must define method as_Pow_args'  #pragma NO COVER
                                  % (self.__class__.__name__))         #pragma NO COVER

    def as_Log_args(self):
        """ Return a 2-tuple such that ``Log(*self.as_Log_args()) == self``
        """
        raise NotImplementedError('%s must define method as_Log_args'  #pragma NO COVER
                                  % (self.__class__.__name__))         #pragma NO COVER

    def as_Terms_args(self):
        """ Return a sequence such that ``Terms(*self.as_Terms_args()) == self``
        """
        raise NotImplementedError('%s must define method as_Terms_args' #pragma NO COVER
                                  % (self.__class__.__name__))          #pragma NO COVER

    def as_Factors_args(self):
        """ Return a sequence such that ``Factors(*self.as_Factors_args()) == self``
        """
        raise NotImplementedError('%s must define method as_Factors_args' #pragma NO COVER
                                  % (self.__class__.__name__))            #pragma NO COVER

    def __pos__(self):
        return self

    def __neg__(self):
        return self.Mul(self, self.convert(-1))

    def __add__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self.Add(self, other)

    def __radd__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self.Add(other, self)

    def __sub__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self + (-other)

    def __rsub__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return other + (-self)

    def __mul__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self.Mul(self, other)

    def __rmul__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self.Mul(other, self)

    def __div__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        if isinstance(other, (int, long)):
            other = self.Number(other)
        return self * other ** (-1)

    def __rdiv__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return other * self ** (-1)

    def __truediv__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self * other ** (-1)

    def __rtruediv__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return other * self ** (-1)

    def __pow__(self, other):
        other = self.convert_exponent(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self.Pow(self, other)

    def __rpow__(self, other):
        other = self.convert(other, False)
        if other is NotImplemented:
            return NotImplemented
        return self.Pow(other,  self.convert_exponent(self))

    def _matches(pattern, expr, repl_dict, wild_info):
        pfunc = pattern.func
        efunc = expr.func
        Add = pattern.Add
        Mul = pattern.Mul
        if pfunc==Add:
            wild_part = []
            exact_part = []
            for a in pattern.args:
                if a.symbols.intersection(wild_info[0]):
                    if a.head is SYMBOL:
                        wild_part.append(a)
                    else:
                        wild_part.insert(0, a)
                else:
                    exact_part.append(a)
            if exact_part:
                expr = expr - Add(*exact_part)
            if len(wild_part)==1:
                return wild_part[0].matches(expr, repl_dict, wild_info)
            expr_args = expr.as_Add_args() + [pattern.zero]
            for i in range(len(wild_part)):
                w = wild_part[i]
                w_rest = Add(*(wild_part[:i]+wild_part[i+1:]))
                for e in expr_args:
                    r = w.matches(e, repl_dict, wild_info)
                    if r is not None:
                        r1 = w_rest.subs(r.items()).matches(expr-e, r, wild_info)
                        if r1 is not None:
                            return r1
        elif pfunc==Mul:
            wild_part = []
            exact_part = []
            for a in pattern.args:
                if a.symbols.intersection(wild_info[0]):
                    if a.head is SYMBOL:
                        wild_part.append(a)
                    else:
                        wild_part.insert(0, a)
                else:
                    exact_part.append(a)
            if exact_part:
                expr = expr / Mul(*exact_part)
            if len(wild_part)==1:
                return wild_part[0].matches(expr, repl_dict, wild_info)
            expr_args = expr.as_Mul_args() + [pattern.one]
            for i in range(len(wild_part)):
                w = wild_part[i]
                w_rest = Mul(*(wild_part[:i]+wild_part[i+1:]))
                for e in expr_args:
                    r = w.matches(e, repl_dict, wild_info)
                    if r is not None:
                        r1 = w_rest.subs(r.items()).matches(expr/e, r, wild_info)
                        if r1 is not None:
                            return r1
        elif pfunc == efunc:
            return pattern._matches_seq(pattern.args, expr.args, repl_dict, wild_info)
        return

from ..utils import NUMBER, SYMBOL
