
from ..core import classes
from ..utils import EQ, NE, LT, LE, GT, GE, SYMBOL, AND, NOT, OR, NUMBER, IN, NOTIN
from ..basealgebra import Algebra, Verbatim

head_mth_map = {
    EQ: Algebra.__eq__,
    NE: Algebra.__ne__,
    LT: Algebra.__lt__,
    LE: Algebra.__le__,
    GT: Algebra.__gt__,
    GE: Algebra.__ge__,
    }

class Logic(Algebra):
    """ Represents n-ary predicate expressions.

    An expression ``predicate(*args)`` is hold in a pair::

      (<predicate constant or function>, args)

    Examples of predicates::

      Logic(NUMBER, True)  - 0-ary predicate, represents truth value True
      Logic(NUMBER, False) - 0-ary predicate, represents truth value False
      Logic(SYMBOL, 'x')   - 0-ary predicate, represents a boolean symbol x
      Logic(NOT, x)        - 1-ary predicate, represents ``not x`` boolean expression
      Logic(AND, frozenset([x,y,z])) - represents ``x and y and z`` boolean expression
      Logic(OR, frozenset([x,y])) - represents boolean expression ``x or y``
      Logic(LT, (x, y))    - represents relation ``x < y``
      Logic(LE, (x, y))    - represents relation ``x <= y``
      Logic(GT, (x, y))    - represents relation ``x > y``
      Logic(GT, (x, y))    - represents relation ``x >= y``
      Logic(EQ, (x, y))    - represents relation ``x == y``
      Logic(NE, (x, y))    - represents relation ``x != y``
      Logic(IN, (x, y))    - represents relation ``x in y``
      Logic(NOTIN, (x, y)) - represents relation ``x not in y``

    """

    @classmethod
    def get_operand_algebra(cls, head, index=0):
        """ Return algebra class for index-th operand of operation with head.
        """
        if head in [AND, OR, NOT]:
            return cls
        if head in [LT, LE, GT, GE, EQ, NE]:
            return classes.Calculus
        if head in [IN, NOTIN]:
            if index==1:
                return classes.Set
            return classes.Verbatim
        raise NotImplementedError('Algebra %s does not support operation %s' % (cls.__name__, head))

    def __nonzero__(self, head_mth_get=head_mth_map.get):
        """ Return lexicographic truth value for relational predicates
        and True for non-numeric predicates.
        """
        head, data = self.pair
        return head.nonzero(type(self), data)
        mth = head_mth_get(head)
        if mth is None:
            if head is NUMBER:
                return data
            return True
        return mth(*data)

    def as_verbatim(self):
        return Verbatim(self.head, self.data)

    @classmethod
    def get_predefined_symbols(cls, name):
        if name=='True': return cls.true
        if name=='False': return cls.false
        return

    @classmethod
    def convert_number(cls, obj, typeerror=True):
        if type(obj) is bool:
            return obj
        return cls.handle_convert_failure('number', obj, typeerror, 'expected bool')

    def convert_operand(self, obj, typeerror=True):
        head = self.head
        if head in [LT, GT, LE, GE, EQ, NE]:
            return classes.Calculus.convert(obj, typeerror)
        if isinstance(obj, bool):
            return type(self)(NUMBER, obj)
        return self.convert(obj, typeerror)

    is_Lt = property(lambda self: self.head is LT)
    is_Le = property(lambda self: self.head is LE)
    is_Gt = property(lambda self: self.head is GT)
    is_Ge = property(lambda self: self.head is GE)
    is_Eq = property(lambda self: self.head is EQ)
    is_Ne = property(lambda self: self.head is NE)
    is_And = property(lambda self: self.head is AND)
    is_Or = property(lambda self: self.head is OR)
    is_Not = property(lambda self: self.head is NOT)

    def _subs(self, subexpr, newexpr):
        head, data = self.pair
        cls = type(self)

        if type(subexpr) is not cls:
            if head is SYMBOL or head is NUMBER:
                if subexpr==data:
                    return self.convert_operand(newexpr)
                return self
        else:
            h, d = subexpr.pair
            if h is head and d==data:
                return self.convert_operand(newexpr)

        if head is SYMBOL or head is NUMBER:            
            return self
        if head is AND:
            return cls.And(*[a._subs(subexpr, newexpr) for a in data])        
        if head is OR:
            return cls.Or(*[a._subs(subexpr, newexpr) for a in data])
        if head is NOT:
            return cls.Not(data._subs(subexpr, newexpr))
        if head is LT:
            return cls.Lt(*[a._subs(subexpr, newexpr) for a in data])
        if head is GT:
            return cls.Gt(*[a._subs(subexpr, newexpr) for a in data])
        if head is LE:
            return cls.Le(*[a._subs(subexpr, newexpr) for a in data])
        if head is GE:
            return cls.Ge(*[a._subs(subexpr, newexpr) for a in data])
        if head is EQ:
            return cls.Eq(*[a._subs(subexpr, newexpr) for a in data])
        if head is NE:
            return cls.Ne(*[a._subs(subexpr, newexpr) for a in data])
        return self

    @classmethod
    def Lt(cls, *seq):
        if len(seq)==2:
            lhs, rhs = seq
            if lhs.head is NUMBER and rhs.head is NUMBER:
                return cls.true if lhs.data < rhs.data else cls.false
            if not lhs < rhs:
                return cls(GT, (rhs, lhs))
        return cls(LT, seq)

    @classmethod
    def Le(cls, *seq):
        if len(seq)==2:
            lhs, rhs = seq
            if lhs.head is NUMBER and rhs.head is NUMBER:
                return cls.true if lhs.data <= rhs.data else cls.false
            if not lhs < rhs:
                return cls(GE, (rhs, lhs))
        return cls(LE, seq)

    @classmethod
    def Gt(cls, *seq):
        if len(seq)==2:
            lhs, rhs = seq
            if lhs.head is NUMBER and rhs.head is NUMBER:
                return cls.true if lhs.data > rhs.data else cls.false
            if not lhs < rhs:
                return cls(LT, (rhs, lhs))
        return cls(GT, seq)

    @classmethod
    def Ge(cls, *seq):
        if len(seq)==2:
            lhs, rhs = seq
            if lhs.head is NUMBER and rhs.head is NUMBER:
                return cls.true if lhs.data >= rhs.data else cls.false
            if not lhs < rhs:
                return cls(LE, (rhs, lhs))
        return cls(GE, seq)

    @classmethod
    def Eq(cls, *seq):
        if len(seq)==2:
            lhs, rhs = seq
            if lhs==rhs:
                return cls.true
            if lhs.head is NUMBER and rhs.head is NUMBER:
                return cls.false
        return cls(EQ, seq)

    @classmethod
    def Ne(cls, *seq):
        if len(seq)==2:
            lhs, rhs = seq
            if lhs==rhs:
                return cls.false
            if lhs.head is NUMBER and rhs.head is NUMBER:
                return cls.true
        return cls(NE, seq)

    @classmethod
    def Element(cls, element, container):
        if hasattr(container, 'contains'):
            return container.contains(element)
        return cls(IN, (element, container))

    @classmethod
    def Contains(cls, container, element):
        return cls.IsElement(element, container)

    #XXX: IsSubset

    @classmethod
    def Or(cls, *seq):
        s = set([])
        for a in seq:
            assert type(a) is cls,`type(a)`
            head, data = a.pair
            if head is NUMBER:
                if data:
                    return cls.true
            elif head is OR:
                s.update(data)
            elif head is NOT:
                if data in s:
                    return cls.true
                s.add(a)
            else:
                if cls(NOT, a) in s:
                    return cls.true
                s.add(a)
        if not s:
            return cls.false
        if len(s)==1:
            return s.pop()
        return cls(OR, frozenset(s))

    @classmethod
    def And(cls, *seq):
        s = set([])
        for a in seq:
            assert type(a) is cls,`type(a)`
            head, data = a.pair
            if head is NUMBER:
                if not data:
                    return cls.false
            elif head is AND:
                s.update(data)
            elif head is NOT:
                if data in s:
                    return cls.false
                else:
                    s.add(a)
            else:
                if cls(NOT,a) in s:
                    return cls.false
                s.add(a)
        if not s:
            return cls.true
        if len(s)==1:
            return s.pop()
        return cls(AND, frozenset(s))

    @classmethod
    def Not(cls, obj):
        assert type(obj) is cls,`type(obj)`
        head, data = obj.pair
        if head is NUMBER: return cls(NUMBER, not data)
        if head is NOT: return data
        if head is NE: return cls(EQ, data)
        if head is EQ: return cls(NE, data)
        if head is LT: return cls(GE, data)
        if head is GT: return cls(LE, data)
        if head is LE: return cls(GT, data)
        if head is GE: return cls(LT, data)
        return cls(NOT, obj)

    #XXX: Xor, Implies, Equiv

Logic.one = Logic.true = Logic(NUMBER, True)
Logic.zero = Logic.false = Logic(NUMBER, False)
classes.Logic = Logic
