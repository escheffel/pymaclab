
__all__ = ['AND', 'OR', 'NOT', 'IN', 'NOTIN', 'IS', 'ISNOT']

from .base import Head, UnaryHead, BinaryHead, NaryHead
from ..core import init_module, classes
init_module.import_heads()

class NotHead(UnaryHead):
    """
    NotHead represents unary boolean NOT operation,
    data is an expression operand.
    """
    op_symbol = 'not '
    def __repr__(self): return 'NOT'

    def new(self, cls, operand):
        assert type(operand) is cls,`type(operand), cls`
        head, data = operand.pair
        if head is NUMBER: return cls(NUMBER, not data)
        if head is NOT: return data
        if head is NE: return cls(EQ, data)
        if head is EQ: return cls(NE, data)
        if head is LT: return cls(GE, data)
        if head is GT: return cls(LE, data)
        if head is LE: return cls(GT, data)
        if head is GE: return cls(LT, data)
        return cls(NOT, operand)

    def reevaluate(self, cls, operand):
        return self.new(cls, operand)

class AndHead(NaryHead):
    """
    AndHead represents n-ary boolean AND operation,
    data is a n-tuple of expression operands.
    """
    op_symbol = ' and '
    def __repr__(self): return 'AND'

    def new(self, cls, operands):
        s = set([])
        for a in operands:
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

class OrHead(NaryHead):
    """
    AndHead represents n-ary boolean OR operation,
    data is a n-tuple of expression operands.
    """
    op_symbol = ' or '
    def __repr__(self): return 'OR'

    def new(self, cls, operands):
        s = set([])
        for a in operands:
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

class IsHead(BinaryHead):
    """
    IsHead represents binary boolean IS operation,
    data is an expression operand.
    """
    op_symbol = ' is '
    def __repr__(self): return 'IS'

class IsnotHead(BinaryHead):
    """
    IsHead represents binary boolean IS NOT operation,
    data is an expression operand.
    """
    op_symbol = ' is not '
    def __repr__(self): return 'ISNOT'

class InHead(BinaryHead):
    """
    InHead represents binary boolean IN operation,
    data is an expression operand.
    """
    op_symbol = ' in '

    def new(self, cls, (lhs, rhs)):
        return cls(self, (lhs, rhs))
    
    def __repr__(self): return 'IN'

    def walk(self, func, cls, data, target):
        lhs, rhs = data
        rcls = cls.get_operand_algebra(self, 1)
        rhs1 = rhs.head.walk(func, rcls, rhs.data, rhs)
        lcls = rhs1.get_element_algebra()
        lhs1 = lhs.head.walk(func, lcls, lhs.data, lhs)
        if lhs1 is lhs and rhs1 is rhs:
            return func(cls, data, target)
        r = self.new(cls, (lhs1, rhs1))
        return func(cls, r.head, r.data, r)

class NotinHead(BinaryHead):
    """
    NotinHead represents binary boolean NOT IN operation,
    data is an expression operand.
    """
    op_symbol = ' not in '
    def __repr__(self): return 'NOTIN'

NOT = NotHead()
OR = OrHead()
AND = AndHead()
IS = IsHead()
ISNOT = IsnotHead()
IN = InHead()
NOTIN = NotinHead()
