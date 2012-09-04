
from __future__ import with_statement

__all__ = ['EQ', 'NE', 'LT', 'LE', 'GT', 'GE']

from ..core import init_module, classes, SymbolicEquality
init_module.import_heads()
init_module.import_numbers()

from .base import Head, BinaryHead, heads_precedence, Expr

class RelationalHead(BinaryHead):

    def walk(self, func, cls, data, target):
        lhs, rhs = data
        lhs1 = lhs.head.walk(func, classes.Calculus, lhs.data, lhs)
        rhs1 = rhs.head.walk(func, classes.Calculus, rhs.data, rhs)
        if lhs1 is lhs and rhs1 is rhs:
            return func(cls, data, target)
        r = self.new(cls, (lhs1, rhs1))
        return func(cls, r.head, r.data, r)

class EqHead(RelationalHead):

    op_symbol = '=='
    op_mth = '__eq__'

    def new(self, cls, (lhs, rhs)):
        if lhs.head is NUMBER and rhs.head is NUMBER:
            return cls(lhs.data == rhs.data)
        return cls(self, (lhs, rhs))

    def reevaluate(self, cls, (lhs, rhs)):
        with SymbolicEquality(type(lhs), type(rhs)):
            return lhs == rhs

    def __repr__(self): return 'EQ'

    def nonzero(self, cls, (lhs, rhs)):
        if isinstance(lhs, Expr):
            lhs = lhs.as_lowlevel()
        if isinstance(rhs, Expr):
            rhs = rhs.as_lowlevel()
        if isinstance(lhs, numbertypes) and isinstance(rhs, numbertypes):
            return lhs == rhs
        return True

class NeHead(RelationalHead):

    op_symbol = '!='
    op_mth = '__ne__'

    def new(self, cls, (lhs, rhs)):
        if lhs.head is NUMBER and rhs.head is NUMBER:
            return cls(lhs.data != rhs.data)
        return cls(self, (lhs, rhs))

    def reevaluate(self, cls, (lhs, rhs)):
        with SymbolicEquality(type(lhs), type(rhs)):
            return lhs != rhs

    def __repr__(self): return 'NE'

    def nonzero(self, cls, (lhs, rhs)):
        if isinstance(lhs, Expr):
            lhs = lhs.as_lowlevel()
        if isinstance(rhs, Expr):
            rhs = rhs.as_lowlevel()
        if isinstance(lhs, numbertypes) and isinstance(rhs, numbertypes):
            return lhs != rhs
        return True
    
class LtHead(RelationalHead):

    op_symbol = '<'
    op_mth = '__lt__'
    def __repr__(self): return 'LT'

    def new(self, cls, (lhs, rhs), evaluate=True):
        if lhs.head is NUMBER and rhs.head is NUMBER:
            return cls(lhs.data < rhs.data)
        return cls(self, (lhs, rhs))

    def reevaluate(self, cls, (lhs, rhs)):
        return lhs < rhs

    def nonzero(self, cls, (lhs, rhs)):
        if isinstance(lhs, Expr):
            lhs = lhs.as_lowlevel()
        if isinstance(rhs, Expr):
            rhs = rhs.as_lowlevel()
        if isinstance(lhs, numbertypes) and isinstance(rhs, numbertypes):
            return lhs < rhs
        return True

class LeHead(RelationalHead):

    op_symbol = '<='
    op_mth = '__le__'

    def new(self, cls, (lhs, rhs)):
        if lhs.head is NUMBER and rhs.head is NUMBER:
            return cls(lhs.data <= rhs.data)
        return cls(self, (lhs, rhs))

    def reevaluate(self, cls, (lhs, rhs)):
        return lhs <= rhs
    
    def __repr__(self): return 'LE'

    def nonzero(self, cls, (lhs, rhs)):
        if isinstance(lhs, Expr):
            lhs = lhs.as_lowlevel()
        if isinstance(rhs, Expr):
            rhs = rhs.as_lowlevel()
        if isinstance(lhs, numbertypes) and isinstance(rhs, numbertypes):
            return lhs <= rhs
        return True

class GtHead(RelationalHead):

    op_symbol = '>'
    op_mth = '__gt__'

    def new(self, cls, (lhs, rhs)):
        if lhs.head is NUMBER and rhs.head is NUMBER:
            return cls(lhs.data > rhs.data)
        return cls(self, (lhs, rhs))


    def reevaluate(self, cls, (lhs, rhs)):
        return lhs > rhs
    
    def __repr__(self): return 'GT'

    def nonzero(self, cls, (lhs, rhs)):
        if isinstance(lhs, Expr):
            lhs = lhs.as_lowlevel()
        if isinstance(rhs, Expr):
            rhs = rhs.as_lowlevel()
        if isinstance(lhs, numbertypes) and isinstance(rhs, numbertypes):
            return lhs > rhs
        return True
    
class GeHead(RelationalHead):

    op_symbol = '>='
    op_mth = '__ge__'

    def new(self, cls, (lhs, rhs)):
        if lhs.head is NUMBER and rhs.head is NUMBER:
            return cls(lhs.data >= rhs.data)
        return cls(self, (lhs, rhs))

    def reevaluate(self, cls, (lhs, rhs)):
        return lhs >= rhs
    
    def __repr__(self): return 'GE'

    def nonzero(self, cls, (lhs, rhs)):
        if isinstance(lhs, Expr):
            lhs = lhs.as_lowlevel()
        if isinstance(rhs, Expr):
            rhs = rhs.as_lowlevel()
        if isinstance(lhs, numbertypes) and isinstance(rhs, numbertypes):
            return lhs >= rhs
        return True

EQ = EqHead()
NE = NeHead()
LT = LtHead()
LE = LeHead()
GT = GtHead()
GE = GeHead()
