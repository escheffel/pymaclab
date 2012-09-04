
__all__ = ['DIV']

from .base import NaryHead
from ..core import init_module

init_module.import_heads()

class DivHead(NaryHead):

    """
    DivHead represents division n-ary operation,
    data is a n-tuple of expression operands.
    """
    op_mth = '__div__'
    op_rmth = '__rdiv__'
    op_symbol = '/'


    def __repr__(self): return 'DIV'

    def new(self, cls, operands, evaluate=True):
        m = len(operands)
        if m==0:
            return cls(NUMBER, 1)
        if m==1:
            return operands[0]
        if not evaluate:
            return cls(self, operands)
        return operands[0] / MUL.new(cls, reversed(operands[1:]))

    def reevaluate(self, cls, operands):
        r = operands[0] if operands else cls(NUMBER, 1)
        for op in operands[1:]:
            r /= op
        return r

    def term_coeff(self, cls, expr):
        return expr, 1

    def commutative_mul(self, cls, lhs, rhs):
        d = {}
        imul = BASE_EXP_DICT.inplace_commutative_data_mul
        for i, op in enumerate(lhs.data):
            if i:
                imul(cls, d, 1/op)
            else:
                imul(cls, d, op)
        imul(cls, d, rhs)
        return BASE_EXP_DICT.new(cls, d)

    def pow(self, cls, lhs, rhs):
        data = lhs.data
        l = [data[0]] + [1/d for d in data[1:]]
        return cls(MUL, l) ** rhs

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        if not rhs:
            return lhs
        return Algebra(ADD, [lhs, Algebra(NUMBER, rhs)])
        
DIV = DivHead()
