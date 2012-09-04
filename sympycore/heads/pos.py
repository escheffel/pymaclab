
__all__ = ['POS']

from .base import ArithmeticHead

def init_module(m):
    from .base import heads
    for n,h in heads.iterNameValue(): setattr(m, n, h)

class PosHead(ArithmeticHead):

    """
    PosHead represents positive unary operation where operand (data)
    is an expression.
    """

    op_mth = '__pos__'
        
    def __repr__(self): return 'POS'

    def data_to_str_and_precedence(self, cls, expr):
        return expr.head.data_to_str_and_precedence(cls, expr.data)

    def to_lowlevel(self, cls, data, pair):
        return data

POS = PosHead()
