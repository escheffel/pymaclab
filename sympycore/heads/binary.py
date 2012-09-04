
__all__ = ['INVERT', 'BOR', 'BXOR', 'BAND', 'LSHIFT', 'RSHIFT']

from .base import Head, UnaryHead, NaryHead

class InvertHead(UnaryHead):
    """
    InvertHead is unary invert operation, data is expression operand.
    """
    op_mth = '__invert__'
    op_symbol = '~'
    def __repr__(self): return 'INVERT'

class BorHead(NaryHead):
    """
    BorHead is a head of n-ary bitwise inclusive (normal) OR operation,
    data is n-tuple of expression operands.
    """
    op_mth = '__or__'
    op_symbol = '|'
    def __repr__(self): return 'BOR'

class BxorHead(NaryHead):
    """
    BxorHead is a head of n-ary bitwise exclusive OR operation,
    data is n-tuple of expression operands.
    """
    op_mth = '__xor__'
    op_symbol = '^'
    def __repr__(self): return 'BXOR'

class BandHead(NaryHead):
    """
    BandHead is a head of n-ary bitwise AND operation,
    data is n-tuple of expression operands.
    """
    op_mth = '__and__'
    op_symbol = '&'
    def __repr__(self): return 'BAND'

class LshiftHead(NaryHead):
    """
    LshiftHead is a head of n-ary bitwise shift left operation,
    data is n-tuple of expression operands.
    """
    op_mth = '__lshift__'
    op_symbol = '<<'
    def __repr__(self): return 'LSHIFT'

class RshiftHead(NaryHead):
    """
    RshiftHead is a head of n-ary bitwise shift right operation,
    data is n-tuple of expression operands.
    """
    op_mth = '__rshift__'
    op_symbol = '>>'
    def __repr__(self): return 'RSHIFT'

INVERT = InvertHead()
BOR = BorHead()
BXOR = BxorHead()
BAND = BandHead()
LSHIFT = LshiftHead()
RSHIFT = RshiftHead()
