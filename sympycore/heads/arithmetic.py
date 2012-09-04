
#obsolete, each head should be defined in a separate file

__all__ = ['FLOORDIV', 'MOD']


from .base import NaryHead

class ModHead(NaryHead):

    """
    ModHead represents module n-ary operation,
    data is a n-tuple of expression operands.
    """
    op_mth = '__mod__'
    op_rmth = '__rmod__'
    op_symbol = '%'
    def __repr__(self): return 'MOD'

class FloordivHead(NaryHead):
    """
    FloordivHead represents floor-division n-ary operation,
    data is a n-tuple of expression operands.
    """
    op_mth = '__floordiv__'
    op_rmth = '__rfloordiv__'
    op_symbol = '//'    

    def __repr__(self): return 'FLOORDIV'



MOD = ModHead()

FLOORDIV = FloordivHead()

