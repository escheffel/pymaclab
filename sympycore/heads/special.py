
from .base import AtomicHead, heads_precedence, Expr

import re
_is_atomic = re.compile(r'\A\w+\Z').match

class SpecialHead(AtomicHead):
    """
    SpecialHead is head for special objects like None, NotImplemented,
    Ellipsis, etc.  Data can be any Python object.
    """

    def __repr__(self): return 'SPECIAL'

    def data_to_str_and_precedence(self, cls, data):
        if isinstance(data, Expr):
            h, d = data.pair
            return h.data_to_str_and_precedence(cls, d)
        if data is Ellipsis:
            return '...', heads_precedence.SYMBOL

        s = str(data)
        if _is_atomic(s):
            return s, heads_precedence.SYMBOL
        return s, 0.0 # force parenthesis

SPECIAL = SpecialHead()
