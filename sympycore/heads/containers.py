
__all__ = ['TUPLE', 'LIST', 'DICT']

from .base import Head, heads_precedence

def init_module(m):
    from .base import heads
    for n,h in heads.iterNameValue(): setattr(m, n, h)

class ContainerHead(Head):

    def data_to_str_and_precedence(self, cls, seq):
        c_p = heads_precedence.COMMA
        s_p = getattr(heads_precedence, repr(self))
        l = []
        for item in seq:
            i, i_p = item.head.data_to_str_and_precedence(cls, item.data)
            if i_p < c_p: i = '('+i+')'
            l.append(i)
        s1, s2 = self.parenthesis
        if self is TUPLE and len(l)==1:
            return s1 + l[0] + ',' + s2, s_p
        return s1 + ', '.join(l) + s2, s_p

class TupleHead(ContainerHead):
    """
    TupleHead represents n-tuple,
    data is n-tuple of expressions.
    """
    parenthesis = '()'

    def __repr__(self): return 'TUPLE'

class ListHead(ContainerHead):
    """
    ListHead represents n-list,
    data is n-tuple of expressions.
    """
    parenthesis = '[]'    
    def __repr__(self): return 'LIST'

class DictHead(ContainerHead):
    """
    DictHead represents n-dict,
    data is n-tuple of expression pairs.
    """

    def data_to_str_and_precedence(self, cls, seq2):
        c_p = heads_precedence.COMMA
        colon_p = heads_precedence.COLON
        s_p = getattr(heads_precedence, repr(self))
        l = []
        for key, value in seq2:
            k, k_p = key.head.data_to_str_and_precedence(cls, key.data)
            if k_p < colon_p: k = '('+k+')'
            v, v_p = value.head.data_to_str_and_precedence(cls, value.data)
            if v_p < colon_p: v = '('+v+')'
            l.append(k + ':' + v)
        return '{' + ', '.join(l) + '}', heads_precedence.DICT

    def __repr__(self): return 'DICT'

TUPLE = TupleHead()
LIST = ListHead()
DICT = DictHead()
