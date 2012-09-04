

from .base import Head, heads_precedence

class IntegralIndefiniteHead(Head):

    def __repr__(self): return 'INTEGRAL_INDEFINITE'

    def data_to_str_and_precedence(self, cls, (integrand, symbol)):
        o_p = heads_precedence.APPLY
        s, s_p = integrand.head.data_to_str_and_precedence(cls, integrand.data)
        return 'Integral(%s, %s)' % (s, symbol), o_p


INTEGRAL_INDEFINITE = IntegralIndefiniteHead()
    
