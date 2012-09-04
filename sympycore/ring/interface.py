
__all__ = ['RingInterface']

from ..core import init_module

@init_module
def _init(m):
    from ..core import heads
    for n,h in heads.iterNameValue(): setattr(m, n, h)

class RingInterface:

    zero = 0
    one = 1

    @classmethod
    def Add(cls, *seq):
        """ Construct a symbolic expression with addition operands
        list representation.
        """
        return ADD.reevaluate(cls, seq)

    @classmethod
    def Sub(cls, *seq):
        return SUB.reevaluate(cls, seq)

    @classmethod
    def Mul(cls, *seq):
        """ Construct a symbolic expression with multiplication
        operands list representation.
        """
        return MUL.reevaluate(cls, seq)

    @classmethod
    def Div(cls, *seq):
        return DIV.reevaluate(cls, seq)

    @classmethod
    def Pow(cls, base, exp):
        return POW.new(cls, (base, exp))

    @classmethod
    def Terms(cls, *seq):
        """ Construct a symbolic expression with term-coefficient
        representation.
        """
        return TERM_COEFF_DICT.reevaluate(cls, dict(seq))

    @classmethod
    def Factors(cls, *seq):
        """ Construct a symbolic expression with base-exponent
        representation.
        """
        return BASE_EXP_DICT.reevaluate(cls, dict(seq))

    @classmethod
    def Polynom(cls, *seq, **kws):
        """ Construct a symbolic expression with exponent-coefficient
        representation:

          Polynom(x,y) -> Algebra(EXP_COEFF_DICT, Pair((x,y), {(1,0):1, (0,1):1}))
          Polynom(x,y, variables=(x,y,z)) -> Algebra(EXP_COEFF_DICT, Pair((x,y,z), {(1,0,0):1, (0,1,0):1}))
        """
        r = cls.Add(*seq)
        return r.head.to_EXP_COEFF_DICT(cls, r.data, r, kws.get('variables'))
