
__all__ = ['NEG']

from .base import heads_precedence, ArithmeticHead
from ..core import init_module

init_module.import_heads()
init_module.import_lowlevel_operations()
init_module.import_numbers()

class NegHead(ArithmeticHead):

    """
    NegHead represents negative unary operation where operand (data)
    is an expression.
    """

    op_mth = '__neg__'

    def is_data_ok(self, cls, data):
        if not isinstance(data, cls):
            return '%s data part must be %s instance but got %s' % (self, cls, type(data))
        
    def __repr__(self): return 'NEG'

    def new(self, cls, expr, evaluate=True):
        if not evaluate:
            return cls(self, expr)
        h, d = expr.pair
        if h is NUMBER:
            return cls(NUMBER, -d)
        if h is TERM_COEFF:
            t, c = d
            return TERM_COEFF.new(cls, (t, -c))
        return cls(NEG, expr)

    def reevaluate(self, cls, expr):
        return -expr

    def data_to_str_and_precedence(self, cls, expr):
        if cls.algebra_options.get('evaluate_addition'):
            if expr.head is NEG:
                expr = expr.data
                return expr.head.data_to_str_and_precedence(cls, expr.data)
        s, s_p = expr.head.data_to_str_and_precedence(cls, expr.data)
        neg_p = heads_precedence.NEG
        if s_p < neg_p:
            return '-(' + s + ')', neg_p
        return '-' + s, neg_p

    def to_lowlevel(self, cls, data, pair):
        if isinstance(data, numbertypes):
            return -data
        if data.head is NUMBER:
            return -data.data
        return cls(TERM_COEFF, (data, -1))

    def term_coeff(self, cls, expr):
        e = expr.data
        t, c = e.head.term_coeff(cls, e)
        return t, -c

    def scan(self, proc, cls, expr, target):
        expr.head.scan(proc, cls, expr.data, target)
        proc(cls, self, expr, target)

    def walk(self, func, cls, operand, target):
        operand1 = operand.head.walk(func, cls, operand.data, operand)
        if operand1 is operand:
            return func(cls, self, operand, target)
        r = self.new(cls, operand1)
        return func(cls, r.head, r.data, r)

    def to_TERM_COEFF_DICT(self, Algebra, data, expr):
        return -data.head.to_TERM_COEFF_DICT(Algebra, data.data, data)

    def to_ADD(self, Algebra, data, expr):
        return -data.head.to_ADD(Algebra, data.data, data)

    def algebra_pos(self, Algebra, expr):
        return expr

    def algebra_neg(self, Algebra, expr):
        if Algebra.algebra_options.get('evaluate_addition'):
            return expr.data
        return Algebra(NEG, expr)

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        return self.algebra_add(Algebra, lhs, Algebra(NUMBER, rhs), inplace)

    def algebra_add(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if rhead is ADD:
            data = [lhs] + rdata
        elif rhead is TERM_COEFF_DICT or rhead is EXP_COEFF_DICT:
            data = [lhs] + rhs.to(ADD).data
        else:
            data = [lhs, rhs]
        if Algebra.algebra_options.get('evaluate_addition'):
            ADD.combine_add_list(Algebra, data)
        return add_new(Algebra, data)

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('is_additive_group_commutative'):
            term, coeff = lhs.head.term_coeff(Algebra, lhs)
            return term_coeff_new(Algebra, (term, coeff * rhs))
        else:
            if Algebra.algebra_options.get('evaluate_addition'):
                if rhs == 0:
                    return Algebra(NUMBER, 0)
                term, coeff = lhs.head.term_coeff(Algebra, lhs)
                return term_coeff_new(Algebra, (term, coeff * rhs))
            return mul_new(Algebra, [lhs, Algebra(NUMBER, rhs)])

    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        ldata = lhs.data
        if Algebra.algebra_options.get('is_additive_group_commutative'):
            return super(type(self), self).algebra_mul(Algebra, lhs, rhs, inplace)
        else:
            if Algebra.algebra_options.get('evaluate_addition'):
                rhead, rdata = rhs.pair
                if rhead is NUMBER:
                    return ldata.head.algebra_mul_number(Algebra, ldata, -rdata, inplace)
                return super(type(self), self).algebra_mul(Algebra, lhs, rhs, inplace)
            return mul_new(Algebra, [lhs, rhs])
    
NEG = NegHead()
