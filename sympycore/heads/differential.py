
from ..core import init_module
from .base import AtomicHead, Expr, heads_precedence

init_module.import_heads()
init_module.import_lowlevel_operations()

class DifferentialHead(AtomicHead):
    
    def base_exp(self, cls, expr):
        return expr, 1

    def term_coeff(self, cls, expr):
        return expr, 1

    def commutative_mul_number(self, cls, lhs, rhs):
        return term_coeff_new(cls, (lhs, rhs))

    def add(self, cls, lhs, rhs):
        h,d = rhs.pair
        if h is NUMBER:
            if d==0: return lhs
            return cls(TERM_COEFF_DICT, {lhs:1, cls(NUMBER,1):d})
        elif h is self:
            if lhs.data == d:
                return cls(TERM_COEFF, (lhs, 2))
            return cls(TERM_COEFF_DICT, {lhs:1, rhs:1})
        elif h is TERM_COEFF:
            t,c = d
            if lhs==t:
                return term_coeff_new(cls, (t, c+1))
            return cls(TERM_COEFF_DICT, {t:c, lhs:1})
        elif h is TERM_COEFF_DICT:
            data = d.copy()
            dict_add_item(cls, data, lhs, 1)
            return term_coeff_dict_new(cls, data)
        raise NotImplementedError(`self, h`)

    def commutative_mul(self, cls, lhs, rhs):
        h, d = rhs.pair
        if h is NUMBER:
            return cls.commutative_mul_number(cls, lhs, d)
        if h is self:
            if lhs.data==d:
                return cls(POW, (lhs, 2))
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        if h is POW:
            base, exp = d
            if base==lhs:
                return pow_new(cls, (base, exp + 1))
            return cls(BASE_EXP_DICT, {base:exp, lhs:1})
            
        raise NotImplementedError(`self, h`)

    def pow(self, cls, base, exp):
        return pow_new(cls, (base, exp))

    pow_number = pow

class DiffHead(DifferentialHead):

    """
    FunctionAlgebra(DIFF, x) - differential with respect to x
    x is symbol
    """
    
    def __repr__(self): return 'DIFF'

    def is_data_ok(self, cls, data):
        if isinstance(data, Expr):
            if data.head is SYMBOL: return
            return 'data must be with SYMBOL head but got %r' % (data.head)
        return 'data must be symbol but got %s' % (type(data))

    def data_to_str_and_precedence(self, cls, data):
        return SUBSCRIPT.data_to_str_and_precedence(cls, (Expr(SYMBOL, 'D'), (data,)))

    def diff_apply(self, cls, data, diff, expr):
        return expr.head.diff(type(expr), expr.data, expr, data.data, 1)

class FDiffHead(DifferentialHead):

    """
    FunctionAlgebra(FDIFF, x) - differential with respect to x-th argument
    x is number or symbol
    """
    
    def __repr__(self): return 'FDIFF'

    def is_data_ok(self, cls, data):
        if isinstance(data, Expr):
            if data.head is SYMBOL or data.head is NUMBER: return
            return 'data must be with SYMBOL|NUMBER head but got %r' % (data.head)
        return 'data must be symbol or number but got %s' % (type(data))

    def data_to_str_and_precedence(self, cls, data):
        return SUBSCRIPT.data_to_str_and_precedence(cls, (Expr(SYMBOL, 'FD'), (data,)))

    def diff_apply(self, cls, data, diff, expr):
        return expr.head.fdiff(type(expr), expr.data, expr, data.data, 1)


DIFF = DiffHead()
FDIFF = FDiffHead()
