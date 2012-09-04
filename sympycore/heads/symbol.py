
__all__ = ['SYMBOL']

import re

from ..core import init_module
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

from .base import AtomicHead, heads_precedence, Expr, Pair, ArithmeticHead, NotImplementedHeadMethod

_is_atomic = re.compile(r'\A\w+\Z').match

class SymbolHead(AtomicHead):
    """
    SymbolHead is a head for symbols, data can be any Python object.
    """

    def new(self, cls, data, evaluate=True):
        return cls(self, data)

    def reevaluate(self, cls, data):
        return cls(self, data)

    def is_data_ok(self, cls, data):
        return

    def __repr__(self): return 'SYMBOL'

    def data_to_str_and_precedence(self, cls, data):
        if isinstance(data, Expr):
            h, d = data.pair
            return h.data_to_str_and_precedence(cls, d)
        s = str(data)
        if _is_atomic(s):
            return s, heads_precedence.SYMBOL
        return s, 0.0 # force parenthesis

    def to_EXP_COEFF_DICT(self, cls, data, expr, variables = None):
        variables = EXP_COEFF_DICT.combine_variables(data, variables)
        exp = EXP_COEFF_DICT.make_exponent(data, variables)
        assert len(exp)==len(variables), `exp, variables, i, data`
        return cls(EXP_COEFF_DICT, Pair(variables, {exp:1}))

    def commutative_mul_number(self, cls, lhs, rhs):
        return term_coeff_new(cls, (lhs, rhs))

    non_commutative_mul_number = commutative_mul_number

    def pow(self, cls, base, exp):
        if type(exp) is cls and exp.head is NUMBER:
            exp = exp.data
        return pow_new(cls, (base, exp))

    def pow_number(self, cls, base, exp):
        if exp==1: return base
        if exp==0: return cls(NUMBER, 1)
        return cls(POW, (base, exp))

    def expand(self, cls, expr):
        return expr

    def expand_intpow(self, cls, expr, intexp):
        if intexp==0:
            return cls(NUMBER, 1)
        return cls(POW, (expr, intexp))

    def diff(self, cls, data, expr, symbol, order, cache={}):
        if order==0:
            return expr
        if data == symbol:
            assert order>0,`order`
            return cls(NUMBER, int(order==1))
        return cls(NUMBER, 0)

    def fdiff(self, cls, data, expr, argument_index, order):
        vcls = cls.get_value_algebra()
        dcls = cls.get_differential_algebra()
        d = dcls(FDIFF, vcls(NUMBER, argument_index))**order
        return cls(APPLY, (d, (expr,)))

    def apply(self, cls, data, func, args):
        return cls(APPLY, (func, args))

    def integrate_indefinite(self, cls, data, expr, x):
        if data==x:
            return cls(TERM_COEFF, (cls(POW, (expr, 2)), mpq((1,2))))
        return cls(BASE_EXP_DICT, {expr:1, cls(SYMBOL, x):1})

    def integrate_definite(self, cls, data, expr, x, a, b):
        if data==x:
            return (b**2-a**2)/2
        return expr*(b-a)

    def to_TERM_COEFF_DICT(self, Algebra, data, expr):
        return expr

    def to_MUL(self, Algebra, data, expr):
        return expr

    def to_ADD(self, cls, data, expr):
        return expr

    def algebra_pos(self, Algebra, expr):
        return expr

    def algebra_neg(self, Algebra, expr):
        if Algebra.algebra_options.get('is_additive_group_commutative'):
            if Algebra.algebra_options.get('evaluate_addition'):
                return Algebra(TERM_COEFF, (expr, -1))
        return Algebra(NEG, expr)

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_addition'):
            if not rhs:
                return lhs            
            if Algebra.algebra_options.get('is_additive_group_commutative'):
                return Algebra(TERM_COEFF_DICT, {Algebra(NUMBER, 1): rhs, lhs:1})
        return Algebra(ADD, [lhs, Algebra(NUMBER, rhs)])

    def algebra_add(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if Algebra.algebra_options.get('is_additive_group_commutative'):
            ldata = lhs.data
            if Algebra.algebra_options.get('evaluate_addition'):
                if rhead is ADD or rhead is EXP_COEFF_DICT or rhead is MUL or rhead is NEG:
                    rhead, rdata = rhs.to(TERM_COEFF_DICT).pair
                if rhead is SYMBOL:
                    if ldata == rdata:
                        return Algebra(TERM_COEFF, (lhs, 2))
                    return Algebra(TERM_COEFF_DICT, {lhs: 1, rhs: 1})
                if rhead is NUMBER:
                    if rdata:
                        return Algebra(TERM_COEFF_DICT, {Algebra(NUMBER, 1): rdata, lhs:1})
                    return lhs
                if rhead is TERM_COEFF:
                    term, coeff = rdata
                    if term==lhs:
                        return term_coeff_new(Algebra, (term, coeff+1))
                    return Algebra(TERM_COEFF_DICT, {term: coeff, lhs:1})
                if rhead is TERM_COEFF_DICT:
                    d = rdata.copy()
                    term_coeff_dict_add_item(Algebra, d, lhs, 1)
                    return term_coeff_dict_new(Algebra, d)
                if rhead is POW or rhead is BASE_EXP_DICT:
                    return Algebra(TERM_COEFF_DICT, {lhs:1, rhs:1})
            else:
                if rhead is TERM_COEFF_DICT or rhead is EXP_COEFF_DICT:
                    rhs = rhs.to(ADD)
                    rhead, rdata = rhs.pair
                if rhead is ADD:
                    data = [lhs] + rdata
                else:
                    data = [lhs, rhs]
                return add_new(Algebra, [lhs, rhs])
            return super(type(self), self).algebra_add(Algebra, lhs, rhs, inplace)
        else:
            if rhead is TERM_COEFF_DICT or rhead is EXP_COEFF_DICT:
                rhs = rhs,to(ADD)
                rhead, rdata = rhs.pair
            if rhead is ADD:
                data = [lhs] + rdata
            else:
                data = [lhs, rhs]
            if Algebra.algebra_options.get('evaluate_addition'):
                ADD.combine_add_list(Algebra, data)
            return add_new(Algebra, data)

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_addition'):
            return term_coeff_new(Algebra, (lhs, rhs))
        return mul_new(Algebra, [lhs, Algebra(NUMBER, rhs)])

    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('is_additive_group'):
            if Algebra.algebra_options.get('evaluate_addition'):
                rhead, rdata = rhs.pair
                if rhead is NUMBER:
                    return term_coeff_new(Algebra, (lhs, rdata))
            else:
                return mul_new(Algebra, [lhs, rhs])
        if Algebra.algebra_options.get('is_multiplicative_group'):
            if Algebra.algebra_options.get('evaluate_multiplication'):
                rhead, rdata = rhs.pair
                if Algebra.algebra_options.get('is_multiplicative_group_commutative'):
                    if rhead is MUL or EXP_COEFF_DICT:
                        rhs = rhead.to_BASE_EXP_DICT(Algebra, rdata, rhs)
                        rhead, rdata = rhs.pair
                    if rhead is NUMBER:
                        return term_coeff_new(Algebra, (lhs, rdata))
                    if rhead is BASE_EXP_DICT:
                        data = rdata.copy()
                        base_exp_dict_add_item(Algebra, data, lhs, 1)
                        return base_exp_dict_new(Algebra, data)
                    if lhs==rhs:
                        return Algebra(POW, (lhs, 2))
                    return Algebra(BASE_EXP_DICT, {lhs:1, rhs:1})
                else:
                    if rhead is MUL:
                        data = [lhs] + rdata
                    else:
                        data = [lhs, rhs]
                    coeff = MUL.combine_mul_list(Algebra, data)
                    assert coeff is 1,`coeff`
                    return mul_new(Algebra, data)
            else:
                pass
        return super(type(self), self).algebra_mul(Algebra, lhs, rhs, inplace)

    def algebra_pow_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('is_multiplicative_group'):
            if Algebra.algebra_options.get('evaluate_multiplication'):
                if Algebra.algebra_options.get('is_multiplicative_group_commutative'):
                    return pow_new(Algebra, (lhs, rhs))
                else:
                    return pow_new(Algebra, (lhs, rhs))
            else:
                return cls(POW, (lhs, rhs))
        return super(type(self), self).algebra_pow_number(Algebra, lhs, rhs, inplace)

SYMBOL = SymbolHead()
