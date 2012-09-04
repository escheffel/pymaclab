
__all__ = ['TERM_COEFF']

from .base import heads_precedence, ArithmeticHead, Expr

from ..core import init_module
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

class TermCoeff(ArithmeticHead):
    """ Expr(TERM_COEFF, (term, coeff)) represents term*coeff
    where term is symbolic expression and coeff is a number or
    symbolic expression.
    """

    def is_data_ok(self, cls, data, allow_number_term=True):
        if type(data) is tuple and len(data)==2:
            term, coeff = data
            if isinstance(term, cls):
                if isinstance(coeff, numbertypes) or type(coeff) is not cls:
                    if term.head is self or (not allow_number_term and term.head is NUMBER):
                        return 'data[0].head cannot be %s' %(term.head)
                    return
                elif isinstance(coeff, cls):
                    if coeff.head is NUMBER:
                        if not isinstance(coeff.data, numbertypes): #pragma: no cover
                            return 'data[1].data must be %s instance for NUMBER head but got %s instance' % (numbertypes, type(coeff.data)) #pragma: no cover
                else:
                    return 'data[1] must be %s instance but got %s instance' % ((cls, numbertypes), type(coeff)) #pragma: no cover
            else:
                return 'data[0] must be %s instance but got %s instance' % (cls, type(term)) #pragma: no cover
        else:
            return 'data must be 2-tuple' #pragma: no cover
        return

    def __repr__(self): return 'TERM_COEFF'

    def new(self, cls, data):
        return term_coeff_new(cls, data)

    def reevaluate(self, cls, (term, coeff)):
        return term * coeff

    def to_EXP_COEFF_DICT(self, cls, (term, coeff), expr, variables=None):
        return term.head.to_EXP_COEFF_DICT(cls, term.data, term, variables) * coeff

    def to_MUL(self, Algebra, (term, coeff), expr):
        term = term.head.to_MUL(Algebra, term.data, term)
        if term.head is MUL:
            data = [Algebra(NUMBER, coeff)] + term.data
        else:
            data = [Algebra(NUMBER, coeff), term]
        return Algebra(MUL, data)

    def to_ADD(self, Algebra, (term, coeff), expr):
        return coeff * term.head.to_ADD(Algebra, term.data, term)

    def data_to_str_and_precedence(self, cls, (term, coeff)):
        neg_p = heads_precedence.NEG
        mul_p = heads_precedence.MUL
        if term==1:
            t, t_p = NUMBER.data_to_str_and_precedence(cls, coeff)
        elif coeff==1:
            t, t_p = term.head.data_to_str_and_precedence(cls, term.data)
        elif coeff==-1:
            t, t_p = term.head.data_to_str_and_precedence(cls, term.data)
            t, t_p = ('-('+t+')' if t_p < mul_p else '-' + t), neg_p
        elif coeff==0:
            t, t_p = '0', heads_precedence.NUMBER
        else:
            t, t_p = term.head.data_to_str_and_precedence(cls, term.data)
            c, c_p = NUMBER.data_to_str_and_precedence(cls, coeff)
            if t=='1':
                return c, c_p
            cs = '('+c+')' if c_p < mul_p else c
            ts = '('+t+')' if t_p < mul_p else t
            t = cs + (ts[1:] if ts.startswith('1/') else '*' + ts)
            t_p = mul_p
        return t, t_p

    def term_coeff(self, cls, expr):
        return expr.data

    def neg(self, cls, expr):
        term, coeff = expr.data
        return term_coeff_new(cls, (term, -coeff))

    def add(self, cls, lhs, rhs):
        term, coeff = lhs.data
        head, data = rhs.pair
        if head is ADD:
            return ADD.new(cls, [lhs] + data)
        if head is TERM_COEFF_DICT:
            data = data.copy()
            dict_add_item(cls, data, term, coeff)
            return term_coeff_dict_new(cls, data)
        if head is NUMBER:
            if data==0:
                return lhs
            return cls(TERM_COEFF_DICT,{term:coeff, cls(NUMBER,1):data})
        if head is SYMBOL:
            if term==rhs:
                return term_coeff_new(cls, (term, coeff + 1))
            return cls(TERM_COEFF_DICT,{term:coeff, rhs:1})
        if head is TERM_COEFF:
            rterm, rcoeff = data
            if rterm==term:
                return term_coeff_new(cls, (term, coeff + rcoeff))
            return cls(TERM_COEFF_DICT,{term:coeff, rterm:rcoeff})
        if head is BASE_EXP_DICT:
            rcoeff = base_exp_dict_get_coefficient(cls, data)
            if rcoeff is not None:
                d = data.copy()
                del d[rcoeff]
                rterm = base_exp_dict_new(cls, d)
                if rterm==term:
                    return term_coeff_new(cls, (term, coeff + rcoeff))
                return cls(TERM_COEFF_DICT,{term:coeff, rterm:rcoeff})
            else:
                if term==rhs:
                    return term_coeff_new(cls, (term, coeff + 1))
                return cls(TERM_COEFF_DICT,{term:coeff, rhs:1})
        if head is POW or head is APPLY or head is DIFF or head is FDIFF or head is SUBSCRIPT:
            if term==rhs:
                return term_coeff_new(cls, (term, coeff + 1))
            return cls(TERM_COEFF_DICT,{term:coeff, rhs:1})
        if head is MUL:
            return cls(ADD, [lhs, rhs])
        raise NotImplementedError(`self, rhs.head`)

    inplace_add = add

    def add_number(self, cls, lhs, rhs):
        if rhs==0:
            return lhs
        term, coeff = lhs.data
        return cls(TERM_COEFF_DICT, {term: coeff, cls(NUMBER, 1): rhs})

    def sub_number(self, cls, lhs, rhs):
        if rhs==0:
            return lhs
        term, coeff = lhs.data
        return cls(TERM_COEFF_DICT, {term: coeff, cls(NUMBER, 1): -rhs})
    
    def sub(self, cls, lhs, rhs):
        return lhs + (-rhs)

    def non_commutative_mul(self, cls, lhs, rhs):
        term, coeff = lhs.data
        head, data = rhs.pair
        if head is NUMBER:
            return term_coeff_new(cls, (term, coeff * data))
        return (term * rhs) * coeff

    commutative_mul = non_commutative_mul

    def commutative_mul_number(self, cls, lhs, rhs):
        if rhs==1:
            return lhs
        if rhs==0:
            return cls(NUMBER, 0)
        term, coeff = lhs.data
        new_coeff = coeff * rhs
        if new_coeff==1:
            return term
        return cls(TERM_COEFF, (term, new_coeff))

    non_commutative_mul_number = commutative_mul_number

    inplace_commutative_mul = commutative_mul

    def commutative_div_number(self, cls, lhs, rhs):
        term, coeff = lhs.data
        r = number_div(cls, coeff, rhs)
        if rhs==0:
            return r * term
        return term_coeff_new(cls, (term, r))

    def commutative_rdiv_number(self, cls, lhs, rhs):
        term, coeff = lhs.data
        r = number_div(cls, rhs, coeff)
        if rhs==0:
            return r * (1/term)
        return term_coeff_new(cls, (1/term, r))

    def commutative_div(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.commutative_div_number(cls, lhs, rdata)
        term, coeff = lhs.data
        if rhead is TERM_COEFF:
            rterm, rcoeff = rdata
            if term==rterm:
                return cls(NUMBER, number_div(cls, coeff, rcoeff))
            return term.head.commutative_div(cls, term, rterm) * number_div(cls, coeff, rcoeff)
        if rhead is SYMBOL or rhead is APPLY:
            if term == rhs:
                return cls(NUMBER, coeff)
            b, e = term.base_exp()
            return cls(TERM_COEFF, (cls(BASE_EXP_DICT, {b:e, rhs:-1}), coeff))
        if rhead is TERM_COEFF_DICT:
            b, e = term.base_exp()
            return cls(TERM_COEFF, (cls(BASE_EXP_DICT, {b:e, rhs:-1}), coeff))
        if term==rhs:
            return cls(NUMBER, coeff)
        return (term / rhs) * coeff
        return term.head.commutative_div(cls, term, rhs) * coeff

    def pow(self, cls, base, exp):
        term, coeff = base.data
        if isinstance(exp, Expr):
            head, data = exp.pair
            if head is NUMBER:
                exp = data
        if isinstance(exp, inttypes):
            if exp<0:
                return term ** exp / coeff ** (-exp)
            return term ** exp * coeff ** exp
        return pow_new(cls, (base, exp))

    def pow_number(self, cls, base, exp):
        term, coeff = base.data
        if isinstance(exp, inttypes):
            if exp<0:
                return term ** exp / coeff ** (-exp)
            return term ** exp * coeff ** exp
        return cls(POW, (base, exp))

    def walk(self, func, cls, data, target):
        term, coeff = data
        h, d = term.pair
        term = h.walk(func, cls, d, term)
        if isinstance(coeff, Expr):
            h, d = coeff.pair
            coeff = h.walk(func, type(coeff), d, coeff)
        s = term * coeff
        h, d = s.pair
        return func(cls, h, d, s)

    def scan(self, proc, cls, data, target):
        term, coeff = data
        term.head.scan(proc, cls, term.data, target)    
        if isinstance(coeff, Expr):
            coeff.head.scan(proc, type(coeff), coeff.data, target)
        proc(cls, self, data, target)

    def expand(self, cls, expr):
        term, coeff  = expr.data
        return term.head.expand(cls, term) * coeff

    def diff(self, cls, data, expr, symbol, order, cache={}):
        term, coeff = data
        return term.head.diff(cls, term.data, term, symbol, order, cache=cache) * coeff


    def diff_apply(self, cls, data, diff, expr):
        term, coeff = data
        return term.head.diff_apply(cls, term.data, term, expr) * coeff

    def apply(self, cls, data, func, args):
        term, coeff = data
        return term.head.apply(cls, term.data, term, args) * coeff

    def integrate_indefinite(self, cls, data, expr, x):
        term, coeff = data
        return term.head.integrate_indefinite(cls, term.data, term, x) * coeff

    def integrate_definite(self, cls, data, expr, x, a, b):
        term, coeff = data
        return term.head.integrate_definite(cls, term.data, term, x, a, b) * coeff

    def to_TERM_COEFF_DICT(self, Algebra, data, expr):
        term, coeff = data
        return term.head.to_TERM_COEFF_DICT(Algebra, term.data, term) * coeff

    def algebra_pos(self, Algebra, expr):
        return expr

    def algebra_neg(self, Algebra, expr):
        if Algebra.algebra_options.get('evaluate_addition'):
            term, coeff = expr.data
            return term_coeff_new(Algebra, (term, -coeff))
        return Algebra(NEG, expr)

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        return self.algebra_add(Algebra, lhs, Algebra(NUMBER, rhs), inplace)

    def algebra_add(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if Algebra.algebra_options.get('is_additive_group_commutative'):
            lterm, lcoeff = lhs.data
            if Algebra.algebra_options.get('evaluate_addition'):
                if rhead is ADD or rhead is MUL or rhead is NEG:
                    rhs = rhead.to_TERM_COEFF_DICT(Algebra, rdata, rhs)
                    rhead, rdata = rhs.pair
                if rhead is NUMBER:
                    if rdata==0:
                        return lhs
                    return Algebra(TERM_COEFF_DICT, {lterm:lcoeff, Algebra(NUMBER, 1): rdata})
                if rhead is SYMBOL:
                    if lterm==rhs:
                        return term_coeff_new(Algebra, (lterm, lcoeff + 1))
                    return Algebra(TERM_COEFF_DICT, {lterm:lcoeff, rhs: 1})
                if rhead is TERM_COEFF:
                    rterm, rcoeff = rdata
                    if lterm == rterm:
                        return term_coeff_new(Algebra, (lterm, lcoeff+rcoeff))
                    return Algebra(TERM_COEFF_DICT, {lterm:lcoeff, rterm:rcoeff})
                if rhead is TERM_COEFF_DICT:
                    d = rdata.copy()
                    term_coeff_dict_add_item(Algebra, d, lterm, lcoeff)
                    return term_coeff_dict_new(Algebra, d)
            else:
                if rhead is TERM_COEFF_DICT or rhead is EXP_COEFF_DICT:
                    rhead, rdata = rhs.to(ADD).pair
                if rhead is ADD:
                    data = [lhs] + rdata
                else:
                    data = [lhs, rhs]
                return add_new(Algebra, data)
            return super(type(self), self).algebra_add(Algebra, lhs, rhs, inplace)
        else:
            if rhead is TERM_COEFF_DICT or rhead is EXP_COEFF_DICT:
                rhead, rdata = rhs.to(ADD).pair
            if rhead is ADD:
                data = [lhs] + rdata
            else:
                data = [lhs, rhs]
            if Algebra.algebra_options.get('evaluate_addition'):
                ADD.combine_add_list(Algebra, data)
            return add_new(Algebra, data)

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_addition'):
            if not rhs:
                return Algebra(NUMBER, 0)
            term, coeff = lhs.data
            return term_coeff_new(Algebra, (term, coeff * rhs))
        return mul_new(Algebra, [lhs, Algebra(NUMBER, rhs)])

    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_addition'):
            term, coeff = lhs.data
            rhead, rdata = rhs.pair
            if rhead is NUMBER:
                return term_coeff_new(Algebra, (term, coeff*rdata))
            return super(type(self), self).algebra_mul(Algebra, lhs, rhs, inplace)
        return mul_new(Algebra, [lhs, rhs])

TERM_COEFF = TermCoeff()
