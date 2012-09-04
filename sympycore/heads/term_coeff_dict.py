
__all__ = ['TERM_COEFF_DICT']

from .base import heads, heads_precedence, ArithmeticHead, Pair

from ..core import init_module, Expr
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

@init_module
def _init(module):
    from ..arithmetic.number_theory import multinomial_coefficients
    module.multinomial_coefficients = multinomial_coefficients

class TermCoeffDictHead(ArithmeticHead):

    def is_data_ok(self, cls, data):
        if type(data) is dict:
            n = len(data)
            #if n<=1:
            #    return 'data dictonary should have more than 1 item'
            for item in data.iteritems():
                msg = TERM_COEFF.is_data_ok(cls, item, allow_number_term=True)
                if msg:
                    return 'TERM_COEFF data=%s: %s' % (item, msg) #pragma: no cover
        else:
            return 'data must be dict instance but got %s' % (type(data)) #pragma: no cover
        return
    
    def __repr__(self): return 'TERM_COEFF_DICT'

    def data_to_str_and_precedence(self, cls, term_coeff_dict):
        r = [cls(TERM_COEFF, tc) for tc in term_coeff_dict.items()]
        return ADD.data_to_str_and_precedence(cls, r)

    def new(self, cls, data, evaluate=True):
        return term_coeff_dict_new(cls, data)

    def reevaluate(self, cls, data):
        r = cls(NUMBER, 0)
        for term, coeff in data.iteritems():
            r += term * coeff
        return r

    def to_ADD(self, Algebra, data, expr):
        return add_new(Algebra, [term_coeff_new(Algebra, term_coeff) for term_coeff in data.iteritems()])

    def to_EXP_COEFF_DICT(self, cls, data, expr, variables=None):
        if variables is None:
            r = cls(EXP_COEFF_DICT, Pair((), {}))
        else:
            r = cls(EXP_COEFF_DICT, Pair(variables, {}))
        for term, coeff in data.iteritems():
            r += term * coeff
        return r

    def term_coeff(self, cls, expr):
        term_coeff_dict = expr.data
        if len(term_coeff_dict)==1:
            return dict_get_item(term_coeff_dict)
        return expr, 1

    def inplace_add(self, cls, lhs, rhs):
        if lhs.is_writable:
            data = lhs.data
        else:
            data = lhs.data.copy()

        self.add(cls, data, rhs, inplace=True)

        return term_coeff_dict_new(cls, data)

    def add(self, cls, lhs, rhs, inplace=False):
        if inplace:
            data = lhs
        else:
            data = lhs.data.copy()
        
        h2, d2 = rhs.pair
        if h2 is NUMBER:
            if d2 != 0:
                dict_add_item(cls, data, cls(NUMBER, 1), d2)
        elif h2 is SYMBOL:
            dict_add_item(cls, data, rhs, 1)
        elif h2 is TERM_COEFF:
            term, coeff = d2
            dict_add_item(cls, data, term, coeff)
        elif h2 is TERM_COEFF_DICT:
            dict_add_dict(cls, data, d2)
        elif h2 is ADD:
            for op in d2:
                term, coeff = op.term_coeff()
                dict_add_item(cls, data, term, coeff)
        elif h2 is SUB or h2 is NEG or h2 is POS:
            raise NotImplementedError(`self, rhs.pair`)
        elif h2 is BASE_EXP_DICT:
            c = base_exp_dict_get_coefficient(cls, d2)
            if c is not None:
                d = d2.copy()
                del d[c]
                t = BASE_EXP_DICT.new(cls, d)
                if t.head is BASE_EXP_DICT:
                    dict_add_item(cls, data, t, c)
                else:
                    self.add(cls, data, t * c, inplace=True)
            else:
                dict_add_item(cls, data, rhs, 1)
        else:
            dict_add_item(cls, data, rhs, 1)

        if inplace:
            return lhs

        return term_coeff_dict_new(cls, data)

    def add_number(self, cls, lhs, rhs):
        if rhs==0:
            return lhs
        data = lhs.data.copy()
        term_coeff_dict_add_item(cls, data, cls(NUMBER, 1), rhs)
        return term_coeff_dict_new(cls, data)

    def sub_number(self, cls, lhs, rhs):
        if rhs==0:
            return lhs
        data = lhs.data.copy()
        term_coeff_dict_add_item(cls, data, cls(NUMBER, 1), -rhs)
        return term_coeff_dict_new(cls, data)

    def sub(self, cls, lhs, rhs):
        d = lhs.data.copy()
        self.add(cls, d, -rhs, inplace=True)
        return term_coeff_dict_new(cls, d)

    def commutative_mul(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            if rdata==0:
                return rhs
            data = lhs.data.copy()
            dict_mul_value(cls, data, rdata)
            return term_coeff_dict_new(cls, data)
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return (lhs * term) * coeff
        if rhead is POW:
            base, exp = rdata
            if lhs==base:
                return POW.new(cls, (lhs, exp + 1))
            return cls(BASE_EXP_DICT, {lhs:1, base:exp})
        if rhead is SYMBOL or rhead is APPLY:
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        if rhead is TERM_COEFF_DICT:
            if rdata==lhs.data:
                return cls(POW, (lhs, 2))
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        if rhead is BASE_EXP_DICT:
            data = rdata.copy()
            dict_add_item(cls, data, lhs, 1)
            return BASE_EXP_DICT.new(cls, data)
        if rhead is ADD:
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        return ArithmeticHead.commutative_mul(self, cls, lhs, rhs)

    def commutative_mul_number(self, cls, lhs, rhs):
        if rhs==0:
            return cls(NUMBER, 0)
        data = lhs.data.copy()
        dict_mul_value(cls, data, rhs)
        return cls(self, data)

    non_commutative_mul_number = commutative_mul_number

    def commutative_rdiv_number(self, cls, lhs, rhs):
        return term_coeff_new(cls, (cls(POW, (lhs, -1)), rhs))

    def commutative_div(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.commutative_div_number(cls, lhs, rdata)
        if rhead is TERM_COEFF_DICT:
            if lhs.data == rdata:
                return cls(NUMBER, 1)
            return cls(BASE_EXP_DICT, {lhs:1, rhs:-1})
        if rhead is SYMBOL or rhead is APPLY:
            return cls(BASE_EXP_DICT, {lhs:1, rhs:-1})
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return number_div(cls, 1, coeff) * (lhs / term)
        if rhead is POW:
            base, exp = rdata
            if lhs==base:
                return pow_new(cls, (lhs, 1-exp))
            return cls(BASE_EXP_DICT, {lhs:1, base:-exp})            
        if rhead is BASE_EXP_DICT:
            data = {lhs:1}
            for base, exp in rdata.iteritems():
                base_exp_dict_add_item(cls, data, base, -exp)
            return base_exp_dict_new(cls, data)
        return ArithmeticHead.commutative_div(self, cls, lhs, rhs)

        
    def pow(self, cls, base, exp):
        if exp==0: return cls(NUMBER, 1)
        if exp==1: return base
        d = base.data
        if len(d)==1:
            t,c = dict_get_item(d)
            t,c = t**exp, c**exp
            if t==1: return cls(NUMBER, c)
            if c==1: return t
            return cls(TERM_COEFF, (t, c))
        return POW.new(cls, (base, exp))

    pow_number = pow

    def neg(self, cls, expr):
        d = expr.data.copy()
        for key in d:
            d[key] = -d[key]
        return cls(TERM_COEFF_DICT, d)

    def expand(self, cls, expr):
        d = {}
        for t, c in expr.data.items():
            self.add(cls, d, t.expand() * c, inplace=True)
        return term_coeff_dict_new(cls, d)

    def expand_intpow(self, cls, expr, intexp):
        if intexp<0:
            return cls(POW, (expr, intexp))
        if intexp==0:
            return cls(NUMBER, 1)
        if intexp==1:
            return expr
        term_coeff_list = [(term.base_exp(), coeff) for term, coeff in expr.data.items()]
        mdata = multinomial_coefficients(len(term_coeff_list), intexp)
        d = {}
        for e,c in mdata.iteritems():
            new_coeff = c
            df = {}
            for e_i, ((base, exp), coeff) in zip(e, term_coeff_list):
                if e_i:
                    if e_i==1:
                        base_exp_dict_add_item(cls, df, base, exp)
                        if coeff is not 1:
                            new_coeff *= coeff
                    else:
                        base_exp_dict_add_item(cls, df, base, exp*e_i)
                        if coeff is not 1:
                            new_coeff *= coeff ** e_i
            new_term = base_exp_dict_new(cls, df)
            term_coeff_dict_add_item(cls, d, new_term, new_coeff)
        return term_coeff_dict_new(cls, d)

    def walk(self, func, cls, data, target):
        d = {}
        flag = False
        add = self.add
        for t, c in data.iteritems():
            t1 = t.head.walk(func, cls, t.data, t)
            if isinstance(c, Expr):
                c1 = c.head.walk(func, cls, c.data, c)
            else:
                c1 = NUMBER.walk(func, cls, c, c)
            if t1 is not t or c1 is not c:
                flag = True
            add(cls, d, t1 * c1, inplace=True)
        if flag:
            r = term_coeff_dict_new(cls, d)
            return func(cls, r.head, r.data, r)
        return func(cls, self, data, target)

    def scan(self, proc, cls, data, target):
        for t, c in data.iteritems():
            t.head.scan(proc, cls, t.data, target)
            if isinstance(c, Expr):
                c.head.scan(proc, cls, c.data, target)
            else:
                NUMBER.scan(proc, cls, c, target)
        proc(cls, self, data, target)

    def diff(self, cls, data, expr, symbol, order, cache={}):
        key = (expr, symbol, order)
        result = cache.get(key)
        if result is not None:
            return result
        if result is None:
            d = {}
            result = cls(NUMBER, 0)
            for term, coeff in data.iteritems():
                result += term.head.diff(cls, term.data, term, symbol, order, cache=cache) * coeff
            key1 = (expr, symbol, 1)
        cache[key] = result
        return result

    def apply(self, cls, data, func, args):
        result = cls(NUMBER, 0)
        for term, coeff in data.iteritems():
            result += term.head.apply(cls, term.data, term, args) * coeff
        return result

    def integrate_indefinite(self, cls, data, expr, x):
        result = cls(NUMBER, 0)
        for term, coeff in data.iteritems():
            result += term.head.integrate_indefinite(cls, term.data, term, x) * coeff
        return result

    def integrate_definite(self, cls, data, expr, x, a, b):
        result = cls(NUMBER, 0)
        for term, coeff in data.iteritems():
            result += term.head.integrate_definite(cls, term.data, term, x, a, b) * coeff
        return result

    def algebra_pos(self, Algebra, expr):
        return expr

    def algebra_neg(self, Algebra, expr):
        if Algebra.algebra_options.get('evaluate_addition'):
            d = expr.data.copy()
            for key in d:
                d[key] = -d[key]
            return Algebra(TERM_COEFF_DICT, d)
        return Algebra(NEG, expr)

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        if not rhs:
            return lhs
        if inplace:
            term_coeff_dict_add_item(Algebra, lhs.data, Algebra(NUMBER, 1), rhs)
            return term_coeff_dict(Algebra, lhs)
        d = lhs.data.copy()
        term_coeff_dict_add_item(Algebra, d, Algebra(NUMBER, 1), rhs)
        return term_coeff_dict_new(Algebra, d)

    def algebra_add(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if Algebra.algebra_options.get('evaluate_addition'):
            ldata = lhs.data
            if rhead is ADD or rhead is EXP_COEFF_DICT or rhead is MUL or rhead is NEG:
                rhs = rhead.to_TERM_COEFF_DICT(Algebra, rdata, rhs)
                rhead, rdata = rhs.pair
            if rhead is NUMBER:
                if not rdata:
                    return lhs
                rterm, rcoeff = Algebra(NUMBER, 1), rdata
            elif rhead is SYMBOL:
                rterm, rcoeff = rhs, 1
            elif rhead is TERM_COEFF:
                rterm, rcoeff = rdata
            elif rhead is TERM_COEFF_DICT:
                if inplace:
                    term_coeff_dict_add_dict(Algebra, ldata, rdata)
                    return term_coeff_dict(Algebra, lhs)
                d = ldata.copy()
                term_coeff_dict_add_dict(Algebra, d, rdata)
                return term_coeff_dict_new(Algebra, d)
            else:
                return super(type(self), self).algebra_add(Algebra, lhs, rhs, inplace)

            if inplace:
                term_coeff_dict_add_item(Algebra, ldata, rterm, rcoeff)
                return term_coeff_dict(Algebra, lhs)
            d = ldata.copy()
            term_coeff_dict_add_item(Algebra, d, rterm, rcoeff)
            return term_coeff_dict_new(Algebra, d)
        else:
            return TERM_COEFF_DICT.to_ADD(Algebra, lhs.data, lhs) + rhs
        return super(type(self), self).algebra_add(Algebra, lhs, rhs, inplace)

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_addition'):
            if not rhs:
                return Algebra(NUMBER, 0)
            if rhs==1:
                return lhs
            if inplace:
                term_coeff_dict_mul_value(Algebra, lhs.data, rhs)
                return lhs
            d = lhs.data.copy()
            term_coeff_dict_mul_value(Algebra, d, rhs)
            return Algebra(TERM_COEFF_DICT, d)
        return Algebra(MUL, [lhs, Algebra(NUMBER, rhs)])

    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.algebra_mul_number(Algebra, lhs, rdata, inplace)
        return super(type(self), self).algebra_mul(Algebra, lhs, rhs, inplace)

    def algebra_div_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_addition'):
            if rhs==1:
                return lhs
            d1 = gcd(*lhs.data.values())
            d2 = gcd(d1, rhs)
            d3 = rhs / d2
            d = {}
            rd = 0
            for t,c in lhs.data.items():
                c /= d2
                q, c = divmod(c, d3)
                if c:
                    d[t] = c
                rd += q
            s = term_coeff_dict_new(Algebra, d)
            if rhs==d2:
                assert rd==0,`lsh, rhs, s,rd`
                return s
            return Algebra(DIV, [s, Algebra(NUMBER, d3)]) + rd
        return Algebra(DIV, [lhs, Algebra(NUMBER, rhs)])

    def algebra_div(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.algebra_div_number(Algebra, lhs, rdata, inplace)
        return super(type(self), self).algebra_div(Algebra, lhs, rhs, inplace)

TERM_COEFF_DICT = TermCoeffDictHead()
