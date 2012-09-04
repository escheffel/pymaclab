
__all__ = ['POW']

from .base import BinaryHead, heads_precedence, Head, Expr, Pair, ArithmeticHead

from ..core import init_module
init_module.import_heads()
init_module.import_lowlevel_operations()
init_module.import_numbers()

@init_module
def _init(module):
    from ..arithmetic.numbers import try_power
    module.try_power = try_power

class PowHead(ArithmeticHead):
    """ PowHead represents exponentiation operation, data is a 2-tuple
    of base and exponent expressions. Both can be number instances or
    algebra instances.
    """
    op_mth = '__pow__'
    op_rmth = '__rpow__'

    def is_data_ok(self, cls, data):
        if type(data) is tuple and len(data)==2:
            base, exp = data
            if isinstance(base, cls):
                if isinstance(exp, numbertypes):
                    return
                if isinstance(exp, cls):
                    if exp.head is NUMBER:
                        if isinstance(exp.data, numbertypes):
                            return 'data[1] must be lowlevel number or non-numeric but got %s' % (type(exp.data))
                    else:
                        return
                else:
                    return 'data[1] must be %s instance but got %s' % ((cls, numbertypes), type(exp))
            else:
                return 'data[0] must be %s instance but got %s' % (cls, type(exp))
        else:
            return 'data must be 2-tuple'
        return

    def __repr__(self): return 'POW'

    def new(self, cls, (base, exp), evaluate=True):
        if exp==1:
            return base
        if exp==0 or base==1:
            return cls(NUMBER, 1)
        if not evaluate:
            return cls(self, (base, exp))
        if type(exp) is cls:
            h, d = exp.pair
            if h is NUMBER:
                exp = d
        if base.head is NUMBER and isinstance(exp, numbertypes):
            b = base.data
            if isinstance(b, numbertypes):
                r, base_exp_list = try_power(b, exp)
                if not base_exp_list:
                    return cls(NUMBER, r)
                if len(base_exp_list)==1:
                    b, e = base_exp_list[0]
                    rest = cls(POW, (cls(NUMBER, b), e))
                else:
                    d = {}
                    for b, e in base_exp_list:
                        d[cls(NUMBER, b)] = e
                    rest = cls(BASE_EXP_DICT, d)
                if r==1:
                    return rest
                return  cls(TERM_COEFF, (rest, r))
        return cls(self, (base, exp))

    def reevaluate(self, cls, (base, exp)):
        return base ** exp

    def data_to_str_and_precedence(self, cls, (base, exp)):
        pow_p = heads_precedence.POW
        div_p = heads_precedence.DIV
        if isinstance(base, Expr):
            b, b_p = base.head.data_to_str_and_precedence(cls, base.data)
        elif isinstance(base, numbertypes):
            b, b_p = NUMBER.data_to_str_and_precedence(cls, base)
        else:
            b, b_p = SYMBOL.data_to_str_and_precedence(cls, base)

        if isinstance(exp, Expr):
            h, d = exp.pair
            if h is NUMBER and isinstance(d, numbertypes):
                exp = d
        if isinstance(exp, numbertypes):
            if exp==0:
                return '1', heads_precedence.NUMBER
            if exp==1:
                return b, b_p
            if exp < 0:
                if exp==-1:
                    s1 = '('+b+')' if b_p <= pow_p else b
                    return '1/' + s1, div_p
                e, e_p = NUMBER.data_to_str_and_precedence(cls, -exp)
                s1 = '('+b+')' if b_p < pow_p else b
                s2 = '('+e+')' if e_p < pow_p else e
                return '1/' + s1 + '**' + s2, div_p
            e, e_p = NUMBER.data_to_str_and_precedence(cls, exp)
        else:
            if isinstance(exp, Expr):
                e, e_p = exp.head.data_to_str_and_precedence(cls, exp.data)
            else:
                e, e_p = str(exp), 0.0
        s1 = '('+b+')' if b_p <= pow_p else b
        s2 = '('+e+')' if e_p < pow_p else e
        return s1 + '**' + s2, pow_p

    def to_ADD(self, Algebra, (base, exp), expr):
        return Algebra(ADD, [expr])

    def to_TERM_COEFF_DICT(self, cls, data, expr):
        return expr

    def to_EXP_COEFF_DICT(self, cls, (base, exp), expr, variables=None):
        if isinstance(exp, Expr):
            if exp.head is NUMBER:
                exp = exp.data
            elif exp.head is TERM_COEFF:
                t. c = exp.data
                return self.to_EXP_COEFF_DICT(cls, (base**t, c), expr, variables)
        if isinstance(exp, inttypes):
            return base.head.to_EXP_COEFF_DICT(cls, base.data, base, variables) ** exp
        if isinstance(exp, rationaltypes):
            numer, denom = exp
            if numer!=1:
                return self.to_EXP_COEFF_DICT(cls, (base**(exp/numer), numer), expr, variables)
        raise NotImplementedError(`base, exp`)

    def non_commutative_mul(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return term_coeff_new(cls, (lhs, rdata))
        if rhead is SYMBOL or rhead is POW:
            return MUL.combine(cls, [lhs, rhs])
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return (lhs * term) * coeff
        if rhead is DIFF:
            base, exp = lhs.data
            if base==rhs:
                return pow_new(cls, (base, exp + 1))
            return MUL.combine(cls, [lhs, rhs])
        if rhead is MUL:
            return MUL.combine(cls, [lhs] + rdata)
        raise NotImplementedError(`self, cls, lhs.pair, rhs.pair`)

    def commutative_mul_number(self, cls, lhs, rhs):
        return term_coeff_new(cls, (lhs, rhs))

    non_commutative_mul_number = commutative_mul_number
    non_commutative_rmul_number = commutative_mul_number

    def commutative_mul(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return term_coeff_new(cls, (lhs, rdata))
        if rhead is SYMBOL or rhead is ADD or rhead is TERM_COEFF_DICT or rhead is APPLY or rhead is DIFF or rhead is FDIFF:
            lbase, lexp = lhs.data
            if lbase == rhs:
                return pow_new(cls, (lbase, lexp + 1))
            return cls(BASE_EXP_DICT, {rhs:1, lbase:lexp})
        if rhead is POW:
            lbase, lexp = lhs.data
            rbase, rexp = rdata
            if lbase==rbase:
                return POW.new(cls, (lbase, lexp + rexp))
            return cls(BASE_EXP_DICT, {lbase:lexp, rbase:rexp})
        if rhead is BASE_EXP_DICT:
            base, exp = lhs.data
            data = rhs.data.copy()
            dict_add_item(cls, data, base, exp)
            return base_exp_dict_new(cls, data)
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return (lhs * term) * coeff
        raise NotImplementedError(`self, cls, lhs.pair, rhs.pair`)

    inplace_commutative_mul = commutative_mul

    def commutative_div_number(self, cls, lhs, rhs):
        r = number_div(cls, 1, rhs)
        if rhs==0:
            return r * lhs
        return term_coeff_new(cls, (lhs, r))

    def commutative_rdiv_number(self, cls, lhs, rhs):
        base, exp = lhs.data
        return term_coeff_new(cls, (pow_new(cls, (base, -exp)), rhs))

    def commutative_div(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.commutative_div_number(cls, lhs, rdata)
        base, exp = lhs.data
        if rhead is POW:
            rbase, rexp = rdata
            if base==rbase:
                return pow_new(cls, (base, exp-rexp))
            return base_exp_dict_new(cls, {base:exp, rbase: -rexp})
        if rhead is BASE_EXP_DICT:
            data = {base:exp}
            for b, e in rdata.iteritems():
                base_exp_dict_add_item(cls, data, b, -e)
            return base_exp_dict_new(cls, data)
        if rhead is TERM_COEFF:
            term, coeff = rhs.term_coeff()
            return (lhs / term) / coeff
        if base==rhs:
            return pow_new(cls, (base, exp-1))
        return base_exp_dict_new(cls, {base:exp, rhs:-1})
    
    def base_exp(self, cls, expr):
        base, exp = expr.data
        return base, exp

    def pow(self, cls, base, exp):
        if exp==0:
            return cls(NUMBER, 1)
        if exp==1:
            return base
        if isinstance(exp, Expr) and exp.head is NUMBER:
            exp = exp.data
        if isinstance(exp, inttypes):
            b, e = base.data
            base, exp = b, e*exp
        
        return POW.new(cls, (base, exp))

    pow_number = pow

    def walk(self, func, cls, data, target):
        base, exp = data
        base1 = base.head.walk(func, cls, base.data, base)
        if isinstance(exp, Expr):
            exp1 = exp.head.walk(func, cls, exp.data, exp)
        else:
            exp1 = NUMBER.walk(func, cls, exp, exp)
        if base1 is base and exp1 is exp:
            return func(cls, self, data, target)
        else:
            r = base1 ** exp1
            return func(cls, r.head, r.data, r)

    def scan(self, proc, cls, data, target):
        base, exp = data
        base.head.scan(proc, cls, base.data, target)
        if isinstance(exp, Expr):
            exp.head.scan(proc, cls, exp.data, target)
        else:
            NUMBER.scan(proc, cls, exp, target)
        proc(cls, self, data, target)

    def expand(self, cls, expr):
        base, exp = expr.data
        if isinstance(exp, Expr):
            exp = exp.expand()
            h, d = exp.pair
            if h is NUMBER and isinstance(d, int):
                exp = d
        if isinstance(base, Expr):
            base = base.expand()
            if isinstance(exp, int):
                return base.head.expand_intpow(cls, base, exp)
        return cls(POW, (base, exp))

    def diff(self, cls, data, expr, symbol, order, cache={}):
        # XXXX needs implementatiin
        key = (expr, symbol, order)
        result = cache.get(key)
        if result is not None:
            return result
        base, exp = data
        texp = type(exp)
        if symbol not in base.symbols_data:
            # constant ** exp
            if texp is cls:
                if exp.head is SYMBOL:
                    if exp.data==symbol:
                        result = expr * cls.Log(base)**order
                        cache[key] = result
                        return result
                    else:
                        return cls(NUMBER, 0)
                if symbol not in exp.symbols_data:
                    return cls(NUMBER, 0)
                key1 = (expr, symbol, 1)
                result = cache.get(key1)
                if result is None:
                    de = exp.head.diff(cls, exp.data, exp, symbol, 1, cache=cache)
                    if symbol not in de.symbols_data:
                        result = expr * de**order * cls.Log(base)**order
                        cache[key] = result
                        return result
                    result = expr * cls.Log(base) * de
                    cache[key1] = result
                if order>1:
                    result = result.head.diff(cls, result.data, result, symbol, order-1, cache=cache)
                    cache[key] = result
                return result
            else:
                return cls(NUMBER, 0)
        elif not (texp is cls and symbol in exp.symbols_data):
            if exp is cls and exp.head is NUMBER:
                exp = exp.data
            # variable ** constant            
            # f(x)**n  -> n*f**(n-1)*f' ->
            db = base.head.diff(cls, base.data, base, symbol, 1, cache=cache)
            if db.head is NUMBER:
                if texp is int and order>exp and exp>0:
                    return cls(NUMBER, 0)            
                p = db.data ** order
                e = exp
                for i in xrange(order):
                    p *= e
                    e -= 1
                result = p * base ** e
                cache[key] = result
                return result

        key1 = (expr, symbol, 1)
        result = cache.get(key1)
        if result is None:
            base, exp = data
            db = base.head.diff(cls, base.data, base, symbol, 1, cache=cache)            
            if isinstance(exp, Expr):
                de = exp.head.diff(cls, exp.data, exp, symbol, 1, cache=cache)
                if de==0:
                    result = base ** (exp-1) * db * exp
                else:
                    result = expr * (db / base * exp + cls.Log(base) * de)
            else:
                result = base ** (exp-1) * db * exp
            cache[key1] = result
        if order>1:
            # todo: if, say, order>400, return unevaluated result to avoid recursion
            result = result.head.diff(cls, result.data, result, symbol, order-1, cache=cache)
            cache[key] = result
        return result

    def diff_apply(self, cls, data, diff, expr):
        base, exp = data
        bhead, bdata = base.pair
        if bhead is DIFF and isinstance(exp, inttypes) and exp>=0:
            return expr.head.diff(type(expr), expr.data, expr, bdata.data, exp)
        if bhead is FDIFF:# and isinstance(exp, inttypes) and exp>=0:
            return expr.head.fdiff(type(expr), expr.data, expr, bdata.data, exp)
        return NotImplemented

    def apply(self, cls, data, func, args):
        base, exp = data
        if isinstance(exp, Expr):
            return NotImplemented
        return base.head.apply(cls, base.data, base, args) ** exp

    def integrate_indefinite(self, cls, data, expr, x):
        base, exp = data
        t = type(exp)
        if t is cls:
            ehead, edata = exp.pair
            exp_symbols_data = exp.symbols_data
        else:
            ehead, edata = NUMBER, exp
            exp_symbols_data = set([])
        bhead, bdata = base.pair
        if bhead is SYMBOL and bdata==x and x not in exp_symbols_data:
            e1 = edata + 1
            return pow_new(cls, (base, e1))/e1
        if bhead is NUMBER or x not in base.symbols_data:
            if ehead is SYMBOL and edata == x:
                return expr / cls.Log(base)
            elif ehead is TERM_COEFF:
                t, c = edata
                if t.head is SYMBOL and t.data==x:
                    return expr / (cls.Log(base) * c)
            if x not in exp_symbols_data:
                return expr * cls(SYMBOL, x)
        raise NotImplementedError("don't know how to integrate %s over %s" % (expr, x))

    def integrate_definite(self, cls, data, expr, x, a, b):
        base, exp = data
        t = type(exp)
        if t is cls:
            ehead, edata = exp.pair
            exp_symbols_data = exp.symbols_data
        else:
            ehead, edata = NUMBER, exp
            exp_symbols_data = set([])
        bhead, bdata = base.pair
        if bhead is SYMBOL and bdata==x and x not in exp_symbols_data:
            e1 = edata + 1
            return (b**e1 - a**e1)/e1
        if bhead is NUMBER or x not in base.symbols_data:
            if ehead is SYMBOL and edata == x:
                return (base**b - base**a) / cls.Log(base)
            elif ehead is TERM_COEFF:
                t, c = edata
                if t.head is SYMBOL and t.data==x:
                    return (base**(b*c) - base**(a*c)) / (cls.Log(base) * c)
            if x not in exp_symbols_data:
                return expr * (b - a)
        raise NotImplementedError("don't know how to integrate %s over %s in [%s, %s]" % (expr, x, a, b))

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        return self.algebra_mul(Algebra, lhs, Algebra(NUMBER, rhs), inplace)
    
    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('is_multiplicative_group'):
            if Algebra.algebra_options.get('evaluate_multiplication'):
                rhead, rdata = rhs.pair
                if Algebra.algebra_options.get('is_multiplicative_group_commutative'):
                    pass
                else:
                    if rhead is MUL:
                        data = [lhs] + rdata
                    else:
                        data = [lhs, rhs]
                    coeff = MUL.combine_mul_list(Algebra, data)
                    assert coeff is 1,`coeff`
                    return mul_new(Algebra, data)
                        
        return super(type(self), self).algebra_mul(Algebra, lhs, rhs, inplace)

    def algebra_pow_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_multiplication'):
            if rhs==1:
                return lhs
            if not rhs:
                return Algebra(NUMBER, 1)
            if rhs==-1:
                base, exp = lhs.data
                return pow_new(Algebra, (base, -exp))
        return super(type(self), self).algebra_pow_number(Algebra, lhs, rhs, inplace)
POW = PowHead()
