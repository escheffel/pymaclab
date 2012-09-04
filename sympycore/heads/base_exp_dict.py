
__all__ = ['BASE_EXP_DICT']

from .base import heads, heads_precedence, ArithmeticHead

from ..core import init_module, Expr
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

class BaseExpDictHead(ArithmeticHead):
    """ BASE_EXP_DICT expression data is a dictionary of base and
    exponent pairs. All base parts must be Expr instances.

    For example, ``Algebra(BASE_EXP_DICT. {x:2, y:a, 3:1, 2:1/2})``
    represents ``3 * 2**(1/2) * x**2 * y**a``.
    """

    def is_data_ok(self, cls, data):
        if type(data) is dict:
            for item in data.iteritems():
                msg = POW.is_data_ok(cls, item)
                if msg:
                    return 'POW data=%s: %s' % (item, msg)
                b, e = item
                if b.head is POW:
                    return 'BASE_EXP_DICT key cannot be POW'
        else:
            return 'data must be dict instance but got %s' % (type(data))
        return


    def __repr__(self): return 'BASE_EXP_DICT'

    def data_to_str_and_precedence(self, cls, base_exp_dict):
        factors = []
        coeff = None
        for base, exp in base_exp_dict.items():
            if exp==1 and base.head is NUMBER:
                coeff = base.data
            else:
                factors.append(cls(POW, (base, exp)))
        if coeff is not None:
            return TERM_COEFF.data_to_str_and_precedence(cls, (cls(MUL, factors), coeff))
        return MUL.data_to_str_and_precedence(cls, factors)

    def reevaluate(self, cls, data):
        r = cls(NUMBER, 1)
        for base, exp in data.iteritems():
            r *= base ** exp
        return r

    def to_ADD(self, Algebra, base_exp_dict, expr):
        return Algebra(ADD, [expr])

    def term_coeff(self, cls, expr):
        data = expr.data
        coeff = base_exp_dict_get_coefficient(cls, data)
        if coeff is not None:
            data = data.copy()
            del data[coeff]
            r = base_exp_dict_new(cls, data)
            t, c = r.head.term_coeff(cls, r)
            return t, c * coeff 
        return expr, 1

    def new(self, cls, base_exp_dict, evaluate=True):
        return base_exp_dict_new(cls, base_exp_dict)

    def neg(self, cls, expr):
        data = expr.data
        coeff = base_exp_dict_get_coefficient(cls, data)
        if coeff is None:
            return cls(TERM_COEFF, (expr, -1))
        data = data.copy()
        del data[coeff]
        return term_coeff_new(cls, (base_exp_dict_new(cls, data), -coeff))

    def inplace_commutative_data_mul(self, cls, data, rhs):
        """
        Multiply base-exp-dictionary with rhs inplace.
        """
        rhead, rdata = rhs.pair
        if rhead is SYMBOL or rhead is ADD or rhead is APPLY or rhead is DIFF or rhead is FDIFF:
            base_exp_dict_add_item(cls, data, rhs, 1)
        elif rhead is NUMBER:
            base_exp_dict_add_item(cls, data, rhs, 1)
        elif rhead is TERM_COEFF:
            term, coeff = rdata
            base_exp_dict_add_item(cls, data, term, 1)
            base_exp_dict_add_item(cls, data, cls(NUMBER, coeff), 1)
        elif rhead is BASE_EXP_DICT:
            base_exp_dict_add_dict(cls, data, rdata)
        elif rhead is POW:
            base, exp = rdata
            base_exp_dict_add_item(cls, data, base, exp)
        elif rhead is TERM_COEFF_DICT:
            base_exp_dict_add_item(cls, data, rhs, 1)
        else:
            raise NotImplementedError(`self, cls, rhs.pair`)
    
    def commutative_mul(self, cls, lhs, rhs):
        data = lhs.data.copy()
        self.inplace_commutative_data_mul(cls, data, rhs)
        return base_exp_dict_new(cls, data)

    def commutative_mul_number(self, cls, lhs, rhs):
        return term_coeff_new(cls, (lhs, rhs))

    def commutative_div_number(self, cls, lhs, rhs):
        r = number_div(cls, 1, rhs)
        if rhs==0:
            return r * lhs
        return term_coeff_new(cls, (lhs, r))

    def commutative_div(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.commutative_div_number(cls, lhs, rdata)
        if rhead is POW:
            data = lhs.data.copy()
            base, exp = rdata
            base_exp_dict_sub_item(cls, data, base, exp)
            return base_exp_dict_new(cls, data)
        if rhead is BASE_EXP_DICT:
            data = lhs.data.copy()
            base_exp_dict_sub_dict(cls, data, rdata)
            return base_exp_dict_new(cls, data)
        if rhead is SYMBOL or rhead is TERM_COEFF_DICT or rhead is APPLY:
            data = lhs.data.copy()
            base_exp_dict_sub_item(cls, data, rhs, 1)
            return base_exp_dict_new(cls, data)
        if rhead is TERM_COEFF:
            term, coeff = rhs.term_coeff()
            return (lhs / term) / coeff
        return ArithmeticHead.commutative_div(self, cls, lhs, rhs)

    def commutative_rdiv_number(self, cls, lhs, rhs):
        data = lhs.data.copy()
        base_exp_dict_mul_value(cls, data, -1)
        return base_exp_dict_new(cls, data) * rhs
    
    def scan(self, proc, cls, data, target):
        for b, e in data.iteritems():
            b.head.scan(proc, cls, b.data, target)
            if isinstance(e, Expr):
                e.head.scan(proc, cls, e.data, target)
            else:
                NUMBER.scan(proc, cls, e, target)
        proc(cls, self, data, target)

    def walk(self, func, cls, data, target):
        d = {}
        flag = False
        for b, e in data.iteritems():
            b1 = b.head.walk(func, cls, b.data, b)
            if isinstance(e, Expr):
                e1 = e.head.walk(func, cls, e.data, e)
            else:
                e1 = NUMBER.walk(func, cls, e, e)
            if b1 is not b or e1 is not e:
                flag = True
            self.inplace_commutative_data_mul(cls, d, b1**e1)
        if flag:
            r = base_exp_dict_new(cls, d)
            return func(cls, r.head, r.data, r)
        return func(cls, self, data, target)

    def pow(self, cls, base, exp):
        if type(exp) is cls:
            h, d = exp.pair
            if h is NUMBER and isinstance(d, numbertypes):
                exp = d
        if isinstance(exp, inttypes):
            if exp:
                data = base.data.copy()
                base_exp_dict_mul_value(cls, data, exp)
                return base_exp_dict_new(cls, data)
            return cls(NUMBER, 1)
        return pow_new(cls, (base, exp))

    pow_number = pow

    def expand(self, cls, expr):
        data = {}
        for b, e in expr.data.items():
            f = pow_new(cls, (b, e)).expand()
            h, d = f.pair
            data1 = {}
            if h is TERM_COEFF_DICT:
                data2 = d
            else:
                t, c = f.term_coeff()
                data2 = {t: c}
            if data:
                term_coeff_dict_mul_dict(cls, data1, data, data2)
                data = data1
            else:
                data = data2
        return term_coeff_dict_new(cls, data)

    def diff(self, cls, data, expr, symbol, order, cache={}):
        key = (expr, symbol, order)
        result = cache.get(key)
        if result is not None:
            return result
        key1 = (expr, symbol, 1)
        result = cache.get(key1)
        if result is None:
            operands = data.items()
            zero = cls(NUMBER, 0)
            result = zero
            for i in range(len(operands)):
                p = pow_new(cls, operands[i])
                d = p.head.diff(cls, p.data, p, symbol, 1, cache=cache)
                if d==zero:
                    continue
                be_dict = data.copy()
                del be_dict[operands[i][0]]
                r = base_exp_dict_new(cls, be_dict)
                result += r * d
            cache[key1] = result
        if order>1:
            result = result.head.diff(cls, result.data, result, symbol, order-1, cache=cache)
            cache[key] = result
        return result

    def apply(self, cls, data, func, args):
        result = cls(NUMBER, 1)
        for base, exp in data.iteritems():
            if isinstance(exp, Expr):
                return NotImplemented
            result *= base.head.apply(cls, base.data, base, args) ** exp
        return result

    def integrate_indefinite(self, cls, data, expr, x):
        d1 = {} # f(x)**g(x)
        d2 = {} # f(x)**const
        d3 = {} # const**g(x)
        d4 = {} # const**const
        for base, exp in data.iteritems():
            if x in base.symbols_data:
                if type(exp) is cls and x in exp.symbols_data:
                    d1[base] = exp
                else:
                    d2[base] = exp                    
            elif type(exp) is cls and x in exp.symbols_data:
                d3[base] = exp
            else:
                d4[base] = exp
        if d1 or (d2 and d3) or (len(d2)>1) or (len(d3)>1):
            raise NotImplementedError("don't know how to integrate %s over %s" % (expr, x))
        if not (d2 or d3):
            return expr * cls(SYMBOL, x)            
        if d4:
            if len(d4)>1:
                const = cls(BASE_EXP_DICT, d4)
            else:
                const = pow_new(cls, dict_get_item(d4))
        else:
            const = 1
        if d2:
            newexpr = pow_new(cls, dict_get_item(d2))
            return newexpr.head.integrate_indefinite(cls, newexpr.data, newexpr, x) * const
        if d3:
            newexpr = pow_new(cls, dict_get_item(d3))
            return newexpr.head.integrate_indefinite(cls, newexpr.data, newexpr, x) * const
        raise NotImplementedError("don't know how to integrate %s over %s" % (expr, x))

    def integrate_definite(self, cls, data, expr, x, a, b):
        d1 = {} # f(x)**g(x)
        d2 = {} # f(x)**const
        d3 = {} # const**g(x)
        d4 = {} # const**const
        for base, exp in data.iteritems():
            if x in base.symbols_data:
                if type(exp) is cls and x in exp.symbols_data:
                    d1[base] = exp
                else:
                    d2[base] = exp                    
            elif type(exp) is cls and x in exp.symbols_data:
                d3[base] = exp
            else:
                d4[base] = exp
        if d1 or (d2 and d3) or (len(d2)>1) or (len(d3)>1):
            raise NotImplementedError("don't know how to integrate %s over %s in [%s, %s]" % (expr, x, a, b))
        if not (d2 or d3):
            return (b-a) * cls(SYMBOL, x)            
        if d4:
            if len(d4)>1:
                const = cls(BASE_EXP_DICT, d4)
            else:
                const = pow_new(cls, dict_get_item(d4))
        else:
            const = 1
        if d2:
            newexpr = pow_new(cls, dict_get_item(d2))
            return newexpr.head.integrate_definite(cls, newexpr.data, newexpr, x, a, b) * const
        if d3:
            newexpr = pow_new(cls, dict_get_item(d3))
            return newexpr.head.integrate_definite(cls, newexpr.data, newexpr, x, a, b) * const
        raise NotImplementedError("don't know how to integrate %s over %s in [%s, %s]" % (expr, x, a, b))
    
BASE_EXP_DICT = BaseExpDictHead()

