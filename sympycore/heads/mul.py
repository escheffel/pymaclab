
__all__ = ['MUL']

from .base import Head, heads_precedence, Pair, Expr, ArithmeticHead

from ..core import init_module
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

class MulHead(ArithmeticHead, Head):
    """
    Algebra(MUL, <list of factors>)

    The list of factors is should not contain numeric expressions
    representing a coefficient, numeric part should be represented
    via TERM_COEFF expression:
    
      Algebra(TERM_COEFF, (Algebra(MUL, <list of factors>), <numeric part>))

    """
    op_mth = '__mul__'
    op_rmth = '__rmul__'

    def is_data_ok(self, cls, data):
        if type(data) in [tuple, list]:
            for a in data:
                if not isinstance(a, cls):
                    return '%s data item must be %s instance but got %s' % (self, cls, type(a))
        else:
            return '%s data part must be a list but got %s' % (self, type(data))

    def __repr__(self): return 'MUL'

    def new(self, cls, operands, evaluate=True):
        operands = [op for op in operands if op!=1]
        n = len(operands)
        if n==0:
            return cls(NUMBER, 1)
        if n==1:
            return operands[0]
        return cls(MUL, operands)

    def reevaluate(self, cls, operands):
        r = operands[0]
        for op in operands[1:]:
            r *= op
        return r

    def base_exp(self, cls, expr):
        return expr, 1

    def data_to_str_and_precedence(self, cls, operands):
        m = len(operands)
        if m==0:
            return '1', heads_precedence.NUMBER
        if m==1:
            op = operands[0]
            return op.head.data_to_str_and_precedence(cls, op.data)
        mul_p = heads_precedence.MUL
        r = ''
        for op in operands:
            f, f_p = op.head.data_to_str_and_precedence(cls, op.data)
            if f=='1': continue
            if not r or r=='-':
                r += '('+f+')' if f_p<mul_p else f
            elif f.startswith('1/'):
                r += '*('+f+')' if f_p<mul_p else f[1:]
            else:
                r += '*('+f+')' if f_p<mul_p else '*'+f
        if not r:
            return '1', heads_precedence.NUMBER
        return r, mul_p

    def term_coeff(self, Algebra, expr):
        data = []
        coeff = 1
        for op in expr.data:
            t, c = op.head.term_coeff(Algebra, op)
            if c is not 1:
                coeff *= c
                if t.head is NUMBER:
                    assert t.data==1,`t`
                else:
                    data.append(t)
            else:
                data.append(op)
        if coeff is 1:
            return expr, 1
        return mul_new(Algebra, data), coeff

    def combine(self, cls, factors_list):
        """ Combine factors in a list and return result.
        """
        lst = []
        compart = 1
        for factor in factors_list:
            if factor.head is NUMBER:
                compart = compart * factor.data
                continue
            r = None
            b2, e2 = factor.head.base_exp(cls, factor)
            if lst and b2.head is MUL:
                l = b2.data
                if len(l)<=len(lst) and lst[-len(l):]==l:
                    # x*a*b*(a*b)**2 -> x*(a*b)**3
                    r = b2 ** (e2 + 1)
                    del lst[-len(l):]
            if lst and r is None:
                b1, e1 = lst[-1].head.base_exp(cls, lst[-1])
                if b1==b2:
                    # x*a**3*a**2 -> x*a**5
                    r = b2 ** (e1 + e2)
                    del lst[-1]
            if r is None:
                lst.append(factor)
                for i in range(2,len(lst)):
                    b1, e1 = lst[-i-1].head.base_exp(cls, lst[-i-1]);
                    if b1.head is MUL:
                        c, l = b1.data
                        if l == lst[-i:] and c==1:
                            # x*(a*b)**2*a*b -> x*(a*b)**3
                            r = b1 ** (e1 + 1)
                            del lst[-i-1:]
                            break
                    if lst[-i:]==lst[-2*i:-i]:
                        # x*a*b * a*b -> x*(a*b)**2
                        r = cls(MUL, lst[-i:])**2
                        del lst[-2*i:]
                        break
            if r is not None:
                if r.head is NUMBER:
                    compart = compart * r
                else:
                    lst.append(r)            
        if not lst:
            return compart
        if compart==1:
            if len(lst)==1:
                return lst[0]
            return cls(MUL, lst)
        return cls(TERM_COEFF, (cls(MUL, lst), compart))

    def non_commutative_mul(self, cls, lhs, rhs):
        head, data = rhs.pair
        if head is NUMBER:
            if data==1:
                return lhs
            return TERM_COEFF.new(cls, (lhs, data))
        if head is SYMBOL or head is POW:
            return self.combine(cls, lhs.data + [rhs])
        if head is TERM_COEFF:
            term, coeff = data
            return (lhs * term) * coeff
        if head is MUL:
            return self.combine(cls, lhs.data + rhs.data)
        raise NotImplementedError(`self, cls, lhs.pair, rhs.pair`)

    def non_commutative_mul_number(self, cls, lhs, rhs):
        return term_coeff_new(cls, (lhs, rhs))

    non_commutative_rmul_number = non_commutative_mul_number

    def pow(self, cls, base, exp):
        if exp==0:
            return cls(NUMBER, 1)
        if exp==1:
            return base
        if exp==-1:
            factors_list = [factor**-1 for factor in base.data]
            factors_list.reverse()
            return cls(MUL, factors_list)
        term, coeff = self.term_coeff(cls, base)
        if coeff!=1:
            return NUMBER.pow(cls, coeff, exp) * term**exp
        if isinstance(exp, Expr):
            h, d = exp.pair
            if h is NUMBER and isinstance(d, inttypes) and d>0:
                factors_list = base.data
                first = factors_list[0]
                last = factors_list[-1]
                a = last * first
                if a is not None and a.head is NUMBER: # todo: or a is commutative
                    compart = NUMBER.pow_number(cls, a, d)
                    rest = factors_list[1:-1]
                    if not rest:
                        return compart
                    if len(rest)==1:
                        middle = rest[0]
                    else:
                        middle = cls(MUL, rest)
                    return compart * first * middle**d * last # could be optimized
            if h is NUMBER:
                exp = d
        return cls(POW, (base, exp))

    def pow_number(self, cls, base, exp):
        return self.pow(cls, base, cls(NUMBER, exp))

    def walk(self, func, cls, data, target):
        l = []
        flag = False
        for op in data:
            o = op.head.walk(func, cls, op.data, op)
            if op is not o:
                flag = True
            l.append(o)
        if flag:
            r = MUL.new(cls, l)
            return func(cls, r.head, r.data, r)
        return func(cls, self, data, target)

    def scan(self, proc, cls, operands, target):
        for operand in operands:
            operand.head.scan(proc, cls, operand.data, target)
        proc(cls, self, operands, target)

    def combine_mul_list(self, Algebra, data):
        """
        Combine mul operands of an multiplicative group in data.
        data will be changed in place.
        """
        commutative = Algebra.algebra_options.get('is_multiplicative_group_commutative')
        coeff = 1
        if commutative:
            d = {}
            for op in data:
                base, exp = op.head.base_exp(Algebra, op)
                base_exp_dict_add_item(Algebra, d, base, exp)
            data[:] = [pow_new(Algebra, base_exp) for base_exp in d.iteritems()]
        else:
            n = len(data)
            i0 = 0
            while 1:
                i = i0
                if i+1 >= n:
                    break
                lhs = data[i]
                if lhs.head is NUMBER:
                    coeff *= lhs.data
                    del data[i]
                    n -= 1
                    continue
                rhs = data[i+1]
                if rhs.head is NUMBER:
                    coeff *= rhs.data
                    del data[i+1]
                    n -= 1
                    continue
                lbase, lexp = lhs.head.base_exp(Algebra, lhs)
                rbase, rexp = rhs.head.base_exp(Algebra, rhs)
                if lbase==rbase:
                    exp = lexp + rexp
                    if exp:
                        del data[i+1]
                        data[i] = pow_new(Algebra, (lbase, exp))
                        i0 = i
                        n -= 1
                    else:
                        del data[i:i+2]
                        i0 = max(i - 1, 0)
                        n -= 2
                elif not rexp:
                    del data[i+1]
                    n -= 1
                    i0 = i
                elif not lexp:
                    del data[i]
                    n -= 1
                    i0 = max(i-1,0)
                else:
                    i0 += 1
        return coeff

    def to_TERM_COEFF_DICT(self, Algebra, data, expr):
        m = data[0].to(TERM_COEFF_DICT)
        for op in data[1:]:
            m *= op.to(TERM_COEFF_DICT)
        return m

    def to_ADD(self, Algebra, data, expr):
        m = data[0].to(ADD)
        for op in data[1:]:
            m *= op.head.to_ADD(Algebra, op.data, op)
        return m

    def algebra_pos(self, Algebra, expr):
        return expr

    def algebra_neg(self, Algebra, expr):
        if Algebra.algebra_options.get('evaluate_addition'):
            return self.algebra_mul_number(Algebra, expr, -1, False)
        return Algebra(NEG, expr)

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        return self.algebra_add(Algebra, lhs, Algebra(NUMBER, rhs), inplace)

    def algebra_add(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_addition'):
            rhead, rdata = rhs.pair
            if rhead is TERM_COEFF_DICT or rhead is EXP_COEFF_DICT:
                rhs = rhs.head.to_ADD(Algebra, rdata, rhs)
                rhead, rdata = rhs.pair
            if rhead is ADD:
                data = [lhs] + rdata
            else:
                data = [lhs, rhs]
            ADD.combine_add_list(Algebra, data)
            return add_new(Algebra, data)
        return Algebra(ADD, [lhs, rhs])

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        return self.algebra_mul(Algebra, lhs, Algebra(NUMBER, rhs), inplace)

    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if rhead is BASE_EXP_DICT or rhead is TERM_COEFF:
            rhs = rhs.to(MUL)
            rhead, rdata = rhs.pair
        if inplace:
            data = lhs.data
        else:
            data = lhs.data[:]
        if rhead is MUL:
            data.extend(rdata)
        else:
            data.append(rhs)
        if Algebra.algebra_options.get('evaluate_multiplication'):
            coeff = self.combine_mul_list(Algebra, data)
            if not coeff:
                return Algebra(NUMBER, 0)
            if coeff != 1:
                if len(data)==1:
                    return Algebra(TERM_COEFF, (data[0], coeff))
                data.insert(0, Algebra(NUMBER, coeff))
        if inplace:
            return mul(Algebra, lhs)
        return mul_new(Algebra, data)

    def algebra_pow_number(self, Algebra, lhs, rhs, inplace):
        if Algebra.algebra_options.get('evaluate_multiplication'):
            if rhs==1:
                return lhs
            if not rhs:
                return Algebra(NUMBER, 1)
            if rhs==-1:
                return Algebra(MUL, [op**-1 for op in lhs.data[::-1]])
        return super(type(self), self).algebra_pow_number(Algebra, lhs, rhs, inplace)

MUL = MulHead()
