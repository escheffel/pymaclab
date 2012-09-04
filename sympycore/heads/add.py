
__all__ = ['ADD']

from .base import heads_precedence, ArithmeticHead

from ..core import init_module
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

@init_module
def _init(module):
    from ..arithmetic.number_theory import multinomial_coefficients
    module.multinomial_coefficients = multinomial_coefficients

class AddHead(ArithmeticHead):

    """
    AddHead represents addition n-ary operation where operands is
    given as a n-sequence of expressions. For example, expression 'a +
    2*b' is 'Expr(ADD, (a, 2*b))' where ADD=AddHead()
    """

    op_mth = '__add__'
    op_rmth = '__radd__'

    def is_data_ok(self, cls, data):
        if type(data) in [tuple, list]:
            for a in data:
                if not isinstance(a, cls):
                    return '%s data item must be %s instance but got %s' % (self, cls, type(a))
        else:
            return '%s data part must be a list but got %s' % (self, type(data))
    
    def __repr__(self): return 'ADD'

    def new(self, cls, operands, evaluate=True):
        if not evaluate:
            n = len(operands)
            if n==1:
                return operands[0]
            if n==0:
                return cls(NUMBER, 0)
            return cls(self, operands)
        d = {}
        l = []
        operands = list(operands)
        while operands:
            op = operands.pop(0)
            if op==0:
                continue
            head, data = op.pair
            if head is ADD:
                operands.extend(data)
                continue
            elif head is SUB:
                operands.append(data[0])
                for o in data[1:]:
                    operands.append(-o)
                continue
            elif head is TERM_COEFF_DICT:
                for term, coeff in data.iteritems():
                    n = len(d)
                    dict_add_item(cls, d, term, coeff)
                    if n < len(d):
                        l.append(term)
                    elif n > len(d):
                        l.remove(term)
            else:
                term, coeff = op.head.term_coeff(cls, op)
                n = len(d)
                dict_add_item(cls, d, term, coeff)
                if n < len(d):
                    l.append(term)
                elif n > len(d):
                    l.remove(term)
        r = []
        one = cls(NUMBER, 1)
        for term in l:
            r.append(TERM_COEFF.new(cls, (term, d[term])))
        m = len(r)
        if m==0:
            return cls(NUMBER, 0)
        if m==1:
            return r[0]
        return cls(self, r)

    def reevaluate(self, cls, operands):
        r = cls(NUMBER, 0)
        for op in operands:
            r += op
        return r

    def data_to_str_and_precedence(self, cls, operands):
        m = len(operands)
        if m==0:
            return '0', heads_precedence.NUMBER
        if m==1:
            op = operands[0]
            return op.head.data_to_str_and_precedence(cls, op.data)
        add_p = heads_precedence.ADD
        r = ''
        evaluate_addition = cls.algebra_options.get('evaluate_addition')
        for op in operands:
            t,t_p = op.head.data_to_str_and_precedence(cls, op.data)
            if not r:
                r += '(' + t + ')' if t_p < add_p else t
            elif evaluate_addition and t.startswith('-'):
                r += ' - ' + t[1:]
            else:
                r += ' + (' + t + ')' if t_p < add_p else ' + ' + t
        return r, add_p

    def term_coeff(self, cls, expr):
        term_list = expr.data
        if not term_list:
            return cls(NUMBER, 0), 1
        if len(term_list)==1:
            expr = term_list[0]
            return expr.head.term_coeff(cls, expr)
        return expr, 1

    def neg(self, cls, expr):
        return cls(ADD, [-term for term in expr.data])

    def add(self, cls, lhs, rhs):
        term_list = lhs.data
        rhead, rdata = rhs.pair
        if rhead is ADD:
            return ADD.new(cls, lhs.data + rdata)
        if rhead is TERM_COEFF_DICT:
            rdata = [t * c for t, c in rdata.items()]
            return ADD.new(cls, lhs.data + rdata)
        if rhead is SUB:
            rdata = rdata[:1] + [-op for op in rdata[1:]]
            return ADD.new(cls, lhs.data + rdata)
        return ADD.new(cls, lhs.data + [rhs])

    inplace_add = add

    def sub(self, cls, lhs, rhs):
        return lhs + (-rhs)

    def commutative_mul(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return ADD.new(cls, [op*rhs for op in lhs.data])
        if rhead is SYMBOL:
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        if rhead is ADD:
            if lhs==rhs:
                return lhs ** 2
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return (lhs * term) * coeff
        if rhead is POW:
            base, exp = rhs.data
            if lhs==base:
                return POW.new(cls, (base, exp + 1))
            return cls(BASE_EXP_DICT, {lhs:1, base:exp})
        if rhead is BASE_EXP_DICT:
            data = rdata.copy()
            dict_add_item(cls, data, lhs, 1)
            return BASE_EXP_DICT.new(cls, data)

        raise NotImplementedError(`self, lhs.pair, rhs.pair`)

    def commutative_mul_number(self, cls, lhs, rhs):
        return ADD.new(cls, [op*rhs for op in lhs.data])

    def pow(self, cls, base, exp):
        return POW.new(cls, (base, exp))

    pow_number = pow

    def walk(self, func, cls, data, target):
        l = []
        flag = False
        for op in data:
            o = op.head.walk(func, cls, op.data, op)
            if op is not o:
                flag = True
            l.append(o)
        if flag:
            r = ADD.new(cls, l)
            return func(cls, r.head, r.data, r)
        return func(cls, self, data, target)

    def scan(self, proc, cls, operands, target):
        for operand in operands:
            operand.head.scan(proc, cls, operand.data, target)
        proc(cls, self, operands, target)

    def expand(self, cls, expr):
        l = []
        for op in expr.data:
            h, d = op.pair
            l.append(h.expand(cls, op))
        return self.new(cls, l)

    def expand_intpow(self, cls, expr, intexp):
        if intexp<=1:
            return POW.new(cls, (expr, intexp))
        operands = expr.data
        mdata = multinomial_coefficients(len(operands), intexp)
        s = cls(NUMBER, 0)
        for exps, n in mdata.iteritems():
            m = cls(NUMBER, n)
            for i,e in enumerate(exps):
                m *= operands[i] ** e
            s += m
        return s

    def to_TERM_COEFF_DICT(self, Algebra, data, expr):
        s = Algebra(NUMBER, 0)
        for op in data:
            s += op.head.to_TERM_COEFF_DICT(Algebra, op.data, op)
        return s

    def to_ADD(self, Algebra, data, expr):
        return expr

    def algebra_pos(self, Algebra, expr):
        if Algebra.algebra_options.get('evaluate_addition'):
            if Algebra.algebra_options.get('is_additive_group_commutative'):
                return +ADD.to_TERM_COEFF_DICT(Algebra, expr.data, expr)
        return expr

    def algebra_neg(self, Algebra, expr):
        if Algebra.algebra_options.get('evaluate_addition'):
            if Algebra.algebra_options.get('is_additive_group_commutative'):
                return -ADD.to_TERM_COEFF_DICT(Algebra, expr.data, expr)
            return add_new(Algebra, [-op for op in expr.data[::-1]])
        return Algebra(NEG, expr)

    def combine_add_list(self, Algebra, data):
        """
        Combine add operands of an additive group in data.
        data will be changed in place.
        """
        commutative = Algebra.algebra_options.get('is_additive_group_commutative')
        if commutative:
            d = {}
            for op in data:
                term, coeff = op.head.term_coeff(Algebra, op)
                term_coeff_dict_add_item(Algebra, d, term, coeff)
            data[:] = [term_coeff_new(Algebra, term_coeff) for term_coeff in d.iteritems()]
        else:
            n = len(data)
            i0 = 0
            while 1:
                i = i0
                if i+1 >= n:
                    break
                lhs = data[i]
                rhs = data[i+1]
                lterm, lcoeff = lhs.head.term_coeff(Algebra, lhs)
                rterm, rcoeff = rhs.head.term_coeff(Algebra, rhs)
                if lterm==rterm:
                    coeff = lcoeff + rcoeff
                    if coeff:
                        del data[i+1]
                        data[i] = term_coeff_new(Algebra, (lterm, coeff))
                        i0 = i
                        n -= 1
                    else:
                        del data[i:i+2]
                        i0 = max(i - 1, 0)
                        n -= 2
                elif not rcoeff:
                    del data[i+1]
                    n -= 1
                    i0 = i
                elif not lcoeff:
                    del data[i]
                    n -= 1
                    i0 = max(i-1,0)
                else:
                    i0 += 1
        return data

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        return self.algebra_add(Algebra, lhs, Algebra(NUMBER, rhs), inplace)
        
    def algebra_add(self, Algebra, lhs, rhs, inplace):
        rhead, rdata = rhs.pair
        if rhead is TERM_COEFF_DICT or rhead is EXP_COEFF_DICT or rhead is MUL or rhead is NEG:
            rhs = rhs.to(ADD)
            rhead, rdata = rhs.pair
        if inplace:
            data = lhs.data
        else:
            data = lhs.data[:]
        if rhead is ADD:
            data.extend(rdata)
        else:
            data.append(rhs)
        if Algebra.algebra_options.get('evaluate_addition'):
            self.combine_add_list(Algebra, data)
        if inplace:
            return add(Algebra, lhs)
        return add_new(Algebra, data)

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        ntype = type(rhs)
        if Algebra.algebra_options.get('is_additive_group_commutative'):
            if not rhs:
                return Algebra(NUMBER, 0)
            return ADD.to_TERM_COEFF_DICT(Algebra, lhs.data, lhs) * rhs
        else:
            if Algebra.algebra_options.get('evaluate_addition'):
                if rhs == 0:
                    return Algebra(NUMBER, 0)
                if rhs == 1:
                    return lhs
                if ntype in inttypes_set:
                    if rhs > 0:
                        # (x+y)*3 = x+y+x+y+x+y
                        # TODO: optimize (x+y+x)*3 = x+y+x+x+y++x+x+y+x = x+y+2*x+y+2*x+y+x
                        data = lhs.data * rhs
                    else:
                        data = [-op for op in (lhs.data * (-rhs))[::-1]]
                    self.combine_add_list(Algebra, data)
                    return add_new(Algebra, data)
                    
            return mul_new(Algebra, [lhs, Algebra(NUMBER, rhs)])

    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        ldata = lhs.data
        if Algebra.algebra_options.get('is_additive_group_commutative'):
            return ADD.to_TERM_COEFF_DICT(Algebra, lhs.data, lhs) * rhs
        else:
            if Algebra.algebra_options.get('evaluate_addition'):
                rhead, rdata = rhs.pair
                if rhead is NUMBER:
                    return ADD.algebra_mul_number(Algebra, lhs, rdata, inplace)
                return super(type(self), self).algebra_mul(Algebra, lhs, rhs, inplace)
            return mul_new(Algebra, [lhs, rhs])
    
ADD = AddHead()
