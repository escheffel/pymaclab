
__all__ = ['EXP_COEFF_DICT']

from .base import ArithmeticHead

from ..core import init_module, Pair, Expr
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

@init_module
def _init(module):
    from ..arithmetic.number_theory import multinomial_coefficients
    module.multinomial_coefficients = multinomial_coefficients

class ExpCoeffDict(ArithmeticHead):
    """
    """

    def is_data_ok(self, cls, data):
        if type(data) is not Pair:
            return 'data must be Pair instance but got %r' % (type(data).__name__)
        variables, exps_coeff_dict = data.pair
        if type(variables) is not tuple:
            return 'data[0] must be tuple but got %r' % (type(variables).__name__)
        if type(exps_coeff_dict) is not dict:
            return 'data[1] must be dict but got %r' % (type(exps_coeff_dict).__name__)
        for exps, coeff in exps_coeff_dict.items():
            if not coeff:
                return 'data[1] contains exp-zero-coeff pair: (%s, %s)' % (exps, coeff)
            if type(exps) is not IntegerList:
                del exps_coeff_dict[exps]
                exps = IntegerList(exps)
                exps_coeff_dict[exps] = coeff
                return 'data[1] keys must be IntegerList instances but got %r' % (type(exps).__name__)
        
    def __repr__(self): return 'EXP_COEFF_DICT'

    def data_to_str_and_precedence(self, cls, data):
        variables, exp_coeff_dict = data
        terms = []
        for exps in sorted(exp_coeff_dict, reverse=True):
            coeff = exp_coeff_dict[exps]
            if type(exps) is not IntegerList:
                # temporary hook for SPARSE_POLY head, remove if block when SPARSE_POLY is gone
                exps = IntegerList(exps)
            factors = []
            for var,exp in zip(variables, exps):
                if not isinstance(var, cls):
                    var = cls(SYMBOL, var)
                factors.append(cls(POW,(var, exp)))
            terms.append(cls(TERM_COEFF, (cls(MUL,factors), coeff)))
        return ADD.data_to_str_and_precedence(cls, terms)

    def to_lowlevel(self, cls, data, pair):
        variables, exp_coeff_dict = data.pair
        n = len(exp_coeff_dict)
        if n==0:
            return 0
        if n==1:
            exps, coeff = dict_get_item(exp_coeff_dict)
            if type(exps) is not IntegerList:
                # temporary hook for SPARSE_POLY head, remove if block when SPARSE_POLY is gone
                exps = IntegerList(exps)
            factors = []
            for var, exp in zip(variables, exps.data):
                if exp==0:
                    continue
                if not isinstance(var, cls):
                    var = cls(SYMBOL, var)
                if exp==1:
                    factors.append(var)
                else:
                    factors.append(cls(POW, (var, exp)))
            if not factors:
                return coeff
            term = MUL.new(cls, factors)
            return term_coeff_new(cls, (term, coeff))
        return pair

    def combine_variables(self, *seq):
        """
        Return a tuple of sorted variables combining given variables.
        """
        variables = set([])
        seq = list(seq)
        while seq:
            s = seq.pop(0)
            if isinstance(s, Expr):
                if s.head is SYMBOL:
                    variables.add(s.data)
                else:
                    variables.add(s)
            elif isinstance(s, str):
                variables.add(s)
            elif isinstance(s, (tuple, list)):
                seq.extend(s)
            elif s is None:
                pass
            else:
                raise TypeError('expected an expression or a sequence of expressions but got %n' % (type(s)))
        return tuple(sorted(variables))

    def make_exponent(self, expr, variables):
        """
        Return exponent list such that expr == variables ** exp_list.
        """
        if expr is 0:
            # shortcut to return exponent of a constant expr
            return IntegerList([0]*len(variables))
        i = list(variables).index(expr)
        if i==-1:
            exp = [0] * len(variables)
        else:
            exp = [0] * (i) + [1] + [0] * (len(variables)-i-1)
        return IntegerList(exp)

    def to_TERM_COEFF_DICT(self, Algebra, data, expr):
        variables, exps_coeff_dict = data
        d = {}
        for exps, coeff in exps_coeff_dict.iteritems():
            d1 = {}
            for exp, var in zip(exps, variables):
                if exp:
                    var = Algebra(var)
                    d1[var] = exp
            term = base_exp_dict_new(Algebra, d1)
            term_coeff_dict_add_item(Algebra, d, term, 1)
        return term_coeff_dict_new(Algebra, d)
        
    def to_EXP_COEFF_DICT(self, cls, data, expr, variables = None):
        if variables is None:
            return expr
        evars, edata = data.pair
        if evars==variables:
            return expr
        variables = self.combine_variables(evars, variables)
        levars = list(evars)
        l = []
        for v in variables:
            try:
                i = levars.index(v)
            except ValueError:
                i = None
            l.append(i)
        d = {}
        for exps, coeff in edata.iteritems():
            new_exps = IntegerList([(exps[i] if i is not None else 0) for i in l])
            d[new_exps] = coeff
        return cls(self, Pair(variables, d))

    def neg(self, cls, expr):
        return self.commutative_mul_number(cls, expr, -1)

    def add_number(self, cls, lhs, rhs, inplace=False):
        if not rhs:
            return lhs
        lvars, ldict = lhs.data.pair
        zero_exp = self.make_exponent(0, lvars)
        if inplace and lhs.is_writable:
            dict_add_item(cls, ldict, zero_exp, rhs)
            return lhs
        d = ldict.copy()
        dict_add_item(cls, d, zero_exp, rhs)
        return cls(self, Pair(lvars, d))

    def add(self, cls, lhs, rhs, inplace=False):
        lvars, ldict = lhs.data.pair
        rhead, rdata = rhs.pair
        if rhead is not EXP_COEFF_DICT:
            rhs = rhead.to_EXP_COEFF_DICT(cls, rdata, rhs, lvars)
            rhead, rdata = rhs.pair
        rvars, rdict = rdata.pair
        if lvars == rvars:
            if inplace and lhs.is_writable:
                dict_add_dict(cls, ldict, rdict)
                return lhs
            d = ldict.copy()
            dict_add_dict(cls, d, rdict)
            return cls(self, Pair(lvars, d))
        variables = tuple(sorted(set(lvars + rvars)))
        lhs = self.to_EXP_COEFF_DICT(cls, lhs.data, lhs, variables)
        rhs = self.to_EXP_COEFF_DICT(cls, rhs.data, rhs, variables)
        d = lhs.data.data.copy()
        dict_add_dict(cls, d, rhs.data.data)
        return cls(self, Pair(variables, d))

    def sub(self, cls, lhs, rhs):
        return self.add(cls, lhs, -rhs)

    def commutative_mul(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.commutative_mul_number(cls, lhs, rdata)
        lvars, ldict = lhs.data.pair
        if rhead is not EXP_COEFF_DICT:
            rhs = rhead.to_EXP_COEFF_DICT(cls, rdata, rhs, lvars)
            rhead, rdata = rhs.pair
        if rhead is EXP_COEFF_DICT:
            rvars, rdict = rdata.pair
            d = {}
            if lvars == rvars:
                exp_coeff_dict_mul_dict(cls, d, ldict, rdict)
                return cls(self, Pair(lvars, d))
            variables = tuple(sorted(set(lvars + rvars)))
            lhs = self.to_EXP_COEFF_DICT(cls, lhs.data, lhs, variables)
            rhs = self.to_EXP_COEFF_DICT(cls, rhs.data, rhs, variables)
            exp_coeff_dict_mul_dict(cls, d, lhs.data[1], rhs.data[1])
            return cls(self, Pair(variables, d))
        raise NotImplementedError(`self, rhs.head`)

    def commutative_mul_number(self, cls, lhs, rhs):
        lvars, ldict = lhs.data.pair
        if rhs==0:
            return cls(self, Pair(lvars, {}))
        if rhs==1:
            return lhs
        d = {}
        for exps, coeff in ldict.iteritems():
            d[exps] = coeff * rhs
        return cls(self, Pair(lvars, d))

    non_commutative_mul_number = commutative_mul_number

    def combine_ncmul_exponents(self, lexps, rexps, variables):
        """
        Return exponents of non-commutative multiplication:
          variables ** ([lexps] + [rexps]) -> variables ** exps
        TODO: move the algorithm to expr.py and implement its C version.
        """
        exps = list(lexps) + list(rexps)
        n = len(exps)
        i0 = 0
        while 1:
            i = i0
            while i < n and not exps[i]:
                i += 1
            if i==n:
                break
            j = i+1
            while j < n and not exps[j]:
                j += 1
            if j==n:
                break
            if variables[i] == variables[j]:
                exps[i] += exps[j]
                exps[j] = 0
                i0 = 0
            else:
                i0 = j
        return IntegerList(exps)

    def eliminate_trivial_exponents(self, Algebra, exps_coeff_dict, variables):
        """
        Eliminate common trivial exponents (=0) from the set of exponents
        and return the corresponding non-trivial list of variables.
        Note: exps_coeff_dict will be changed in place.
        """
        has_non_zeros = [0] * len(variables)
        for exps in exps_coeff_dict:
            n = 0
            for i, exp in enumerate(exps):
                if has_non_zeros[i]:
                    n += 1
                    continue
                if exp:
                    n += 1
                    has_non_zeros[i] = 1
            if n==len(variables): # no common trivial exponents
                return variables
        non_trivial_indices = [i for i in range(len(variables)) if has_non_zeros[i]]
        assert len(non_trivial_indices) < len(variables),`non_trivial_indices, variables`
        for exps in exps_coeff_dict.keys():
            coeff = exps_coeff_dict[exps]
            del exps_coeff_dict[exps]
            exps = IntegerList([exps[i] for i in non_trivial_indices])
            dict_add_item(Algebra, exps_coeff_dict, exps, coeff)
        return tuple([variables[i] for i in non_trivial_indices])
        
    def non_commutative_mul(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return self.non_commutative_mul_number(cls, lhs, rdata)
        lvars, ldict = lhs.data.pair
        if rhead is not EXP_COEFF_DICT:
            rhs = rhead.to_EXP_COEFF_DICT(cls, rdata, rhs, lvars)
            rhead, rdata = rhs.pair
        if rhead is EXP_COEFF_DICT:
            rvars, rdict = rdata.pair
            variables = lvars + rvars
            lhs = self.to_EXP_COEFF_DICT(cls, lhs.data, lhs, variables)
            rhs = self.to_EXP_COEFF_DICT(cls, rhs.data, rhs, variables)
            d = {}
            # TODO: move the following double-loop to expr.py and implement it in C.
            for lexp, lcoeff in ldict.iteritems():
                for rexp, rcoeff in rdict.iteritems():
                    exp = self.combine_ncmul_exponents(lexp, rexp, variables)
                    dict_add_item(cls, d, exp, lcoeff * rcoeff)
            variables = self.eliminate_trivial_exponents(cls, d, variables)
            return cls(self, Pair(variables, d))
        raise NotImplementedError(`self, rhs.head`)
    
    def commutative_div(self, cls, lhs, rhs):
        return lhs * rhs ** -1
    
    def pow(self, cls, base, exp):
        variables, exp_coeff_dict = base.data.pair
        if isinstance(exp, Expr) and exp.head is NUMBER and isinstance(exp.data, inttypes):
            exp = exp.data
        if isinstance(exp, inttypes):
            return self.pow_number(cls, base, exp)
        expr = cls(POW, (base, exp))
        variables = self.combine_variables(variables, expr)
        exps = self.make_exponent(expr, variables)
        return cls(self, Pair(variables, {exps:1}))

    def pow_number(self, cls, base, exp):
        variables, exp_coeff_dict = base.data.pair
        if isinstance(exp, inttypes):
            if exp==0:
                return cls(self, Pair(variables, {(0,)*len(variables):1}))
            if exp==1:
                return base
            if exp>1:
                exps_coeff_list = base.data.data.items()
                m = len(variables)
                mdata = multinomial_coefficients(len(exps_coeff_list), exp)
                d = {}
                for e,c in mdata.iteritems():
                    new_exps = IntegerList([0]*m)
                    new_coeff = c
                    for e_i, (exps,coeff) in zip(e, exps_coeff_list):
                        new_exps += exps * e_i
                        new_coeff *= coeff ** e_i
                    dict_add_item(cls, d, new_exps, new_coeff)
                return cls(self, Pair(variables, d))
            # exp is negative integer
            if len(exp_coeff_dict)==1:
                exps, coeff = dict_get_item(exp_coeff_dict)
                inv_coeff =  number_div(cls, 1, coeff)
                return cls(self, Pair(variables, {-exps: inv_coeff}))
        expr = cls(POW, (base, exp))
        variables = self.combine_variables(variables, expr)
        exps = self.make_exponent(expr, variables)
        return cls(self, Pair(variables, {exps:1}))

EXP_COEFF_DICT = ExpCoeffDict()
