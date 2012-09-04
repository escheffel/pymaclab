
__all__ = ['SPARSE_POLY', 'DENSE_POLY']

from .base import Head, heads

from ..core import init_module, Pair, Expr
init_module.import_heads()
init_module.import_numbers()

class SparsepolyHead(Head):
    """
    SparsepolyHead is a head for sparse polynomials represented as a
    dictionary of exponent and coefficient pairs, data can be dict or
    frozenset.
    """

    def __repr__(self): return 'SPARSE_POLY'

    def data_to_str_and_precedence(self, cls, data):
        return EXP_COEFF_DICT.data_to_str_and_precedence(cls, Pair(cls.variables, data))

    def reevaluate(self, cls, data):
        return cls(self, data)

    def to_lowlevel(self, cls, data, pair):
        return EXP_COEFF_DICT.to_lowlevel(cls, Pair(cls.variables, data), pair)

    def add(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is not self:
            rhs = rhead.to_SPARSE_POLY(cls, rdata, rhs)
            return lhs + rhs
        return NotImplemented

    inplace_add = add

    def sub_number(self, cls, lhs, rhs):
        return lhs + (-rhs)

    def sub(self, cls, lhs, rhs):
        return lhs + (-rhs)

    def commutative_mul(self, cls, lhs, rhs):
        rhead, rdata = rhs.pair
        if rhead is not self:
            rhs = rhead.to_SPARSE_POLY(cls, rdata, rhs)
            return lhs * rhs
        return NotImplemented

    inplace_commutative_mul = commutative_mul

    def term_coeff(self, cls, expr):
        return expr, 1

    def integrate_indefinite_index(self, cls, data, expr, index):
        """
        Return indefinite integral of expr with respect to index-th variable.
        data is expr.data, index is integer.
        """
        new_data = {}
        for exp, coeff in data.iteritems():
            new_exp = exp.copy()
            new_exp[index] += 1
            new_coeff = number_div(cls.ring, coeff, new_exp[index])
            new_data[new_exp] = new_coeff
        return cls(new_data)

    def diff_index(self, cls, data, expr, index, order):
        """
        Return the order-th derivative of expr with respect to index-th variable.
        data is expr.data, index is integer.
        """
        if order==0:
            return expr
        new_data = {}
        for exp, coeff in data.iteritems():
            if order > exp[index]:
                continue
            new_exp = exp.copy()
            new_exp[index] -= order
            new_coeff = coeff
            for i in range(order):
                new_coeff = new_coeff * (exp[index] - i)
            new_data[new_exp] = new_coeff
        return cls(new_data)

    def expand(self, cls, expr):
        new_data = {}
        for exp, coeff in expr.data.iteritems():
            if isinstance(coeff, Expr):
                new_data[exp] = coeff.expand()
            else:
                new_data[exp] = coeff
        return cls(new_data)

    def evalf(self, cls, expr, n):
        new_data = {}
        for exp, coeff in expr.data.iteritems():
            if isinstance(coeff, Expr):
                new_data[exp] = coeff.evalf(n)
            else:
                new_data[exp] = coeff
        return cls(new_data)
        

class DensepolyHead(Head):
    """
    DensepolyHead is a head for dense polynomials represented
    as n-tuple of coefficients, data is a 2-tuple (symbol, coeffseq).
    """
    def __repr__(self): return 'DENSE_POLY'

    def data_to_str_and_precedence(self, cls, (symbol, data)):
        if not isinstance(symbol, cls):
            symbol = cls(SYMBOL, symbol)
        terms = []
        for exp, coeff in enumerate(data):
            if coeff:
                terms.append(cls(TERM_COEFF, (cls(POW, (symbol, exp)), coeff)))
        return ADD.data_to_str_and_precedence(cls, terms)
                     
SPARSE_POLY = SparsepolyHead()
DENSE_POLY = DensepolyHead()
