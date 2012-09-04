
__all__ = ['Head', 'UnaryHead', 'BinaryHead', 'NaryHead', 'HEAD']

not_implemented_error_msg = '%s.%s method\n  Please report this issue to http://code.google.com/p/sympycore/issues/ .'

from ..core import Expr, heads, heads_precedence, Pair, MetaCopyMethodsInheritDocs

from ..core import init_module
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

class NotImplementedHeadMethod(NotImplementedError):

    def __init__(self, head, signature):
        NotImplementedError.__init__(self, not_implemented_error_msg % (head, signature))

class Head(object):

    """
    Head - base class for expression heads.

    Recall that any expression is represented as a pair: Expr(head,
    data).  The head part defines how the data part should be
    interpreted.  Various operations on expressions are defined as
    Head methods taking the data part as an argument and returning the
    result of the operation.
    """

    """
    MetaCopyMethodsInheritDocs performs the following tasks:

      * If a class defines a method foo and sets ``rfoo = foo`` then
        make a copy of ``foo`` to ``rfoo``.

      * Append documentation strings of base class methods to class
        method documentation string.
    """
    __metaclass__ = MetaCopyMethodsInheritDocs
    
    """
    precedence_map defines the order of operations if they appear in
    the same expression. Operations with higher precedence are applied
    first. Operations with equal precedence are applied from left to
    right.

    For example, multiplication has higher precedence than addition.

    In the following, the precendence value is a floating point number
    in a range [0.0, 1.0]. Lowest precedence value 0.0 is assigned
    for atomic expressions.
    """
    precedence_map = dict(
        LAMBDA = 0.0, ASSIGN = 0.0,
        ARG = 0.01, KWARG = 0.01, COLON = 0.01, COMMA = 0.01,
        OR = 0.02, AND = 0.03, NOT = 0.04,
        LT = 0.1, LE = 0.1, GT = 0.1, GE = 0.1, EQ = 0.09, NE = 0.09,
        IN = 0.1, NOTIN = 0.1, IS = 0.1, ISNOT = 0.1,
        BOR = 0.2, BXOR = 0.21, BAND = 0.22, 
        LSHIFT = 0.3, RSHIFT = 0.3,
        
        ADD = 0.4, SUB = 0.4, TERM_COEFF_DICT = 0.4,

        TERM_COEFF = 0.5, NCMUL = 0.5, MUL = 0.5, DIV = 0.5,
        MOD = 0.5, FLOORDIV = 0.5,
        FACTORS = 0.5, BASE_EXP_DICT = 0.5,

        POS = 0.6, NEG = 0.6, INVERT = 0.6,
        POW = 0.7, POWPOW = 0.71,
        # a.foo(), a[]()
        APPLY = 0.8,
        ATTR = 0.81, SUBSCRIPT = 0.82, SLICE = 0.83, 
        TUPLE = 0.91, LIST = 0.92, DICT = 0.93,
        CALLABLE = 0.85,
        DOT = 0.9,
        SYMBOL = 1.0, NUMBER = 1.0, SPECIAL = 1.0
        )

    """
    op_mth and op_rmth contain the names of the corresponding Python
    left and right, respectively, operation methods.
    """
    op_mth = None
    op_rmth = None
    is_singleton = True

    """
    _cache contains unique instances of Head classes.
    """
    _cache = {}
    
    def __new__(cls, *args):
        """ Creates a Head instance.

        In most cases Head instances are singletons so that it can be
        efficiently compared using Python ``is`` operator.
        """
        if len(args)==1:
            arg = args[0]
            key = '%s(%s:%s)' % (cls.__name__, arg, id(arg))
        else:
            key = '%s%r' % (cls.__name__, args)
        obj = cls._cache.get(key)
        if obj is None:
            obj = object.__new__(cls)
            obj._key = key
            obj.init(*args)
            if cls.is_singleton:
                cls._cache[key] = obj
                setattr(heads, repr(obj), obj)
        return obj

    def init(self, *args):
        """ Initialize Head instance attributes.
        """
        pass #pragma NO COVER

    def as_unique_head(self):
        """
        Return unique instance of the head.

        The method is used by pickler support to make Head instances
        unique and ensure that unpickled heads are singletons.
        """
        return self._cache.get(self._key, self)

    def new(self, Algebra, data, evaluate=True):
        """
        Return a new Algebra expression instance containing data.  If
        evaluate is True then return canonized expression.
        """
        return Algebra(self, data)

    def reevaluate(self, Algebra, data):
        """
        Return reevaluated expression.
        """
        return self.new(Algebra, data, evaluate=True)

    def data_to_str_and_precedence(self, cls, data):
        """
        Return (<str>, <number>) where <str> is repr representation of
        an expression cls(self, data) and <number> is precedence
        number of the expression.
        """
        return '%s(%r, %r)' % (cls.__name__, self, data), 1.0

    def to_lowlevel(self, cls, data, pair):
        """
        Return a low-level representation of an expression pair.  It
        is used in object comparison and hash computation methods.
        """
        return pair

    def scan(self, proc, Algebra, data, target):
        """
        Apply proc function proc(Algebra, <head>, <data>, target) to the
        content of data and then to data.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'scan(proc, cls, data, target)')) #pragma NO COVER

    def walk(self, func, Algebra, data, target):
        """
        Apply func function func(Algebra, <head>, <data>, target) to
        the operands of Algebra(self, data) expression and form a new
        expression from the results of func calls.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'walk(func, cls, data, target)')) #pragma NO COVER        

    def is_data_ok(self, cls, data):
        """
        Check if data is valid object for current head.

        If data object is valid, return None.  If data object is
        invalid then return a string containing information which
        validity test failed.

        The ``is_data_ok`` method is called when creating expressions
        in pure Python mode and an error will be raised when data is
        invalid. When creating expressions in C mode then for
        efficiency it is assumed that data is always valid and
        ``is_data_ok`` method will be never called during expression
        creation.
        """
        return #pragma NO COVER

    def nonzero(self, cls, data):
        """
        Return truth value of relational expression.  If the truth
        value cannot be determined, return True,
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'nonzero(cls, data)')) #pragma NO COVER

    #def __todo_repr__(self):
    #    # TODO: undefined __repr__ should raise not implemented error
    #    raise NotImplementedError('Head subclass must implement __repr__ method returning singleton name')

    def base_exp(self, cls, expr):
        """
        Return (base, exponent) form of an expression cls(self, expr)
        so that

          base ** exponent == expr

        where exponent is a number.
        
        This method is used to collect powers of multiplication
        operands.
        """
        return expr, 1

    def term_coeff(self, cls, expr):
        """
        Return (term, coeff) form of an expression cls(self, expr)
        so that

          coeff * term == expr

        where coeff is a number.
        
        This method is used to collect coefficients of addition
        operands.
        """
        return expr, 1

    # We are not using ``rop_mth = op_mth`` because when derived class
    # redefines ``op_mth`` then it must also set ``rop_mth =
    # op_mth``. Not doing so may lead to bugs that are difficult to
    # track down. However, final derived classes may define ``rop_mth
    # = op_mth`` for efficiency.

    def neg(self, cls, expr):
        """
        Return the result of negation on given expression: -expr.
        """
        return cls(TERM_COEFF, (expr, -1))

    def add_number(self, cls, lhs, rhs):
        """
        Return the sum of expressions: lhs + rhs, where rhs is a number.
        """
        return cls(TERM_COEFF_DICT, {lhs:1, cls(NUMBER,1):rhs}) if rhs else lhs

    def add(self, cls, lhs, rhs):
        """ Return the sum of expressions: lhs + rhs.
        """
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            if rdata==0:
                return lhs
            return cls(TERM_COEFF_DICT, {lhs:1, cls(NUMBER,1):rdata})
        if rhead is self:
            if lhs==rhs:
                return cls(TERM_COEFF, (lhs, 2))
            return cls(TERM_COEFF_DICT, {lhs:1, rhs:1})
        if rhead is TERM_COEFF:
            t,c = rdata
            if lhs==t:
                return term_coeff_new(cls, (t, c+1))
            return cls(TERM_COEFF_DICT, {t:c, lhs:1})
        if rhead is TERM_COEFF_DICT:
            data = rdata.copy()
            term_coeff_dict_add_item(cls, data, lhs, 1)
            return term_coeff_dict_new(cls, data)            
        if rhead is SYMBOL or rhead is APPLY or rhead is CALLABLE\
               or rhead is BASE_EXP_DICT or rhead is POW or rhead is SUBSCRIPT:
            return cls(TERM_COEFF_DICT, {lhs:1, rhs:1})
        if rhead is MUL:
            return cls(ADD, [lhs, rhs])
        raise NotImplementedError(\
            not_implemented_error_msg % \
            (self, 'add(cls, <%s expression>, <%s expression>)' \
             % (self, rhs.head))) #pragma NO COVER

    def inplace_add(self, cls, lhs, rhs):
        """ Return the sum of expressions: lhs + rhs. If lhs is
        writable then lhs += rhs can be executed and lhs returned.
        """
        return self.add(cls, lhs, rhs)

    def sub_number(self, cls, lhs, rhs):
        """ Return the subtract of expressions: lhs - rhs, where rhs
        is number.
        """
        return cls(TERM_COEFF_DICT, {lhs:1, cls(NUMBER,1):-rhs}) if rhs else lhs

    def sub(self, cls, lhs, rhs):
        """ Return the subtract of expressions: lhs - rhs.
        """
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            if rdata==0:
                return lhs
            return cls(TERM_COEFF_DICT, {lhs:1, cls(NUMBER,1):-rdata})
        if rhead is self:
            if lhs==rhs:
                return cls(NUMBER, 0)
            return cls(TERM_COEFF_DICT, {lhs:1, rhs:-1})
        if rhead is TERM_COEFF:
            t,c = rdata
            if lhs==t:
                return term_coeff_new(cls, (t, 1-c))
            return cls(TERM_COEFF_DICT, {t:-c, lhs:1})
        if rhead is TERM_COEFF_DICT:
            data = rdata.copy()
            term_coeff_dict_mul_value(cls, data, -1)
            term_coeff_dict_add_item(cls, data, lhs, 1)
            return term_coeff_dict_new(cls, data)            
        if rhead is SYMBOL or rhead is APPLY or rhead is CALLABLE\
               or rhead is BASE_EXP_DICT or rhead is POW:
            return cls(TERM_COEFF_DICT, {lhs:1, rhs:-1})
        if rhead is MUL:
            return cls(ADD, [lhs, -rhs])
        raise NotImplementedError(\
            not_implemented_error_msg % \
            (self, 'sub(cls, <%s expression>, <%s expression>)' \
             % (self, rhs.head))) #pragma NO COVER

    def inplace_sub(self, cls, lhs, rhs):
        """
        Return the subtract of expressions: lhs - rhs. If lhs is
        writable then lhs -= rhs can be executed and lhs returned.
        """
        return self.sub(cls, lhs, rhs)

    def commutative_mul_number(self, cls, lhs, rhs):
        """
        Return the commutative product of expressions: lhs * rhs
        where rhs is a number.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'commutative_mul_number(Algebra, lhs, rhs)')) #pragma NO COVER
        return term_coeff_new(cls, (lhs, rhs))

    def commutative_rmul_number(self, cls, rhs, lhs):
        """
        Return the commutative product of expressions: lhs * rhs
        where lhs is a number.
        """
        return self.commutative_mul_number(cls, rhs, lhs)

    def commutative_mul(self, cls, lhs, rhs):
        """
        Return the commutative product of expressions: lhs * rhs.
        """
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return term_coeff_new(cls, (lhs, rdata))
        if rhead is self:
            if lhs.data==rdata:
                return cls(POW, (lhs, 2))
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return (lhs * term) * coeff
        if rhead is POW:
            rbase, rexp = rdata
            if rbase==lhs:
                return pow_new(cls, (lhs, rexp+1))
            return cls(BASE_EXP_DICT, {lhs:1, rbase:rexp})
        if rhead is BASE_EXP_DICT:
            data = rdata.copy()
            base_exp_dict_add_item(cls, data, lhs, 1)
            return base_exp_dict_new(cls, data)
        if rhead is SYMBOL or rhead is CALLABLE or rhead is APPLY \
           or rhead is TERM_COEFF_DICT or rhead is ADD or rhead is SUBSCRIPT:
            return cls(BASE_EXP_DICT, {lhs:1, rhs:1})
        if rhead is EXP_COEFF_DICT:
            return lhs * rhs.to(TERM_COEFF_DICT)
        raise NotImplementedError(\
            not_implemented_error_msg % \
            (self, 'commutative_mul(cls, <%s expression>, <%s expression>)' \
             % (self, rhs.head))) #pragma NO COVER

    def inplace_commutative_mul(self, cls, lhs, rhs):
        """
        Return the commutative product of expressions: lhs * rhs. If
        lhs is writable then lhs *= rhs can be executed and lhs
        returned.
        """
        return self.commutative_mul(cls, lhs, rhs)

    def commutative_div_number(self, cls, lhs, rhs):
        """
        Return the commutative division of expressions: lhs / rhs
        where rhs is a number.
        """
        r = number_div(cls, 1, rhs)
        if rhs==0:
            return r * lhs
        return self.commutative_mul_number(cls, lhs, r)

    def commutative_rdiv_number(self, cls, rhs, lhs):
        """
        Return the commutative division of expressions: lhs / rhs
        where lhs is a number.
        """
        # ensure that rhs is such that rhs ** -1 == cls(POW,(rhs,-1)).
        return term_coeff_new(cls, (cls(POW, (rhs, -1)), lhs))

    def commutative_div(self, cls, lhs, rhs):
        """
        Return the commutative division of expressions: lhs / rhs.
        """
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            r = number_div(cls, 1, rdata)
            if rdata==0:
                return r * lhs
            return term_coeff_new(cls, (lhs, r))
        if rhead is self:
            if lhs.data==rdata:
                return cls(NUMBER, 1)
            return cls(BASE_EXP_DICT, {lhs:1, rhs:-1})
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return number_div(cls, 1, coeff) * (lhs / term) 
        if rhead is POW:
            rbase, rexp = rdata
            if lhs==rbase:
                return pow_new(cls, (lhs, 1-rexp))
            return cls(BASE_EXP_DICT, {lhs:1, rbase:-rexp, })
        if rhead is BASE_EXP_DICT:
            data = {lhs:1}
            base_exp_dict_sub_dict(cls, data, rdata)
            return base_exp_dict_new(cls, data)
        if rhead is SYMBOL or rhead is CALLABLE or rhead is APPLY \
           or rhead is TERM_COEFF_DICT or head is ADD:
            return cls(BASE_EXP_DICT, {lhs:1, rhs:-1})
        
        raise NotImplementedError(\
            not_implemented_error_msg % \
            (self, 'commutative_div(cls, <%s expression>, <%s expression>)' \
             % (self, rhs.head))) #pragma NO COVER

    def non_commutative_mul_number(self, cls, lhs, rhs):
        """
        Return the non-commutative product of expressions: lhs * rhs
        where rhs is a number (which is assumed to be commutator).
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'non_commutative_mul_number(Algebra, lhs, rhs)')) #pragma NO COVER

    def non_commutative_rmul_number(self, cls, rhs, lhs):
        """
        Return the non-commutative product of expressions: lhs * rhs
        where rhs is a number (which is assumed to be commutator).
        """
        return self.non_commutative_mul_number(cls, rhs, lhs)

    def non_commutative_mul(self, cls, lhs, rhs):
        """
        Return the non-commutative product of expressions: lhs * rhs.
        """
        rhead, rdata = rhs.pair
        if rhead is NUMBER:
            return term_coeff_new(cls, (lhs, rdata))
        if rhead is self:
            if lhs.data == rdata:
                return cls(POW, (lhs, 2))
            return cls(MUL, [lhs, rhs])
        if rhead is TERM_COEFF:
            term, coeff = rdata
            return (lhs * term) * coeff
        if rhead is POW:
            return MUL.combine(cls, [lhs, rhs])
        if rhead is MUL:
            return MUL.combine(cls, [lhs] + rdata)
        raise NotImplementedError(\
            not_implemented_error_msg % \
            (self, 'non_commutative_mul(cls, <%s expression>, <%s expression>)' \
             % (self, rhs.head))) #pragma NO COVER

    def non_commutative_div_number(self, cls, lhs, rhs):
        """
        Return the non-commutative division of expressions: lhs / rhs
        where rhs is a number (which is assumed to be commutator).
        """
        r = number_div(cls, 1, rhs)
        if rhs==0:
            # lhs/0 -> zoo * lhs
            return r * lhs
        return self.non_commutative_mul_number(cls, lhs, r)

    def non_commutative_div(self, cls, lhs, rhs):
        """
        Return the non-commutative division of expressions: lhs / rhs.
        """
        return lhs * (rhs**-1)

    def pow_number(self, cls, base, exp):
        """
        Return the exponentiation: base ** exp, where exp is number.
        """
        return pow_new(cls, (base, exp))

    def pow(self, cls, base, exp):
        """
        Return the exponentiation: base ** exp.
        """
        return pow_new(cls, (base, exp))

    def diff(self, Algebra, data, expr, symbol, order):
        """
        Return the order-th derivative of expr with respect to symbol.
        data is expr.data.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'diff(Algebra, data, expr, symbol, order)')) #pragma NO COVER

    def fdiff(self, Algebra, data, expr, argument_index, order):
        """
        Return the order-th derivative of a function expr with respect
        to argument_index-th argument. data is expr.data.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'fdiff(Algebra, data, expr, argument_index, order)')) #pragma NO COVER

    def integrate_indefinite(self, Algebra, data, expr, x):
        """
        Return indefinite integral of expr with respect to x.
        data is expr.data, x is str object.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'integrate_indefinite(Algebra, data, expr, x)')) #pragma NO COVER

    def integrate_definite(self, Algebra, data, expr, x, a, b):
        """
        Return definite integral of expr with respect to x in the
        interval [a, b].  data is expr.data, x is str object.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'integrate_definite(Algebra, data, expr, x, a, b)')) #pragma NO COVER

    def apply(self, cls, data, func, args):
        """
        Return unevaluated function applied to arguments.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'apply(Algebra, data, func, args)')) #pragma NO COVER

    def diff_apply(self, cls, data, diff, expr):
        """
        Return unevaluated derivative applied to expr.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'diff_apply(Algebra, data, diff, expr)')) #pragma NO COVER
    
    def expand(self, Algebra, expr):
        """
        Return the expanded expression of expr, i.e. open parenthesis.
        """
        return expr

    def expand_intpow(self, Algebra, base, exp):
        """
        Return the expanded expression of base ** exp, where exp is
        integer.
        """
        return Algebra(POW, (base, exp))

    def to_ADD(self, Algebra, data, expr):
        """
        Convert expr to an expression with ADD head.
        data is expr.data.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'to_ADD(Algebra, data, expr)')) #pragma NO COVER

    def to_MUL(self, Algebra, data, expr):
        """
        Convert expr to an expression with MUL head.
        data is expr.data.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'to_MUL(Algebra, data, expr)')) #pragma NO COVER

    def to_EXP_COEFF_DICT(self, Algebra, data, expr, variables = None):
        """
        Convert expr to an expression with EXP_COEFF_DICT head.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'to_EXP_COEFF_DICT(Algebra, data, expr, variables=)')) #pragma NO COVER

    def to_TERM_COEFF_DICT(self, Algebra, data, expr):
        """
        Convert expr to an expression with TERM_COEFF_DICT head.
        data is expr.data. Note that the returned result may have
        actual head NUMBER, SYMBOL, TERM_COEFF, POW, BASE_EXP_DICT
        instead of TERM_COEFF_DICT.
        """
        raise NotImplementedError(not_implemented_error_msg % (self, 'to_TERM_COEFF_DICT(Algebra, data, expr)')) #pragma NO COVER

    def algebra_pos(self, Algebra, expr):
        """
        Return the position of an expression: +expr
        """
        raise NotImplementedHeadMethod(self, "algebra_pos(Algebra, expr)") #pragma NO COVER

    def algebra_neg(self, Algebra, expr):
        """
        Return the negation of an expression: -expr
        """
        raise NotImplementedHeadMethod(self, "algebra_neg(Algebra, expr)") #pragma NO COVER

    def algebra_add_number(self, Algebra, lhs, rhs, inplace):
        """
        Return the sum of expressions: lhs + rhs, where rhs is a number.
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_add_number(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

    def algebra_add(self, Algebra, lhs, rhs, inplace):
        """
        Return the sum of expressions: lhs + rhs.
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_add(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

    def algebra_mul_number(self, Algebra, lhs, rhs, inplace):
        """
        Return the product of expressions: lhs * rhs, where rhs is a number.
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_mul_number(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

    def algebra_mul(self, Algebra, lhs, rhs, inplace):
        """
        Return the product of expressions: lhs * rhs.
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_mul(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

    def algebra_div_number(self, Algebra, lhs, rhs, inplace):
        """
        Return the division of expressions: lhs / rhs, where rhs is a number.
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_div_number(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

    def algebra_div(self, Algebra, lhs, rhs, inplace):
        """
        Return the division of expressions: lhs / rhs.
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_div(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

    def algebra_pow_number(self, Algebra, lhs, rhs, inplace):
        """
        Return the exponentiation of expressions: lhs ** rhs, where rhs is a number
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_pow_number(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

    def algebra_pow(self, Algebra, lhs, rhs, inplace):
        """
        Return the exponentiation of expressions: lhs ** rhs.
        """
        t = getattr(rhs, 'head', type(rhs).__name__)
        raise NotImplementedHeadMethod(self, "algebra_pow(Algebra, lhs, rhs, inplace)<=%s" % (t)) #pragma NO COVER

class AtomicHead(Head):
    """
    AtomicHead is a base class to atomic expression heads.
    """

    def reevaluate(self, cls, data):
        return cls(self, data)

    def scan(self, proc, cls, data, target):
        proc(cls, self, data, target)

    def walk(self, func, cls, data, target):
        return func(cls, self, data, target)

class UnaryHead(Head):
    """
    UnaryHead is base class for unary operation heads,
    data is an expression operand.
    """

    # Derived class must define a string member:
    op_symbol = None
    
    def data_to_str_and_precedence(self, cls, operand):
        u_p = getattr(heads_precedence, repr(self))
        o, o_p = operand.head.data_to_str_and_precedence(cls, operand.data)
        if o_p < u_p: o = '(' + o + ')'
        return self.op_symbol + o, u_p

    def scan(self, proc, cls, expr, target):
        expr.head.scan(proc, cls, expr.data, target)
        proc(cls, self, expr, target)

    def walk(self, func, cls, operand, target):
        operand1 = operand.head.walk(func, cls, operand.data, operand)
        if operand1 is operand:
            return func(cls, self, operand, target)
        r = self.new(cls, operand1)
        return func(cls, r.head, r.data, r)
        
class BinaryHead(Head):
    """
    BinaryHead is base class for binary operation heads,
    data is a 2-tuple of expression operands.
    """
    def data_to_str_and_precedence(self, cls, (lhs, rhs)):
        rel_p = getattr(heads_precedence, repr(self))
        if isinstance(lhs, Expr):
            l, l_p = lhs.head.data_to_str_and_precedence(cls, lhs.data)
        else:
            l, l_p = str(lhs), 0.0
        if isinstance(rhs, Expr):
            r, r_p = rhs.head.data_to_str_and_precedence(cls, rhs.data)
        else:
            r, r_p = str(rhs), 0.0
        if l_p < rel_p: l = '(' + l + ')'
        if r_p < rel_p: r = '(' + r + ')'
        return l + self.op_symbol + r, rel_p

    def reevaluate(self, cls, (lhs, rhs)):
        return self.new(cls, (lhs, rhs))

    def scan(self, proc, cls, data, target):
        lhs, rhs = data
        lhs.head.scan(proc, cls, lhs.data, target)
        rhs.head.scan(proc, cls, rhs.data, target)
        proc(cls, self, data, target)

    def walk(self, func, cls, data, target):
        lhs, rhs = data
        lhs1 = lhs.head.walk(func, cls, lhs.data, lhs)
        rhs1 = rhs.head.walk(func, cls, rhs.data, rhs)
        if lhs1 is lhs and rhs1 is rhs:
            return func(cls, data, target)
        r = self.new(cls, (lhs1, rhs1))
        return func(cls, r.head, r.data, r)

class NaryHead(Head):
    """
    NaryHead is base class for n-ary operation heads,
    data is a n-tuple of expression operands.
    """
    def new(self, cls, data, evaluate=True):
        return cls(self, data)

    def reevaluate(self, cls, operands):
        return self.new(cls, operands)
    
    def data_to_str_and_precedence(self, cls, operand_seq):
        op_p = getattr(heads_precedence, repr(self))
        l = []
        for operand in operand_seq:
            o, o_p = operand.head.data_to_str_and_precedence(cls, operand.data)
            if o_p < op_p: o = '(' + o + ')'
            l.append(o)
        return self.op_symbol.join(l), op_p

    def scan(self, proc, cls, operands, target):
        for operand in operands:
            operand.head.scan(proc, cls, operand.data, target)
        proc(cls, self, operands, target)

    def walk(self, func, cls, operands, target):
        l = []
        flag = False
        for operand in operands:
            o = operand.head.walk(func, cls, operand.data, operand)
            if o is not operand:
                flag = True
            l.append(o)
        if flag:
            r = self.new(cls, l)
            return func(cls, r.head, r.data, r)
        return func(cls, self, operands, target)

class ArithmeticHead(Head):
    """ Base class for heads representing arithmetic operations.
    """
    
for k, v in Head.precedence_map.items():
    setattr(heads_precedence, k, v)

HEAD = Head()
