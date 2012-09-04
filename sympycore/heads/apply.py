
__all__ = ['APPLY']

from .base import Head, heads_precedence
from .functional import FunctionalHead

from ..core import init_module, Expr
init_module.import_heads()
init_module.import_numbers()
init_module.import_lowlevel_operations()

class ApplyHead(FunctionalHead):
    """
    ApplyHead is a head for n-ary apply operation,
    data is a (1+n)-tuple of expression operands
    """
    op_mth = '__call__'
    parenthesis = '()'

    def is_data_ok(self, cls, data):
        if isinstance(data, tuple) and len(data)==2:
            func, args = data
            fcls = cls.get_function_algebra()
            if not isinstance(func, fcls):
                return '%s data[0] must be %s instance but got %s' % (self, fcls, type(func))
            msg = func.head.is_data_ok(type(func), func.data)
            if msg:
                return '%s data=%r: %s' % (func.head, func.pair, msg)
            if type(args) is tuple:
                for i,a in enumerate(args):
                    acls = func.get_argument_algebra(i)
                    if not isinstance(a, Expr) or type(a) is not acls:
                        return '%s data[1][%s] must be %s instance but got %s' % (self, i, acls, type(a))
            else:
                return '%s data[1] must be tuple but got %s' % (self, type(args))
        else:
            return '%s data instance must be 2-tuple' % (self)

    def __repr__(self): return 'APPLY'

    def new(self, cls, (func, args), evaluate=True):
        if not isinstance(func, Expr):
            fcls = cls.get_function_algebra()
            func = fcls(CALLABLE, func)
        return cls(APPLY, (func, args))

    def commutative_mul_number(self, cls, lhs, rhs):
        return term_coeff_new(cls, (lhs, rhs))

    def pow(self, cls, base, exp):
        return POW.new(cls, (base, exp))

    pow_number = pow

    def scan(self, proc, cls, data, target):
        f, args = data
        f.head.scan(proc, cls, f.data, target)
        for arg in args:
            arg.head.scan(proc, cls, arg.data, target)
        proc(cls, self, data, target)

    def walk(self, func, cls, data, target):
        f, args = data
        new_f = f.head.walk(func, cls, f.data, f)
        new_args = []
        flag = f is not new_f
        for arg in args:
            h, d = arg.pair
            new_arg = arg.head.walk(func, cls, arg.data, arg)
            if arg is not new_arg:
                flag = True
            new_args.append(new_arg)
        if flag:
            new_args = tuple(new_args)
            if new_f.head is CALLABLE:
                r = new_f.data(*new_args)
            else:
                r = cls(APPLY, (new_f, new_args))
            return func(cls, r.head, r.data, r)
        return func(cls, self, data, target)

    def expand(self, cls, expr):
        f, args = expr.data
        new_f = f
        new_args = [a.expand() for a in args]
        return new_f (*new_args)
        return expr

    def expand_intpow(self, cls, base, exp):
        return cls(POW, (base, exp))

    def diff(self, cls, data, expr, symbol, order, cache={}):
        key = (expr, symbol, order)
        result = cache.get(key)
        if result is not None:
            return result
        if isinstance(order, inttypes):
            pass            
        key1 = (expr, symbol, 1)
        result = cache.get(key1)
        if result is None:
            func, args = data
            fcls = cls.get_function_algebra()
            dcls = fcls.get_differential_algebra()
            if len(args)==1:
                arg = args[0]
                da = arg.head.diff(cls, arg.data, arg, symbol, 1, cache=cache)
                if symbol not in da.symbols_data:
                    # argument is linear with respect to symbol
                    df = func.head.fdiff(fcls, func.data, func, 0, order)
                    if df is NotImplemented:
                        if isinstance(order, int):
                            df = func.head.fdiff(fcls, func.data, func, 0, 1)
                            result = df(*args) * da
                            if order>1:
                                result = result.head.diff(cls, result.data, result, symbol, order-1, cache=cache)
                        else:
                            d = dcls(FDIFF, cls(NUMBER, 0))**order
                            df = fcls(APPLY, (d, (func,)))
                            result = df(*args) * da**order
                    else:
                        result = df(*args) * da**order
                elif isinstance(order, int):
                    df = func.head.fdiff(fcls, func.data, func, 0, 1)
                    result = df(*args) * da
                    if order>1:
                        result = result.head.diff(cls, result.data, result, symbol, order-1, cache=cache)
                else:
                    d = dcls(DIFF, cls(NUMBER, 0))**order
                    result = cls(APPLY, (d, (expr,)))
                cache[key] = result
                return result
            else:
                result = cls(NUMBER, 0)
                for i in range(len(args)):
                    arg = args[i]
                    da = arg.head.diff(cls, arg.data, arg, symbol, 1, cache=cache)
                    df = func.head.fdiff(fcls, func.data, func, i, 1)
                    result += df(*args) * da
            cache[key1] = result
        if order>1:
            result = result.head.diff(cls, result.data, result, symbol, order-1, cache=cache)
            cache[key] = result
        return result

    def apply(self, cls, data, func, args):
        return cls(APPLY, (func, args))

    def integrate_indefinite(self, cls, data, expr, x):
        if x not in expr.symbols_data:
            return expr * cls(SYMBOL, x)
        return cls(INTEGRAL_INDEFINITE, (expr, x))

    def integrate_definite(self, cls, data, expr, x, a, b):
        if x not in expr.symbols_data:
            return expr * (b - a)
        return cls(INTEGRAL_DEFINITE, (expr, x, a, b))

    

APPLY = ApplyHead()
