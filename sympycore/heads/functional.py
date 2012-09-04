
__all__ = ['SUBSCRIPT', 'SLICE', 'LAMBDA', 'ATTR', 'KWARG']

from .base import Head, heads_precedence

from ..core import init_module
init_module.import_heads()

class FunctionalHead(Head):

    # Derived classes must define `parenthesis` data member
    
    def data_to_str_and_precedence(self, cls, (func, args)):
        o_p = getattr(heads_precedence, repr(self))
        arg_p = heads_precedence.ARG        
        f, f_p = func.head.data_to_str_and_precedence(type(func), func.data)
        if f_p < o_p: f = '('+f+')'
        l = []
        assert type(args) in [tuple], `args`
        for arg in args:
            a, a_p = arg.head.data_to_str_and_precedence(cls, arg.data)
            if a_p < arg_p: a = '('+a+')'
            l.append(a)
        p1,p2 = self.parenthesis
        f = f + p1 + ', '.join(l) + p2
        return f,o_p


class SubscriptHead(FunctionalHead):
    """
    SubscriptHead is a head for n-ary subscript operation,
    data is a (1+n)-tuple of expression operands
    """

    op_mth = '__getitem__'
    parenthesis = '[]'

    def is_data_ok(self, cls, data):
        if type(data) is tuple:
            if len(data)<2:
                return '%s data part must be a 2-tuple but got %s-tuple' % (self, len(data))
            base, items = data
            if not isinstance(base, cls):
                return '%s data 1st item must be %s instance but got %s' % (self, cls, type(a))
            if type(items) is not tuple:
                return '%s data 2nd item must be a tuple but got %s' % (self, type(items))
            for i,a in enumerate(items):
                if not isinstance(a, cls):
                    return '%s data 2nd items %sth item must be %s instance but got %s' % (self, i, cls, type(a))
        else:
            return '%s data part must be a tuple but got %s' % (self, type(data))

    def __repr__(self): return 'SUBSCRIPT'

    def walk(self, func, cls, (base, items), target):
        new_base = base.head.walk(func, cls, base.data, base)
        flag = new_base is not base
        new_items = []
        for op in items:
            o = op.head.walk(func, cls, op.data, op)
            if op is not o:
                flag = True
            new_items.append(o)
        if flag:
            r = self.new(cls, (new_base, tuple(new_items)))
            return func(cls, r.head, r.data, r)
        return func(cls, self, (base, items), target)



class SliceHead(Head):
    """
    SliceHead is a head for 3-ary slice operation,
    data is a 3-tuple of expression operands
    """

    def data_to_str_and_precedence(self, cls, (start, stop, step)):
        slice_p = heads_precedence.SLICE
        colon_p = heads_precedence.COLON
        h, d = start.pair
        if h is SPECIAL and d is None:
            start = None
        else:
            start, p = h.data_to_str_and_precedence(cls, d)
            if p < colon_p: start = '('+start+')'
        h, d = stop.pair
        if h is SPECIAL and d is None:
            stop = None
        else:
            stop, p = h.data_to_str_and_precedence(cls, d)
            if p < colon_p: stop = '('+stop+')'
        h, d = step.pair
        if h is SPECIAL and d is None:
            step = None
        else:
            step, p = h.data_to_str_and_precedence(cls, d)
            if p < colon_p: step = '('+step+')'
        if start is None:
            if stop is None:
                if step is None: r = ':'
                else: r = '::%s' % step
            else:
                if step is None: r = ':%s' % stop
                else: r = ':%s:%s' % (stop, step)
        else:
            if stop is None:
                if step is None: r = '%s:' % start
                else: r = '%s::%s' % (start, step)
            else:
                if step is None: r ='%s:%s' % (start, stop)
                else: r = '%s:%s:%s' % (start, stop, step)
        return r, slice_p


    def __repr__(self): return 'SLICE'

class LambdaHead(Head):
    """
    LambdaHead is a head for lambda expression,
    data is a 2-tuple of arguments and an expression.
    """
    def data_to_str_and_precedence(self, cls, (args, expr)):
        arg_p = heads_precedence.ARG
        colon_p = heads_precedence.COLON
        l = []
        assert type(args) in [tuple],`args`
        for arg in args:
            a, a_p = arg.head.data_to_str_and_precedence(cls, arg.data)
            if a_p < arg_p: a = '('+a+')'
            l.append(a)
        e, e_p = expr.head.data_to_str_and_precedence(cls, expr.data)
        if e_p < colon_p: e = '('+e+')'
        return 'lambda %s: %s' % (', '.join(l), e), heads_precedence.LAMBDA

    def __repr__(self): return 'LAMBDA'

class AttrHead(Head):
    """
    AttrHead is a head for attribute operation,
    data is a 2-tuple of expression operands.
    """

    def data_to_str_and_precedence(self, cls, (expr, attr)):
        dot_p = heads_precedence.DOT
        e, e_p = expr.head.data_to_str_and_precedence(cls, expr.data)
        if e_p < dot_p: e = '('+e+')'
        a, a_p = attr.head.data_to_str_and_precedence(cls, attr.data)
        if a_p < dot_p: a = '('+a+')'
        return e + '.' + a, heads_precedence.ATTR
            
    def __repr__(self): return 'ATTR'

class KwargHead(Head):
    """
    Kwarg is a head for keyword argument,
    data is a 2-tuple of expression operands.
    """
    def data_to_str_and_precedence(self, cls, (name, expr)):
        kw_p = heads_precedence.KWARG
        assign_p = heads_precedence.ASSIGN
        n, n_p = name.head.data_to_str_and_precedence(cls, name.data)
        if n_p < assign_p: n = '('+n+')'
        e, e_p = expr.head.data_to_str_and_precedence(cls, expr.data)
        if e_p < assign_p: e = '('+e+')'
        return n + '=' + e, kw_p

    def __repr__(self): return 'KWARG'


SUBSCRIPT = SubscriptHead()
SLICE = SliceHead()
LAMBDA = LambdaHead()
ATTR = AttrHead()
KWARG = KwargHead()
