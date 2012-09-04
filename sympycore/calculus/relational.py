#
# Created in February 2008 by Fredrik Johansson.
#
""" Provides some basic implementation of assumptions support.
"""

__docformat__ = "restructuredtext"
__all__ = ['Assumptions', 'is_positive']

from ..core import init_module
init_module.import_heads()
init_module.import_numbers()

from .algebra import Calculus
from ..logic import Logic, Lt, Le

@init_module
def init(module):
    module.no_assumptions = Assumptions([])
    module.globalctx = GlobalContext()
    module.globalctx.assumptions = module.no_assumptions
    from .functions import pi, E
    module.pi = pi
    module.E = E


def is_number(x):
    return isinstance(x, numbertypes) or (isinstance(x, Calculus) \
        and x.head is NUMBER)

class Assumptions:

    def __init__(self, statements=[]):
        self.pos_values = []
        self.nonneg_values = []
        for stmt in statements:
            if stmt is True:
                continue
            if stmt is False:
                raise ValueError("got False as assumption")
            if isinstance(stmt, Logic):
                head, (lhs, rhs) = stmt.pair
                if head is LT:
                    self.pos_values.append(rhs - lhs)
                elif head is LE:
                    self.nonneg_values.append(rhs - lhs)
                elif head is GT:
                    self.pos_values.append(lhs - rhs)
                elif head is GE:
                    self.nonneg_values.append(lhs - rhs)
                else:
                    raise ValueError("unknown assumption: " + repr(stmt))
            else:
                raise ValueError("unknown assumption type: " + repr(stmt))

    def __repr__(self):
        ps = [repr(Lt(0, a)) for a in self.pos_values]
        ns = [repr(Le(0, a)) for a in self.nonneg_values]
        return "Assumptions([%s])" % ", ".join(ps + ns)

    def __enter__(self):
        globalctx.assumptions = self

    def __exit__(self, *args):
        globalctx.assumptions = no_assumptions

    def check(self, cond):
        if isinstance(cond, (bool, type(None))):
            return cond
        if isinstance(cond, Logic):
            head, (lhs, rhs) = cond.pair
            if head is LT:
                return self.positive(rhs - lhs)
            if head is LE:
                return self.nonnegative(rhs - lhs)
            if head is GT:
                return self.positive(lhs - rhs)
            if head is GE:
                return self.nonnegative(lhs - rhs)
        raise ValueError(`cond`)

    def eq(s, a, b): return s.zero(a-b)
    def ne(s, a, b): return s.nonzero(a-b)
    def lt(s, a, b): return s.positive(b-a)
    def le(s, a, b): return s.nonnegative(b-a)
    def gt(s, a, b): return s.positive(a-b)
    def ge(s, a, b): return s.nonnegative(a-b)

    def negative(s, x):
        t = s.positive(x)
        if t is None:
            return t
        return not t

    def nonpositive(s, x):
        t = s.nonnegative(x)
        if t is None:
            return t
        return not t

    def zero(s, x):
        if is_number(x):
            return x == 0
        if s.positive(x) or s.negative(x):
            return False
        return None

    def nonzero(s, x):
        if is_number(x):
            return x != 0
        if s.positive(x) or s.negative(x):
            return True
        return None

    def positive(s, x):
        x = Calculus.convert(x)
        if x.head is NUMBER:
            val = x.data
            if isinstance(x.data, realtypes):
                return val > 0
        elif x.head is ADD:
            args = x.data
            if any(s.positive(a) for a in args) and all(s.nonnegative(a) for a in args): return True
            if any(s.negative(a) for a in args) and all(s.nonpositive(a) for a in args): return False
        elif x.head is TERM_COEFF:
            return s.positive(Calculus(MUL, map(Calculus, x.data)))
        elif x.head is TERM_COEFF_DICT:
            l = []
            for t,c in x.data.iteritems():
                l.append(t * c)
            return s.positive(Calculus(ADD, l))
        elif x.head is BASE_EXP_DICT:
            l = []
            for b, e in x.data.iteritems():
                l.append(b * e)
            return s.positive(Calculus(MUL, l))
        elif x.head is MUL:
            args = x.data
            if any(not s.nonzero(a) for a in args):
                return None
            neg = sum(s.negative(a) for a in args)
            return (neg % 2) == 0
        elif x.head is POW:
            b, e = x.data
            if s.positive(b) and s.positive(e):
                return True
        elif x == pi or x == E:
            return True
        if s.pos_values:
            # NOTE: this should check both x-p and x+p, i.e. bounds from both directions
            t1 = any(no_assumptions.nonnegative(x-p) for p in s.pos_values)
            t2 = any(no_assumptions.nonpositive(x-p) for p in s.pos_values)
            if t1 and not t2: return True
        return None

    def nonnegative(s, x):
        x = Calculus.convert(x)
        if x.head is NUMBER:
            val = x.data
            if isinstance(x.data, realtypes):
                return val >= 0
        elif x.head is ADD:
            args = x.data
            if all(s.nonnegative(a) for a in args): return True
            if all(s.negative(a) for a in args): return False
        elif x.head is MUL:
            args = x.data
            if all(s.nonnegative(a) for a in args): return True
        elif x.head is TERM_COEFF:
            return s.nonnegative(Calculus(MUL, map(Calculus, x.data)))            
        elif x.head is POW:
            b, e = x.data
            if s.nonnegative(b) and s.nonnegative(e):
                return True
        elif x.head is TERM_COEFF_DICT:
            l = []
            for t,c in x.data.iteritems():
                l.append(t * c)
            return s.nonnegative(Calculus(ADD, l))
        elif x.head is BASE_EXP_DICT:
            l = []
            for b, e in x.data.iteritems():
                l.append(b * e)
            return s.nonnegative(Calculus(MUL, l))
        elif x == pi or x == E:
            return True
        if s.nonneg_values:
            # NOTE: this should check both x-p and x+p, i.e. bounds from both directions
            t1 = any(no_assumptions.nonnegative(x-p) for p in s.nonneg_values)
            t2 = any(no_assumptions.negative(x-p) for p in s.nonneg_values)
            if t1 and not t2: return True
        return None


class GlobalContext(object):
    pass

def is_positive(e):
    return globalctx.assumptions.positive(e)
