#
# Created January 2008 by Fredrik Johansson
#
""" Provides elementary calculus functions sqrt, exp, log, sin, etc and constants pi, E.
"""

__all__ = ['Sqrt', 'Exp', 'Log', 'Sin', 'Cos', 'Tan', 'Cot', 'Sign', 'Mod', 'Factorial',
           'E', 'pi', 'gamma', 'Ln', 'ArcSin']
__docformat__ = "restructuredtext"

from ...core import init_module
init_module.import_heads()
init_module.import_numbers()

from ..algebra import Calculus, I
from ..infinity import oo, zoo, undefined, CalculusInfinity
from ..constants import const_pi, const_E, const_gamma
from ...utils import str_SUM, str_SYMBOL

from ...core import DefinedFunction, get_nargs
from ...arithmetic.evalf import evalf
from ...arithmetic.numbers import complextypes, realtypes, inttypes
from ...arithmetic.number_theory import factorial
from ...arithmetic import infinity

from .algebra import CalculusFunctionRing

Apply = Calculus.Apply

import math

E = const_E.as_algebra(Calculus)
pi = const_pi.as_algebra(Calculus)
gamma = const_gamma.as_algebra(Calculus)

@init_module
def _init(m):
    m.zero = Calculus.zero
    m.one = Calculus.one
    m.half = one/2
    sqrt2 = Calculus.convert('2**(1/2)')
    sqrt3 = Calculus.convert('3**(1/2)')

    m.Ipi = Ipi = I*pi
    m.Ipi2 = Ipi2 = Ipi/2

    m.log_number_table = {
        m.zero.data : -oo,
        m.one.data : m.zero,
        I.data : Ipi2,
        (-I).data : -Ipi2
        }

    # Tabulated values of sin(x) at multiples of pi/12
    C0 = (sqrt3-1)/(2*sqrt2)
    C1 = one/2
    C2 = sqrt2/2
    C3 = sqrt3/2
    C4 = (sqrt3+1)/(2*sqrt2)
    
    # Replace entries with None to prevent from evaluating
    m.sine_table = [ \
  Calculus.zero, C0, C1, C2, C3, C4, Calculus.one, C4, C3, C2, C1,C0,
  Calculus.zero,-C0,-C1,-C2,-C3,-C4,-Calculus.one,-C4,-C3,-C2,-C1,-C0]

class CalculusDefinedFunction(DefinedFunction):

    @classmethod
    def get_argument_algebras(cls):
        return (Calculus,)*(get_nargs(cls))

class Sign(CalculusDefinedFunction):
    def __new__(cls, arg):
        if not isinstance(arg, Calculus):
            arg = Calculus.convert(arg)
        return Apply(cls, arg)

class Mod(CalculusDefinedFunction):
    def __new__(cls, x, y):
        if not isinstance(x, Calculus):
            x = Calculus.convert(x)
        if not isinstance(y, Calculus):
            y = Calculus.convert(y)
        xh,xd = x.pair
        yh,yd = y.pair
        if xh is NUMBER and yh is NUMBER:
            return Calculus.convert(xd % yd)
        return Apply(cls, (x, y))

class Factorial(CalculusDefinedFunction):
    
    def __new__(cls, x):
        if not isinstance(x, Calculus):
            x = Calculus.convert(x)
        if x.head is NUMBER:
            n = x.data
            if isinstance(n, (int, long)):
                if n >= 0:
                    return Calculus.convert(factorial(n))
                return zoo
        return Apply(cls, x)

#---------------------------------------------------------------------------#
#                                  Exponentials                             #
#---------------------------------------------------------------------------#

class Sqrt(CalculusDefinedFunction):
    def __new__(cls, arg):
        return arg ** half

class Exp(CalculusDefinedFunction):
    def __new__(cls, arg):
        return E ** arg

    @staticmethod
    def derivative(arg):
        return E ** arg

    @classmethod
    def fdiff(cls, FAlgebra, argument_index, order):
        if argument_index!=0:
            if isinstance(argument_index, inttypes):
                raise TypeError('%s.fdiff argument_index must be 0 but got %r' % (cls.__name__, argument_index))
            return NotImplemented
        return FAlgebra(CALLABLE, cls)

class Log(CalculusDefinedFunction):
    def __new__(cls, arg, base=E):
        if type(arg) is not Calculus:
            if isinstance(arg, CalculusInfinity):
                if arg == oo:
                    return oo
                if arg == undefined:
                    return undefined
                return Apply(cls, arg)
            else:
                arg = Calculus.convert(arg)
        if base != E:
            base = Calculus.convert(base)
            bd = base.data
            ad = arg.data
            if base.head is NUMBER and isinstance(bd, inttypes) and \
                arg.head is NUMBER and isinstance(ad, inttypes) and \
                ad > 0 and bd > 1:
                l = int(math.log(ad, bd) + 0.5)
                if bd**l == ad:
                    return Calculus.convert(l)
            return cls(arg) / cls(base)
        head = arg.head
        data = arg.data
        if head is NUMBER:
            v = log_number_table.get(data)
            if v is not None:
                return v
            if isinstance(data, realtypes) and data < 0:
                return Ipi + Log(-arg)
            if isinstance(data, complextypes) and data.real == 0:
                im = data.imag
                if im > 0: return Apply(cls, Calculus(NUMBER, im)) + Ipi2
                if im < 0: return Apply(cls, Calculus(NUMBER, -im)) - Ipi2
            return Apply(cls, arg)
        if arg == E:
            return one
        from ..relational import is_positive
        if head is POW:
            base, expt = data
            if is_positive(base) and isinstance(expt, realtypes):
                return Apply(cls, base) * expt
        if head is TERM_COEFF:
            term, coeff = data
            if (isinstance(coeff, realtypes) and coeff < 0) and is_positive(base):
                return Ipi + Log(-arg)
        return Apply(cls, arg)

    @classmethod
    def derivative(cls, arg):
        return 1/arg

    @classmethod
    def nth_derivative(cls, arg, n=1):
        return (-1)**(n-1) * Factorial(n-1) * arg**(-n)

    @classmethod
    def fdiff(cls, FAlgebra, argument_index, order):
        if argument_index!=0:
            if isinstance(argument_index, inttypes):
                raise TypeError('%s.fdiff argument_index must be 0 but got %r' % (cls.__name__, argument_index))
            return NotImplemented
        if order==1:
            func = lambda x: 1/x
        else:
            func = lambda x: (-1)**(order-1) * Factorial(order-1) * x**(-order)
        return FAlgebra(CALLABLE, func)

class Ln(Log):
    def __new__(cls, arg):
        return Log(arg)

#---------------------------------------------------------------------------#
#                          Trigonometric functions                          #
#---------------------------------------------------------------------------#



def get_pi_shift(arg, N):
    """Parse as x, n where arg = x + n*pi/N"""
    if arg.head is TERM_COEFF:
        e, c = arg.data
        if e == pi:
            c *= N
            if isinstance(c, inttypes):
                return zero, c
    elif arg.head is TERM_COEFF_DICT:
        c = arg.data.get(pi)
        if c is not None:
            c *= N
            if isinstance(c, inttypes):
                d = arg.data.copy()
                d.pop(pi)
                return TERM_COEFF_DICT.new(type(arg), d), c
    elif arg.head is ADD:
        rest = []
        x1 = None
        for op in arg.data:
            x, n = get_pi_shift(op, N)
            if n:
                x1, n1 = x, n
                assert x1==0,`x1,n1`
            else:
                rest.append(op)
        if x1 is not None:
            return ADD.new(type(arg), rest), n1
    elif arg == pi:
        return zero, N
    return arg, 0

def is_odd (expr):
    if isinstance(expr, int) and expr % 2:
        return True
    return False

def has_leading_sign(arg):
    """Detect symmetry cases for odd/even functions."""
    if arg.head is NUMBER:
        if arg < 0:
            return True
    if arg.head is TERM_COEFF:
        if arg.data[1] < 0:
            return True
        return has_leading_sign(arg.data[0])
    if arg.head is TERM_COEFF_DICT:
        for t,c in arg.data.iteritems():
            if not c < 0:
                return False
        return True
    if arg.head is BASE_EXP_DICT:
        for b,e in arg.data.iteritems():
            if is_odd(e):
                if has_leading_sign(b):
                    return True
    return None

class TrigonometricFunction(CalculusDefinedFunction):

    parity = None   # 'even' or 'odd'
    period = None   # multiple of pi

    def __new__(cls, arg):
        if type(arg) is not Calculus:
            if isinstance(arg, CalculusInfinity):
                if arg == undefined:
                    return undefined
                arg = Calculus(SPECIAL, arg)
                return Apply(cls, arg)
            else:
                arg = Calculus.convert(arg)
        head, data = arg.pair
        if head is APPLY:
            func, args = data
            if func == cls.get_inverse_function():
                assert len (args)==1,`args`
                return args[0]

        x, m = get_pi_shift(arg, 12)
        m %= (12*cls.period)
        if x == zero:
            v = cls.eval_direct(m)
            if v is not None:
                return v
        period = cls.period
        negate_result = False
        # Full-period symmetry
        if not m % (12*period):
            arg = x
        else:
            # Half-period symmetry
            if not m % 12:
                arg = x
                negate_result = True
            # Quarter-period symmetry
            elif not m % 6:
                f = conjugates[cls]
                sign = (-1)**((((m-6)//12) % period) + (f.parity == 'odd'))
                return sign * f(x)
        if has_leading_sign(arg):
            arg = -arg
            negate_result ^= (cls.parity == 'odd')
        if negate_result:
            return -Apply(cls, arg)
        else:
            return Apply(cls, arg)

    @classmethod
    def get_inverse_function(cls):
        return

class Sin(TrigonometricFunction):
    parity = 'odd'
    period = 2

    @classmethod
    def eval_direct(cls, m):
        return sine_table[m % 24]

    @classmethod
    def derivative(cls, arg, n=1):
        if n == 1:
            return Cos(arg)

    @classmethod
    def nth_derivative(cls, arg, n):
        return Sin(arg + n*pi/2)

    @classmethod
    def fdiff(cls, FAlgebra, argument_index, order):
        if argument_index!=0:
            if isinstance(argument_index, inttypes):
                raise TypeError('%s.fdiff argument_index must be 0 but got %r' % (cls.__name__, argument_index))
            return NotImplemented
        if isinstance(order, inttypes):
            n = order % 4
            if n==0:
                return FAlgebra(CALLABLE, Sin)
            if n==1:
                return FAlgebra(CALLABLE, Cos)
            if n==2:
                return -FAlgebra(CALLABLE, Sin)
            return -FAlgebra(CALLABLE, Cos)
        order = order.data
        func = lambda x: Sin(x + order * pi/2)
        return FAlgebra(CALLABLE, func)

    @classmethod
    def get_inverse_function(cls):
        return CalculusFunctionRing(CALLABLE, ArcSin)

class Cos(TrigonometricFunction):
    parity = 'even'
    period = 2

    @classmethod
    def eval_direct(cls, m):
        return sine_table[(m+6) % 24]

    @classmethod
    def derivative(cls, arg, n=1):
        return -Sin(arg)

    @classmethod
    def nth_derivative(cls, arg, n=1):
        return Cos(arg + n*pi/2)

    @classmethod
    def fdiff(cls, FAlgebra, argument_index, order):
        if argument_index!=0:
            if isinstance(argument_index, inttypes):
                raise TypeError('%s.fdiff argument_index must be 0 but got %r' % (cls.__name__, argument_index))
            return NotImplemented
        if isinstance(order, inttypes):
            n = (order+1) % 4
            if n==0:
                return FAlgebra(CALLABLE, Sin)
            if n==1:
                return FAlgebra(CALLABLE, Cos)
            if n==2:
                return -FAlgebra(CALLABLE, Sin)
            return -FAlgebra(CALLABLE, Cos)
        func = lambda x: Cos(x + order * pi/2)
        return FAlgebra(CALLABLE, func)

class Tan(TrigonometricFunction):
    parity = 'odd'
    period = 1

    @classmethod
    def eval_direct(cls, m):
        a = sine_table[m % 24]
        b = sine_table[(m+6) % 24]
        if a is None or b is None:
            return
        return a / b

    @classmethod
    def derivative(cls, arg):
        return 1+Tan(arg)**2

    @classmethod
    def fdiff(cls, FAlgebra, argument_index, order):
        if argument_index!=0:
            if isinstance(argument_index, inttypes):
                raise TypeError('%s.fdiff argument_index must be 0 but got %r' % (cls.__name__, argument_index))
            return NotImplemented
        if order==1:
            tan = FAlgebra(CALLABLE, cls)
            return tan**2 + 1
        return NotImplemented

class Cot(TrigonometricFunction):
    parity = 'odd'
    period = 1

    @classmethod
    def eval_direct(cls, m):
        a = sine_table[m % 24]
        b = sine_table[(m+6) % 24]
        if a is None or b is None:
            return
        return b / a

    @classmethod
    def derivative(cls, arg, n=1):
        return -(1+Cot(arg)**2)

    @classmethod
    def fdiff(cls, FAlgebra, argument_index, order):
        if argument_index!=0:
            if isinstance(argument_index, inttypes):
                raise TypeError('%s.fdiff argument_index must be 0 but got %r' % (cls.__name__, argument_index))
            return NotImplemented
        if order==1:
            cot = FAlgebra(CALLABLE, cls)
            return -cot**2 - 1
        return NotImplemented

# pi/2-x symmetry
conjugates = {
  Sin : Cos,
  Cos : Sin,
  Tan : Cot,
  Cot : Tan
}

class ArcSin(CalculusDefinedFunction):
    def __new__(cls, x):
        if not isinstance(x, Calculus):
            x = Calculus.convert(x)
        return Apply(cls, x)
