#
# Created January 2008 by Pearu Peterson
#
""" Provides Constant class.
"""
__docformat__ = "restructuredtext"
__all__ = ['Constant', 'const_pi', 'const_E', 'const_gamma']

from ..core import classes
from ..utils import SYMBOL, NUMBER
from ..arithmetic.evalf import evalf
from ..arithmetic import mpmath, setdps

class Constant(str):

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, str(self))

    def evalf(self, n=None):
        if n:
            setdps(n)
        if self == 'pi':
            return +mpmath.pi
        if self == 'E':
            return +mpmath.e
        if self == 'gamma':
            return +mpmath.euler
        raise NotImplementedError('%s(%r).evalf'
                                  % (self.__class__.__name__, self))

    def as_algebra(self, cls):
        if cls is classes.Calculus:
            return cls.Symbol(self)
        return cls.convert(self)

    def get_direction(self):
        if self in ['pi', 'E']:
            return 1
        return NotImplemented

    @property
    def is_bounded(self):
        if self in ['pi', 'E']:
            return True
        return None

const_pi = Constant('pi')
const_E = Constant('E')
const_gamma = Constant('gamma')
