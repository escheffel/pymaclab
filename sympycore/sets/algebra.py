""" Implements Set algebra support.
"""
#
# Author: Pearu Peterson
# Created: April, 2008
#

__all__ = ['Set', 'Integers']

from ..utils import NUMBER, SYMBOL, IN
from ..core import classes
from ..basealgebra import Algebra
from ..arithmetic.numbers import inttypes_set, realtypes_set, complextypes_set, numbertypes_set


class Set(Algebra):
    """ Set algebra.

    Set algebra 'numbers' are Python sets, for example, ``Set(NUMBER,
    frozenset([1,2]))`` represents a set ``{1, 2}``.

    The following Set algebra symbols have special meanings::

      Integers = Set(SYMBOL, 'Integers')
      Reals = Set(SYMBOL, 'Reals')
      ...
    
    """

    def get_element_algebra(self):
        head, data = self.pair
        if head is SYMBOL:
            if data in ['Integers', 'Reals']:
                return classes.Calculus
        return classes.Verbatim

    @classmethod
    def convert_number(cls, obj, typeerror=True):
        t = type(obj)
        if t in [set, frozenset, list, tuple]:
            return frozenset(obj)
        return cls.handle_convert_failure('number', obj, typeerror)

    def contains(self, obj):
        head, data = self.pair
        if head is NUMBER:
            if obj in data:
                return classes.Logic.true
            return classes.Logic.false
        if head is SYMBOL:
            t = type(obj)
            if data=='Integers':
                if t in inttypes_set:
                    return classes.Logic.true
                elif t in numbertypes_set:
                    return classes.Logic.false
                elif isinstance(obj, Algebra):
                    h, d = obj.pair
                    if h is NUMBER:
                        return self.contains(d)
        return classes.Logic(IN, (obj, self))

    def __str__(self):
        head, data = self.pair
        if head is NUMBER:
            return '{%s}' % (', '.join(map(str, data)))
        return Algebra.__str__(self)

classes.Set = Set
Integers = Set(SYMBOL, 'Integers')
