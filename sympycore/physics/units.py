""" Provides Unit class.
"""

__docformat__ = "restructuredtext"
__all__ = [\
    'Unit', 'meter', 'kilogram', 'second', 'ten', 'yotta', 'zetta',
    'exa', 'peta', 'tera', 'giga', 'mega', 'kilo', 'deca', 'deci',
    'centi', 'milli', 'micro', 'nano', 'pico', 'femto', 'atto',
    'zepto', 'yocto' ]

from ..core import classes

from ..calculus.algebra import Calculus, algebra_numbers
from ..ring import CommutativeRing

class Unit(CommutativeRing):
    """ Represents an algebra of physical units.

    Elements of the units algebra are unit symbols.
    Coefficients of the units algebra are Calculus instances.
    """

    def as_algebra(self, cls, typeerror=True):
        """ Convert algebra to another algebra.
        """
        if cls is classes.Verbatim:
            return self.as_verbatim()
        if cls is Calculus:
            return NotImplemented
        return self.as_verbatim().as_algebra(cls)

classes.Unit = Unit

meter = Unit.Symbol('m')
kilogram = Unit.Symbol('kg')
second = Unit.Symbol('s')

def init_module(m):


    m.ten = ten = Calculus.Number(10)
    m.yotta = ten**24
    m.zetta = ten**21
    m.exa   = ten**18
    m.peta  = ten**15
    m.tera  = ten**12
    m.giga  = ten**9
    m.mega  = ten**6
    m.kilo  = ten**3
    m.deca  = ten**1
    m.deci  = ten**-1
    m.centi = ten**-2
    m.milli = ten**-3
    m.micro = ten**-6
    m.nano  = ten**-9
    m.pico  = ten**-12
    m.femto = ten**-15
    m.atto  = ten**-18
    m.zepto = ten**-21
    m.yocto = ten**-24

