""" Implementes functions ring support.
"""
#
# Author: Pearu Peterson
# Created: April, 2008
#

__all__ = ['FunctionRing']

from ..core import classes, objects, init_module
from ..basealgebra import Verbatim, Algebra
from ..ring import CommutativeRing

init_module.import_heads()


class FunctionRing(CommutativeRing):
    """ Base class to functions ring classes.

    Use ``Function`` function to construct instances.
    """

    argument_algebras = None
    nargs = None

    @classmethod
    def get_value_algebra(cls):
        return CommutativeRing

    def get_argument_algebra(self, index):
        return self.get_value_algebra()

    @classmethod
    def get_function_algebra(cls):
        return classes.OperatorRing

    @classmethod
    def get_differential_algebra(cls):
        return classes.DifferentialRing

    @classmethod
    def get_predefined_symbols(cls, name):
        if name=='D': return D
        return

    @classmethod
    def convert(cls, obj, typeerror=True):
        tobj = type(obj)
        if tobj is cls:
            return obj        
        if isinstance(obj, cls.get_value_algebra()):
            return cls(NUMBER, obj)
        return super(CommutativeRing, cls).convert(obj, typeerror=typeerror)

    def as_algebra(self, cls, typeerror=True):
        if cls is classes.Verbatim:
            return self.as_verbatim()
        if type(self) is cls:
            return self
        #if isinstance(self, cls):
        #    return self.as_verbatim().as_algebra(cls)
        if typeerror:
            raise TypeError('Cannot convert %s to %s instance' % (type(self).__name__, cls.__name__))
        return NotImplemented

    def __call__(self, *args, **options):
        cls = self.get_value_algebra()
        #cls = classes.Calculus
        evaluate = options.get('evaluate', True)
        if evaluate:
            result = self.head.apply(cls, self.data, self, args)
            if result is not NotImplemented:
                return result
        return cls(APPLY, (self, args))

classes.FunctionRing = FunctionRing

