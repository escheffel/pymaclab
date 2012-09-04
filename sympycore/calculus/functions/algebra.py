
from ...core import init_module, classes
from ...functions import FunctionRing, DifferentialRing, OperatorRing

from ..algebra import Calculus

init_module.import_heads()

class CalculusFunctionRing(FunctionRing):

    """
    Functions ring with Calculus value algebra.
    """

    @classmethod
    def get_value_algebra(cls):
        return Calculus

    @classmethod
    def get_function_algebra(cls):
        return CalculusOperatorRing

    @classmethod
    def get_differential_algebra(cls):
        return CalculusDifferentialRing

class CalculusOperatorRing(OperatorRing):

    @classmethod
    def get_value_algebra(cls):
        return CalculusFunctionRing

    @classmethod
    def get_function_algebra(cls):
        raise NotImplementedError(`cls`)

    @classmethod
    def get_differential_algebra(cls):
        raise NotImplementedError(`cls`)
    
class CalculusDifferentialRing(CalculusOperatorRing, DifferentialRing):

    pass

classes.CalculusFunctionRing = CalculusFunctionRing
classes.CalculusOperatorRing = CalculusOperatorRing
classes.CalculusDifferentialRing = CalculusDifferentialRing
