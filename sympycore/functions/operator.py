
__all__ = ['OperatorRing']

from ..core import classes, init_module
from .algebra import FunctionRing

init_module.import_heads()
init_module.import_numbers()

class OperatorRing(FunctionRing):

    @classmethod
    def convert(cls, obj, typeerror=True):
        if isinstance(obj, cls):
            return obj
        fcls = cls.get_value_algebra()
        if isinstance(obj, fcls):
            return cls(NUMBER, obj)
        if isinstance(obj, fcls.get_value_algebra()):
            return cls(NUMBER, obj)
        if isinstance(obj, (str, unicode, classes.Verbatim)):
            return classes.Verbatim.convert(obj).as_algebra(cls)
        return cls.handle_convert_failure('algebra', obj, typeerror)

    def as_algebra(self, cls, typeerror=True):
        raise NotImplementedError(`cls`)

    @classmethod
    def get_value_algebra(cls):
        return classes.FunctionRing

    @classmethod
    def get_function_algebra(cls):
        raise NotImplementedError(`cls`)

    @classmethod
    def get_differential_algebra(cls):
        raise NotImplementedError(`cls`)

classes.OperatorRing = OperatorRing
