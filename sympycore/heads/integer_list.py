
__all__ = ['INTEGER_LIST']

from .base import Head

class IntegerListHead(Head):

    def __repr__(self): return 'INTEGER_LIST'

    def is_data_ok(self, cls, data):
        if type(data) is not list:
            return 'data must be list but got %r' % (type(data).__name__) #pragma: no cover
        for i, item in enumerate(data):
            if not isinstance(item, (int, long)): #pragma: no cover
                return 'all data items must be integers but %s-th item is %r' % (i, type(item).__name__) #pragma: no cover

    def to_lowlevel(self, cls, data, pair):
        if len(data)==1:
            return data[0]
        return pair

INTEGER_LIST = IntegerListHead()
