"""Provides various implementation specific constants (expression heads, etc).
"""

__docformat__ = 'restructuredtext'

import inspect

class HEAD(object):
    """ Base class to head constants.

    Head constants are singletons.
    """
    _cache = {}
    def __new__(cls, *args):
        key = '%s%s' % (cls.__name__, args)
        obj = cls._cache.get(key)
        if obj is None:
            obj = object.__new__(cls)
            cls._cache[key] = obj
            obj._key = key
            obj.init(*args)
        return obj

    def init(self, *args):
        # derived class may set attributes here
        pass #pragma NO COVER

    def as_unique_head(self):
        # used by the pickler support to make HEAD instances unique
        try:
            return self._cache[self._key]
        except KeyError:
            self._cache[self._key] = self
            return self

    def is_data_ok(self, cls, expr):
        return

# The following constants define both the order of operands
# as well as placing parenthesis for classes deriving from
# CollectingField:

str_SUM = -1
str_PRODUCT = -2
str_POWER = -3
str_APPLY = -4
str_SYMBOL = -5
str_NUMBER = -6

# The following constants are used by Verbatim and
# CollectingField classes.



#POLY = HEAD('POLY')
#DENSE_POLY = HEAD('DENSE_POLY')

DIFF = HEAD('DIFF')

MATRIX_DICT = intern('MATRIX_DICT')
MATRIX_DICT_T = intern('MATRIX_DICT_T')
MATRIX_DICT_A = intern('MATRIX_DICT_A')
MATRIX_DICT_TA = intern('MATRIX_DICT_TA')
MATRIX_DICT_D = intern('MATRIX_DICT_D')
MATRIX_DICT_TD = intern('MATRIX_DICT_TD')

class MATRIX(HEAD):
    """ Matrix head singleton class.

    Usage::

      MATRIX(<rows>, <cols>, <storage>)

    where
    
      ``<rows>``     - number of matrix rows
      ``<cols>``     - number of matrix columns
      ``<strorage>`` - constant describing data storage properties:
                       MATRIX_DICT, MATRIX_DICT_T, MATRIX_DICT_A, MATRIX_DICT_TA,
                       MATRIX_DICT_D, MATRIX_DICT_DT
    """

    def __str__ (self):
        return '%s(%r, %r, %r)' % (self.__class__.__name__, self.rows, self.cols, self.storage)

    __repr__ = __str__

    def to_lowlevel(self, cls, data, pair):
        return pair
    
    def init(self, rows, cols, storage):
        self.rows = rows
        self.cols = cols
        self.shape = (rows, cols)
        self.storage = storage

        self.is_transpose = is_transpose = storage in [MATRIX_DICT_T, MATRIX_DICT_TA, MATRIX_DICT_TD]
        self.is_array = storage in [MATRIX_DICT_A, MATRIX_DICT_TA]
        self.is_diagonal = storage in [MATRIX_DICT_D, MATRIX_DICT_TD]

        if storage==MATRIX_DICT:
            self.T = type(self)(cols, rows,  MATRIX_DICT_T)
            self.A = type(self)(rows, cols,  MATRIX_DICT_A)
            self.M = self
            self.D = type(self)(rows, cols,  MATRIX_DICT_D)
        elif storage==MATRIX_DICT_T:
            self.T = type(self)(cols, rows,  MATRIX_DICT)
            self.A = type(self)(rows, cols,  MATRIX_DICT_TA)
            self.M = self
            self.D = type(self)(rows, cols,  MATRIX_DICT_TD)
        elif storage==MATRIX_DICT_A:
            self.T = type(self)(cols, rows,  MATRIX_DICT_TA)
            self.A = self
            self.M = type(self)(rows, cols,  MATRIX_DICT)
            self.D = type(self)(rows, cols,  MATRIX_DICT_D)
        elif storage==MATRIX_DICT_TA:
            self.T = type(self)(cols, rows,  MATRIX_DICT_A)
            self.A = self
            self.M = type(self)(rows, cols,  MATRIX_DICT_T)
            self.D = type(self)(rows, cols,  MATRIX_DICT_TD)
        elif storage==MATRIX_DICT_D:
            self.T = type(self)(cols, rows,  MATRIX_DICT_T)
            self.A = type(self)(rows, cols,  MATRIX_DICT_A)
            self.M = type(self)(rows, cols,  MATRIX_DICT)
            self.D = self
        elif storage==MATRIX_DICT_TD:
            self.T = type(self)(cols, rows,  MATRIX_DICT)
            self.A = type(self)(rows, cols,  MATRIX_DICT_TA)
            self.M = type(self)(rows, cols,  MATRIX_DICT_T)
            self.D = self
        else:
            raise NotImplementedError(`storage`) #pragma NO COVER

def totree(obj, tab=''):
    from .core import Expr
    if isinstance(obj, Expr):
        head, data = obj.pair
        s = '%s%s head:%r\n' % (tab, type(obj).__name__, head)
        s += '%s' % (totree(data, tab+'  '))
        return s
    elif isinstance(obj, dict):
        return '%s%s:\n%s' % (tab, type(obj).__name__, '\n'.join([totree(item, tab+'  ') for item in obj.iteritems()]))
    elif isinstance(obj, (list, tuple)):
        return '%s%s:\n%s' % (tab, type(obj).__name__, '\n'.join([totree(item, tab+'  ') for item in obj]))
    else:
        return '%s%s:%r' % (tab, type(obj).__name__, obj)


def test_operations(operands, expected_results, unary_operations, binary_operations):
    from sympycore import Expr

    results = {}
    for line in expected_results.split('\n'):
        line = line.strip()
        if ':' not in line: continue
        expr, result = line.split(':')
        for e in expr.split(';'):
            e = e.strip()
            results[e] = [r.strip() for r in result.split(';')]

    if isinstance(operands, tuple):
        operands1, operands2 = operands
    else:
        operands1 = operands2 = operands

    for op1 in operands1:
        if isinstance(op1, Expr):
            for op in unary_operations:
                expr = '%s(%s)' % (op, op1) 

                try:
                    result_obj = eval('%s(op1)' % (op), dict(op1=op1))
                    result = str(result_obj)
                except Exception, msg:
                    print  expr,'failed with %s' % (msg)
                    raise
                
                if expr not in results:
                    print '%s:%s' % (expr, result)
                    continue
                try:
                    assert result in results[expr], `results[expr], result`
                except AssertionError:
                    print 
                    print 'op1: %s' % (op1)
                    print totree(op1, '  ')
                    print '%r result: %s' % (op, result_obj)
                    print totree(result_obj, '  ')
                    raise
        for op2 in operands2:
            if not (isinstance(op1, Expr) or isinstance(op2, Expr)):
                continue
            for op in binary_operations:
                expr = '(%s)%s(%s)' % (op1, op, op2)

                result_obj = None
                try:
                    result_obj = eval('(op1)%s(op2)' % op, dict(op1=op1, op2=op2))
                    result = str(result_obj)
                except Exception, msg:
                    result = str(msg)
                    if not result.startswith('unsupported'):
                        print  expr,`op1,op,op2`,'failed with %s' % (msg)
                        raise
                if expr not in results:
                    print '%s:%s' % (expr, result)
                    continue
                if result.startswith('unsupported'):
                    assert 'unsupported' in results[expr], `results[expr], result, op1, op2, expr`
                else:
                    try:
                        assert result in results[expr], `results[expr], result`
                    except AssertionError:
                        print
                        print 'op1: %s' % (op1)
                        print totree(op1, '  ')
                        print 'op2: %s' % (op2)
                        print totree(op2, '  ')
                        print '%r result: %s' % (op, result_obj)
                        print totree(result_obj, '  ')
                        raise


from heads import *
