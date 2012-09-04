""" Provides matrix utility functions.
"""
#
# Author: Pearu Peterson
# Created: March 2008
#

__docformat__ = "restructuredtext"
__all__ = ['eye', 'concatenate', 'jacobian']

from ..utils import MATRIX, MATRIX_DICT
from .algebra import MatrixDict, Matrix, MatrixBase

def jacobian(expr_list, var_list):
    """ Return a jacobian matrix of functions in expr_list with
    respect to variables in var_list.
    """
    m, n = len(expr_list), len(var_list)
    jac = Matrix(m, n)
    for i, e in enumerate(expr_list):
        for j, v in enumerate(var_list):
            jac[i,j] = e.diff(v)
    return jac

def eye(m, n=None, k=0):
    """ Return n x m matrix where the k-th diagonal is all ones,
    everything else is zeros.
    """
    if n is None:
        n = m
    d = {}
    if k<0:
        for i in range(min(m,n)+k):
            d[i-k,i] = 1
    else:
        for i in range(min(m,n)-k):
            d[i,i+k] = 1
    return MatrixDict(MATRIX(m, n, MATRIX_DICT), d)

def concatenate(*args, **kws):
    """ Join matrices together.

    concatenate(m1, m2, ..) - join matrices m1, m2, .. together along rows
    concatenate(m1, m2, .., axis=1) - join matrices m1, m2, .. together along columns
    concatenate(m1, m2, .., diagonal=True) - join matrices m1, m2, .. together along diagonal
    """
    assert args,`args`
    assert set(kws).issubset(['axis','diagonal']) and len(kws)<=1,`kws`
    axis = kws.get('axis',0)
    diagonal = kws.get('diagonal',False)
    d = {}
    rows, cols = 0, 0
    if axis:
        for a in args:
            head, data = Matrix(a)[:].pair
            assert not (head.is_transpose or head.is_diagonal)
            m, n = head.shape
            for (i,j),x in data.items():
                d[rows + i,j] = x
            rows += m
            cols = max(cols, n)
    elif diagonal:
        for a in args:
            head, data = Matrix(a)[:].pair
            assert not (head.is_transpose or head.is_diagonal)
            m, n = head.shape
            for (i,j),x in data.items():
                d[rows + i, cols + j] = x
            rows += m
            cols += n
    else:
        for a in args:
            head, data = Matrix(a)[:].pair
            assert not (head.is_transpose or head.is_diagonal)
            m, n = head.shape
            for (i,j),x in data.items():
                d[i,cols+j] = x
            cols += n
            rows = max(rows, m)
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)
