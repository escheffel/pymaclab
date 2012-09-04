""" Implements matrix determinant algoritms.
"""
# Author: Pearu Peterson
# Created: March 2008
from ..arithmetic.numbers import div
from .linalg import swap_rows_MATRIX

def MATRIX_DICT_determinant(self, overwrite=False):
    """ Determinant of a n x n matrix.

    Using fraction-free Gaussian elimination method when n>3.
    For n<=3, direct expressions are used.
    """
    head, data = self.pair

    m, n = head.shape
    assert m==n,`m,n`
    if head.is_diagonal:
        dget = data.get
        d = 1
        for i in range (m):
            v = dget((i,i))
            if v is None:
                return 0
            d *= v
        return d
    if not (overwrite and self.is_writable):
        data = dict(data)

    return determinant_MATRIX(m, data)

def get_minor_MATRIX(p, q, data):
    d = {}
    for (i,j),x in data.items():
        if p==i or q==j:
            continue
        if i<p:
            if j<q:
                d[i,j] = x
            else:
                d[i,j-1] = x
        else:
            if j<q:
                d[i-1,j] = x
            else:
                d[i-1,j-1] = x
    return d

def _det2(a_00, a_01, a_10, a_11):
    if a_00 and a_11:
        if a_01 and a_10:
            return a_00*a_11 - a_01*a_10
        return a_00*a_11
    if a_01 and a_10:
        return -a_01*a_10
    return 0

def _det3(a_00, a_01, a_02, a_10, a_11, a_12, a_20, a_21, a_22):
    if a_00:
        if a_01:
            if a_02:
                return a_00 * _det2(a_11, a_12, a_21, a_22) \
                       - a_01 * _det2(a_10, a_12, a_20, a_22) \
                       + a_02 * _det2(a_10, a_11, a_20, a_21)
            return a_00 * _det2(a_11, a_12, a_21, a_22) \
                   - a_01 * _det2(a_10, a_12, a_20, a_22)
        if a_02:
            return a_00 * _det2(a_11, a_12, a_21, a_22) \
                   + a_02 * _det2(a_10, a_11, a_20, a_21)
        return a_00 * _det2(a_11, a_12, a_21, a_22)
    if a_01:
        if a_02:
            return - a_01 * _det2(a_10, a_12, a_20, a_22) \
                   + a_02 * _det2(a_10, a_11, a_20, a_21)
        return - a_01 * _det2(a_10, a_12, a_20, a_22)
    if a_02:
        return + a_02 * _det2(a_10, a_11, a_20, a_21)
    return 0

def determinant_MATRIX(m, data):
    data_get = data.get
    if m==1:
        return data_get((0,0),0)
    if m==2:
        a_00 = data_get((0,0))
        a_01 = data_get((0,1))
        a_10 = data_get((1,0))
        a_11 = data_get((1,1))
        return _det2(a_00, a_01, a_10, a_11)
    if m==3:
        a_00 = data_get((0,0))
        a_01 = data_get((0,1))
        a_02 = data_get((0,2))
        a_10 = data_get((1,0))
        a_11 = data_get((1,1))
        a_12 = data_get((1,2))
        a_20 = data_get((2,0))
        a_21 = data_get((2,1))
        a_22 = data_get((2,2))
        return _det3(a_00, a_01, a_02, a_10, a_11, a_12, a_20, a_21, a_22)
    
    newdata = {}
    d = 1
    sign = 1
    for k in range(m-1):
        a_kk = data_get((k,k))
        k1 = k+1
        if not a_kk:
            for i in range(k1,m):
                a_kk = data_get((i,k))
                if a_kk:
                    break
            if not a_kk:
                return 0
            sign = -sign
            swap_rows_MATRIX(data, i, k)
        for i in range(k1,m):
            a_ik = data_get((i,k))
            if a_ik:
                for j in range(k1, m):
                    ij = i,j
                    a_ij = data_get(ij)
                    if a_ij:
                        a_kj = data_get((k,j))
                        if a_kj:
                            x = a_kk * a_ij - a_ik * a_kj
                            if x:
                                if d is 1:
                                    newdata[ij] = x
                                else:
                                    newdata[ij] = div(x, d)
                        else:
                            if d is 1:
                                newdata[ij] = a_kk * a_ij
                            else:
                                newdata[ij] = div(a_kk * a_ij, d)
                    else:
                        a_kj = data_get((k,j))
                        if a_kj:
                            if d is 1:
                                newdata[ij] = -a_ik * a_kj
                            else:
                                newdata[ij] = div(-a_ik * a_kj, d)
            else:
                for j in range(k1, m):
                    ij = i,j
                    a_ij = data_get(ij)
                    if a_ij:
                        if d is 1:
                            newdata[ij] = a_kk * a_ij
                        else:
                            newdata[ij] = div(a_kk * a_ij, d)
        d = a_kk
        data = newdata
        data_get = data.get
        newdata = {}
    return sign * data.get((m-1,m-1),0)
