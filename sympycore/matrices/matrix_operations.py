
from ..utils import MATRIX, MATRIX_DICT
from .algebra import Matrix, MatrixBase, MatrixDict

from ..core import init_module
init_module.import_lowlevel_operations()

def MATRIX_DICT_iadd(self, other):
    """ Inplace matrix add.
    """
    t = type(other)
    if t is list or t is tuple:
        other = Matrix(other)
        t = type(other)
    if t is MatrixDict:
        if self.is_writable:
            ret = self
        else:
            ret = self.copy()
        head1, data1 = ret.pair
        head2, data2 = other.pair
        assert head1.shape==head2.shape,`head1, head2`
        if head1.is_transpose:
            if head2.is_transpose:
                iadd_MATRIX_MATRIX_TT(data1, data2)
            else:
                iadd_MATRIX_MATRIX_TA(data1, data2)
        elif head2.is_transpose:
            iadd_MATRIX_MATRIX_AT(data1, data2)
        else:
            iadd_MATRIX_MATRIX_AA(data1, data2)
        if head2.is_array:
            return ret.A
        return ret.M
    elif isinstance(other, MatrixBase):
        raise NotImplementedError(`type(other)`)
    else:
        if other:
            if self.is_writable:
                ret = self
            else:
                ret = self.copy()
            head, data = ret.pair
            rows, cols = head.shape
            if head.is_transpose:
                iadd_MATRIX_T_SCALAR(rows, cols, data, other)
            else:
                iadd_MATRIX_SCALAR(rows, cols, data, other)
            return ret
        else:
            return self

def MATRIX_DICT_imul(self, other):
    """ Inplace matrix multiplication.
    """
    t = type(other)
    if t is list or t is tuple:
        other = Matrix(other)
        t = type(other)
    if t is MatrixDict:
        head1, data1 = self.pair
        head2, data2 = other.pair
        if head1.is_array or head2.is_array:
            assert head1.shape==head2.shape,`head1, head2`
            if self.is_writable:
                ret = self
            else:
                ret = self.copy()
            data1 = ret.data
            if head1.is_transpose:
                if head2.is_transpose:
                    imul_MATRIX_MATRIX_ATT(data1, data2)
                else:
                    imul_MATRIX_MATRIX_TA(data1, data2)
            elif head2.is_transpose:
                imul_MATRIX_MATRIX_AT(data1, data2)
            else:
                imul_MATRIX_MATRIX_AA(data1, data2)
            if head2.is_array:
                return ret.A
            return ret.M
        else:
            assert head1.cols==head2.rows,`head1, head2`
            args = data1, data2, head1.rows, head2.cols, head1.cols
            if head1.is_transpose:
                if head2.is_transpose:
                    ret = mul_MATRIX_MATRIX_MTT(*args)
                else:
                    ret = mul_MATRIX_MATRIX_TM(*args)
            elif head2.is_transpose:
                ret = mul_MATRIX_MATRIX_MT(*args)
            else:
                ret = mul_MATRIX_MATRIX_MM(*args)
            return ret
    elif isinstance(other, MatrixBase):
        raise NotImplementedError(`type(other)`)
    else:
        if self.is_writable:
            ret = self
        else:
            ret = self.copy()
        if other:
            head, data = ret.pair
            for key in data:
                data[key] *= other
        else:
            head, data = ret.pair
            data.clear()
        return ret

def iadd_MATRIX_SCALAR(rows, cols, data, value):
    col_indices = range(cols)
    for i in xrange(rows):
        for j in col_indices:
            key = (i,j)
            b = data.get(key)
            if b is None:
                data[key] = value
            else:
                b += value
                if b:
                    data[key] = b
                else:
                    del data[key]

def iadd_MATRIX_T_SCALAR(rows, cols, data, value):
    col_indices = range(cols)
    for i in xrange(rows):
        for j in col_indices:
            key = (j,i)
            b = data.get(key)
            if b is None:
                data[key] = value
            else:
                b += value
                if b:
                    data[key] = b
                else:
                    del data[key]

def iadd_MATRIX_MATRIX_AA(data1, data2):
    for key,x in data2.items():
        b = data1.get(key)
        if b is None:
            data1[key] = x
        else:
            b += x
            if b:
                data1[key] = b
            else:
                del data1[key]

def iadd_MATRIX_MATRIX_AT(data1, data2):
    for (j,i),x in data2.items():
        key = i,j
        b = data1.get(key)
        if b is None:
            data1[key] = x
        else:
            b += x
            if b:
                data1[key] = b
            else:
                del data1[key]

def iadd_MATRIX_MATRIX_TA(data1, data2):
    for (i,j),x in data2.items():
        key = j,i
        b = data1.get(key)
        if b is None:
            data1[key] = x
        else:
            b += x
            if b:
                data1[key] = b
            else:
                del data1[key]

iadd_MATRIX_MATRIX_TT = iadd_MATRIX_MATRIX_AA

def imul_MATRIX_MATRIX_AA(data1, data2):
    for key in data1:
        b = data2.get(key)
        if b is None:
            del data1[key]
        else:
            data1[key] *= b

def imul_MATRIX_MATRIX_AT(data1, data2):
    for key in data1:
        i,j = key
        b = data2.get((j,i))
        if b is None:
            del data1[key]
        else:
            data1[key] *= b

imul_MATRIX_MATRIX_TA = imul_MATRIX_MATRIX_AT
imul_MATRIX_MATRIX_ATT = imul_MATRIX_MATRIX_AA

def mul_MATRIX_MATRIX_AA(data1, data2, rows, cols):
    indices = range(n)
    d = {}
    data1_get = data1.get
    data2_get = data2.get
    for i in xrange(rows):
        for j in xrange(cols):
            key = i,j
            a_ij = data1_get(key)
            if a_ij is None:
                continue
            b_ij = data2_get(key)
            if b_ij is None:
                continue
            d[key] = a_ij * b_ij
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

def mul_MATRIX_MATRIX_AT(data1, data2, rows, cols):
    indices = range(n)
    d = {}
    data1_get = data1.get
    data2_get = data2.get
    for i in xrange(rows):
        for j in xrange(cols):
            key = i,j
            a_ij = data1_get(key)
            if a_ij is None:
                continue
            b_ij = data2_get((j,i))
            if b_ij is None:
                continue
            d[key] = a_ij * b_ij
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

def mul_MATRIX_MATRIX_TA(data1, data2, rows, cols):
    indices = range(n)
    d = {}
    data1_get = data1.get
    data2_get = data2.get
    for i in xrange(rows):
        for j in xrange(cols):
            key = i,j
            a_ij = data1_get((j,i))
            if a_ij is None:
                continue
            b_ij = data2_get(key)
            if b_ij is None:
                continue
            d[key] = a_ij * b_ij
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

def mul_MATRIX_MATRIX_ATT(data1, data2, rows, cols):
    indices = range(n)
    d = {}
    data1_get = data1.get
    data2_get = data2.get
    for i in xrange(rows):
        for j in xrange(cols):
            ikey = j,i
            a_ij = data1_get(ikey)
            if a_ij is None:
                continue
            b_ij = data2_get(ikey)
            if b_ij is None:
                continue
            d[i,j] = a_ij * b_ij
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

def mul_MATRIX_MATRIX_MM(data1, data2, rows, cols, n):
    d = {}
    left_indices = {}
    right_indices = {}
    for i,j in data1.keys():
        s = left_indices.get(i)
        if s is None:
            s = left_indices[i] = set()
        s.add(j)
    for j,k in data2.keys():
        s = right_indices.get(k)
        if s is None:
            s = right_indices[k] = set()
        s.add(j)

    for i,ji in left_indices.items ():
        for k,jk in right_indices.items ():
            for j in ji.intersection(jk):
                dict_add_item(None, d, (i,k), data1[(i,j)] * data2[(j,k)])
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

def mul_MATRIX_MATRIX_TM(data1, data2, rows, cols, n):
    d = {}
    left_indices = {}
    right_indices = {}
    for j,i in data1.keys():
        s = left_indices.get(i)
        if s is None:
            s = left_indices[i] = set()
        s.add(j)
    for j,k in data2.keys():
        s = right_indices.get(k)
        if s is None:
            s = right_indices[k] = set()
        s.add(j)

    for i,ji in left_indices.items ():
        for k,jk in right_indices.items ():
            for j in ji.intersection(jk):
                dict_add_item(None, d, (i,k), data1[(j,i)] * data2[(j,k)])
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

def mul_MATRIX_MATRIX_MT(data1, data2, rows, cols, n):
    d = {}
    left_indices = {}
    right_indices = {}
    for i,j in data1.keys():
        s = left_indices.get(i)
        if s is None:
            s = left_indices[i] = set()
        s.add(j)
    for k,j in data2.keys():
        s = right_indices.get(k)
        if s is None:
            s = right_indices[k] = set()
        s.add(j)

    for i,ji in left_indices.items ():
        for k,jk in right_indices.items ():
            for j in ji.intersection(jk):
                dict_add_item(None, d, (i,k), data1[(i,j)] * data2[(k,j)])
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

def mul_MATRIX_MATRIX_MTT(data1, data2, rows, cols, n):
    d = {}
    left_indices = {}
    right_indices = {}
    for j,i in data1.keys():
        s = left_indices.get(i)
        if s is None:
            s = left_indices[i] = set()
        s.add(j)
    for k,j in data2.keys():
        s = right_indices.get(k)
        if s is None:
            s = right_indices[k] = set()
        s.add(j)

    for i,ji in left_indices.items ():
        for k,jk in right_indices.items ():
            for j in ji.intersection(jk):
                dict_add_item(None, d, (i,k), data1[(j,i)] * data2[(k,j)])
    return MatrixDict(MATRIX(rows, cols, MATRIX_DICT), d)

