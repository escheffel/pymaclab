#
# Author: Pearu Peterson
# Created: February 2008
#
""" Provides Matrix support.
"""
from __future__ import division

__docformat__ = "restructuredtext"
__all__ = ['Matrix', 'MatrixBase', 'MatrixDict']

import random

from ..core import classes
from ..utils import SYMBOL, NUMBER, ADD, MUL, POW, MATRIX, MATRIX_DICT
from ..basealgebra import Algebra
from ..basealgebra.verbatim import Verbatim
from ..arithmetic.numbers import div, mpq

def is_sequence(obj):
    t = type(obj)
    if t is list or t is tuple:
        return True
    elif t is mpq:
        return False
    try:
        len(obj)
        return True
    except TypeError:
        return False

def is_integer(obj):
    t = type(obj)
    if t is int or t is long:
        return True
    try:
        int(obj)
        return True
    except TypeError:
        pass
    return False

def Matrix(*args, **kws):
    """ Construct a matrix instance.

    The following signatures are supported:

      Matrix(n)           - n x n matrix of zeros
      Matrix(n, random=True) - n x n matrix of random values belonging to [-10,10]
      Matrix(n, random=[l, u]) - n x n matrix of random values belonging to [l,u]
      Matrix(m, n)        - m x n matrix of zeros
      Matrix(n, m, random=True) - m x n matrix of random values belonging to [-10,10]
      Matrix(n, m, random=[l, u]) - m x n matrix of random values belonging to [l,u]
      Matrix(seq)         - n x 1 column matrix where n = len(seq)
      Matrix(seq, diagonal=True)    - n x n diagonal matrix where n = len(seq)
      Matrix(seq, permutation=True) - n x n permutation matrix where n = len(seq)
      Matrix(array)       - m x n matrix where m=len(array), n=max(len(array[:]))
      Matrix(m, n, dict)  - m x n matrix with nonzero elements defined by dict,
                            dict keys must be 2-tuples of integers.
      Matrix(m, n, seq)   - m x n matrix with elements in seq row-wise.
    """
    if len(args)==1:
        a = args[0]
        if isinstance(a, MatrixBase):
            return a
        elif isinstance(a, int):
            m, n = a, a
            data = {}
            interval = kws.get('random')
            if interval is not None:
                if interval==True:
                    interval = (-10, 10)
                for i in range(m):
                    for j in range(n):
                        num = random.randint(*interval)
                        if num:
                            data[i,j] = num                
        elif is_sequence(a):
            m = len(a)
            data = {}
            if True in map(is_sequence, a):
                n = 0
                for i,row in enumerate(a):
                    if is_sequence(row):
                        l = len(row)
                        n = max(n, l)
                        for j, x in enumerate(row):
                            if x:
                                data[i,j] = x
                    else:
                        n = max(n, 1)
                        if row:
                            data[i,0] = row
            else:
                n = m
                if kws.get('diagonal'):
                    for i, x in enumerate(a):
                        if x:
                            data[i, i] = x
                elif kws.get('permutation'):
                    for i,x in enumerate(a):
                        assert 0<=x<m,`x,m`
                        data[i, x] = 1
                else: # single column matrix
                    n = 1
                    for i,x in enumerate(a):
                        if x:
                            data[i, 0] = x
        else:
            raise TypeError('Matrix call with 1 argument: argument must be integer or a sequence, got %r' % (type(a)))
    elif len(args)==2:
        m, n = args
        data = {}
        if not (isinstance(m, int) and isinstance(n, int)):
            raise TypeError('Matrix call with 2 arguments: arguments must be integers, got %r, %r' % (type(m), type(n)))
        interval = kws.get('random')
        if interval is not None:
            if interval==True:
                interval = (-10, 10)
            for i in range(m):
                for j in range(n):
                    num = random.randint(*interval)
                    if num:
                        data[i,j] = num                
    elif len(args)==3:
        m, n, data = args
        flag = True
        if isinstance(data, dict):
            pass
        elif is_sequence(data):
            k = 0
            d = {}
            if data and is_sequence(data[k]):
                for i in range(m):
                    for j in xrange(n):
                        x = data[i][j]
                        if x:
                            d[i,j] = x
            else:
                for i in range(m):
                    for j in xrange(n):
                        x = data[k]
                        if x:
                            d[i,j] = x
                        k += 1
            data = d
        else:
            flag = False
        if not (isinstance(m, int) and isinstance(n, int) and flag):
            raise TypeError('Matrix call with 3 arguments: arguments must be integers and a dictionary|sequence, got %r, %r, %r'\
                            % (type(m), type(n), type(data)))
        
        
    else:
        raise TypeError('Matrix takes 1, 2, or 3 arguments, got %r' % (len(args)))

    return MatrixDict(MATRIX(m, n, MATRIX_DICT), data)


class MatrixBase(Algebra):
    """ Base class to matrix classes.

    Currently the following matrix classes are implemented:

      MatrixDict - matrix content is saved in a dictonary
    """

    @classmethod
    def convert(cls, data, typeerror=True):
        if isinstance(data, list):
            return Matrix(data)
        return super(CommutativeRing, cls).convert(data, typeerror=typeerror)

    @classmethod
    def random(cls, n, m, interval=(-10,10)):
        return Matrix(n, m, random=interval)

    def _str_column_fill(self, columns, j, rowiter1, sep = None, rowiter2 = None):
        col = []
        if rowiter1 is not None:
            for i in rowiter1:
                s = str(self[i,j])
                col.append(s)
        if sep is not None:
            col.append (sep)
        if rowiter2 is not None:
            for i in rowiter2:
                s = str(self[i,j])
                col.append(s)
        width = max(map(len,col)+[0])
        fmt = ' %'+str(width)+'s '
        col = [fmt % (s) for s in col]
        columns.append(col)

    def __str__(self, max_nrows=20, max_ncols=20):
        """ Return pretty string representation of a matrix.
        
        Parameters
        ----------
        max_nrows, max_ncols : int
          Specify the maximum number of rows and columns, respectively, that
          are included in the result.

        """
        rows, cols = self.rows, self.cols
        if self.head.is_diagonal:
            # XXX: what should be the output of diagonal view?
            self = self.M
        columns = []
        if cols>max_ncols and rows>max_nrows:
            for j in range((max_ncols+1)//2):
                self._str_column_fill (columns, j, xrange((max_nrows+1)//2), ':', xrange(rows-max_nrows//2, rows))
            columns.append([' .. ']*rows)
            for j in range(cols-(max_ncols)//2,cols):
                self._str_column_fill (columns, j, xrange((max_nrows+1)//2), ':', xrange(rows-max_nrows//2, rows))
        elif cols>max_ncols:
            for j in range((max_ncols+1)//2):
                self._str_column_fill(columns, j, xrange(rows))
            columns.append([' .. ']*rows)
            for j in range(cols-(max_ncols)//2,cols):
                self._str_column_fill (columns, j, xrange(rows))
        elif rows>max_nrows:
            for j in xrange(cols):
                self._str_column_fill(columns, j, xrange((max_nrows+1)//2), ':', xrange(rows-max_nrows//2, rows))
        else:
            for j in xrange(cols):
                self._str_column_fill (columns, j, xrange(rows), None, None)
        return '\n'.join([''.join(row).rstrip() for row in zip(*columns)])

    def __repr__ (self):
        rows, cols = self.rows, self.cols
        return 'Matrix (%s, %s, %r)' % (rows, cols, self.tolist())

    def __iadd__(self, other):
        raise NotImplementedError('%s must implement __iadd__' % (type(self)))

    def __imul__(self, other):
        raise NotImplementedError('%s must implement __imul__' % (type(self)))

    def __pos__(self):
        return self

    def __add__(self, other):
        ret = self.copy()
        ret += other
        return ret

    def __radd__(self, other):
        t = type(other)
        if t is list or t is tuple:
            return Matrix(other) + self
        # other is scalar
        ret = self.copy()
        ret += other
        return ret

    def __isub__(self, other):
        self += -other
        return self

    def __sub__(self, other):
        ret = self.copy()
        ret -= other
        return ret

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        ret = self.copy()
        ret *= other
        return ret

    def __rmul__(self, other):
        t = type(other)
        if t is list or t is tuple:
            return Matrix(other) * self
        # other is scalar
        ret = self.copy()
        ret *= other
        return ret

    def __idiv__(self, other):
        self *= div(1, other)
        return self

    def __div__(self, other):
        return self * div(1, other)

    def __rdiv__(self, other):
        head, data = self.pair
        if head.is_array:
            return type(self)(head, dict([(ij,div(other,x)) for ij,x in data.iteritems()]))
        return other * self.inv()

    def __pow__(self, other):
        head, data = self.pair
        m, n = head.shape
        t = type(other)
        if t is int or t is long:
            if head.is_array:
                if other<0:
                    return type(self)(head, dict([(ij,div(1,x)**(-other)) for ij,x in data.iteritems()]))
                return type(self)(head, dict([(ij,x**other) for ij,x in data.iteritems()]))
            if head.is_diagonal:
                raise NotImplementedError(`head, str(head)`)
            if other < 0:
                return self.inv() ** (-other)
            if other==1:
                return self
            if other==0:
                return Matrix([1]*min(m,n), diagonal=True)
            if other==2:
                return self * self
            r = 1
            x = self
            while 1:
                if other & 1:
                    r *= x
                    other -= 1
                    if not other:
                        break
                x = x*x
                other //= 2
            return r
        return NotImplemented

class MatrixDict(MatrixBase):
    """ Implementation of matrix where elements are stored in a dictionary.
    """

    rows = property(lambda self: self.head.rows)
    cols = property(lambda self: self.head.cols)
    shape = property(lambda self: self.head.shape)
    is_square = property(lambda self: self.head.rows==self.head.cols)

    @property
    def is_lower(self):
        head, data = self.pair
        if head.is_transpose:
            return all(i<=j for i, j in data)
        return all(i>=j for i, j in data)

    @property
    def is_upper(self):
        head, data = self.pair
        if head.is_transpose:
            return all(i>=j for i, j in data)
        return all(i<=j for i, j in data)

    @property
    def is_orthogonal(self):
        if not self.is_square:
            return False
        return (self*self.T).is_identity

    @property
    def is_identity(self):
        head, data = self.pair
        m, n = head.shape
        if m!=n or len(data)!=n:
            return False
        return all(i==j and x==1 for (i,j),x in data.iteritems())

    @property
    def is_empty (self):
        head, data = self.pair
        m, n = head.shape
        return m == n == 0

    @property
    def is_zero(self):
        head, data = self.pair
        if not data:
            return True
        m, n = head.shape
        for (i,j),x in data.iteritems ():
            if (i<0 or i>=m): continue
            if (j<0 or j>=n): continue
            return False
        return True

    @property
    def is_row_echelon_form(self):
        """ Return True if matrix is in row echelon form, that is,

        (i) All nonzero rows (rows with at least one nonzero element)
        are above any rows of all zeroes, and
        
        (ii) The leading coefficient (the first nonzero number from
        the left, also called the pivot) of a nonzero row is always
        strictly to the right of the leading coefficient of the row
        above it.
        """
        head, data = self.pair
        m, n = head.shape
        expect_i, pivot_i, pivot_j = 0, -1, -1
        for (i,j) in sorted(data.keys()):
            if pivot_i==i:
                continue
            if expect_i!=i or pivot_j>=j:
                return False
            pivot_i, pivot_j, expect_i = i,j,i+1
        return True

    @property
    def is_row_canonical_form(self):
        """ Return True if matrix is in reduced row echelon form, that
        is,

        (i) All nonzero rows (rows with at least one nonzero element)
        are above any rows of all zeroes, and
        
        (ii) The leading coefficient (the first nonzero number from
        the left, also called the pivot) of a nonzero row is always
        strictly to the right of the leading coefficient of the row
        above it.

        (iii) Every leading coefficient is 1 and is the only nonzero
        entry in its column.
        """
        head, data = self.pair
        m, n = head.shape
        expect_i, pivot_i, pivot_j = 0,-1, -1
        for (i,j) in sorted(data.keys()):
            if pivot_i==i:
                continue
            if expect_i!=i or pivot_j>=j:
                return False
            x = data[i,j]
            if x!=1:
                return False
            if not self[:i, j].is_zero:
                return False
            pivot_i, pivot_j, expect_i = i,j,i+1
        return True

    @property
    def T(self):
        """ Return transposed view of a matrix.
        """
        head, data = self.pair
        return type(self)(head.T, data)

    @property
    def A(self):
        """ Return array view of a matrix.
        """
        head, data = self.pair
        if head.is_array:
            return self
        newhead = head.A
        if newhead is head:
            return self
        return type(self)(newhead, data)

    @property
    def M(self):
        """ Return matrix view of a matrix.
        """
        head, data = self.pair
        if not (head.is_array or head.is_diagonal):
            return self
        newhead = head.M
        if newhead is head:
            return self
        return type(self)(newhead, data)

    @property
    def D(self):
        """ Return diagonal view of a matrix.
        """
        head, data = self.pair
        if head.is_diagonal:
            return self
        newhead = head.D
        if newhead is head:
            return self
        return type(self)(newhead, data)

    @property
    def I(self):
        """ Return inverse matrix.
        """
        return self.inv()

    def resize(self, m, n):
        """ Return resized view of a matrix.

        Notes:

         - Use crop() to remove elements out-of-matrix bounds.
         - If matrix is not writable then matrix data is copied.
        """
        head, data = self.pair
        p, q = head.shape
        if p==m and q==n:
            return self
        newhead = MATRIX(m, n, head.storage)
        if not self.is_writable:
            data = dict(data)
        return type(self)(newhead, data)

    def tolist(self):
        """Convert matrix to a list of lists."""
        rows, cols = self.head.shape
        return [[self[i,j] for j in xrange(cols)] for i in xrange(rows)]
    
    def flatten(self):
        rows, cols = self.head.shape
        for i in xrange (rows):
            for j in xrange (cols):
                yield self[i,j]

    def _get_diagonal(self, key):
        head, data = self.pair
        rows, cols = head.shape
        t = type(key)
        if t is int or t is long:
            if key<0:
                m = min(rows, cols) + key
            else:
                m = min(rows, cols) - key
            d = {}
            if head.is_transpose:
                key = -key
            for (i,j),x in data.items():
                if j-i==key:
                    k = min(i,j)
                    d[k,0] = x
            return Matrix(m, 1, d)
        elif is_integer(key):
            return self._get_diagonal[int(key)]
        raise NotImplementedError(`head, key, t`)

    def __getitem__(self, key):
        tkey = type(key)
        head, data = self.pair
        if head.is_diagonal:
            return self._get_diagonal(key)
        if tkey is tuple:
            i, j = key
            ti, tj = type(i), type(j)
            if ti is int and i<0:
                i += head.rows
            if tj is int and j<0:
                j += head.cols
            if ti is int and tj is int:
                if head.is_transpose:
                    key = j, i
                elif head.is_diagonal:
                    raise NotImplementedError(`head, str(head)`)
                else:
                    key = i, j
                if i>=head.rows or j>=head.cols:
                    raise IndexError (`i,j,head.cols,head.rows`)
                return data.get(key, 0)
            if ti is slice:
                row_indices = dict([(i0,k) for k,i0 in enumerate(xrange(*i.indices(head.rows)))])
            elif ti is int:
                row_indices = {i:0}
            elif ti is tuple:
                row_indices = dict([(i0,k) for k,i0 in enumerate(i)])                
            elif is_integer(i):
                return self[int(i), j]
            else:
                raise IndexError('first index must contain int or slice or tuple, got %s' % ((`ti`)))
            if tj is slice:
                col_indices = dict([(j0,k) for k,j0 in enumerate(xrange(*j.indices(head.cols)))])
            elif tj is int:
                col_indices = {j:0}
            elif tj is tuple:
                col_indices = dict([(j0,k) for k,j0 in enumerate(j)])
            elif is_integer(j):
                return self[i, int (j)]
            else:
                raise IndexError('second index must contain int or slice or tuple, got %s' % ((`tj`)))
            newdata = {}
            if head.is_transpose:
                for (j,i), x in data.items():
                    ki = row_indices.get(i)
                    if ki is not None:
                        kj = col_indices.get(j)
                        if kj is not None:
                            newdata[ki, kj] = x
            elif head.is_diagonal:
                raise NotImplementedError(`head, str(head)`)
            else:
                for (i,j), x in data.items():
                    ki = row_indices.get(i)
                    if ki is not None:
                        kj = col_indices.get(j)
                        if kj is not None:
                            newdata[ki, kj] = x
            return Matrix(len(row_indices), len(col_indices), newdata)
        elif tkey is int:
            if key>=head.rows:
                raise IndexError (`key, head.rows`)
            return self[key, :]
        elif tkey is slice:
            return self[key, :]
        elif tkey is tuple:
            return self[key, :]
        elif is_integer(key):
            return self[int(key)]

        raise IndexError('index must be int, slice, or tuple, got %s' % (tkey))

    def _set_diagonal(self, key, value):
        head, data = self.pair
        assert type(key) in [int, long],`type(key)`
        rows, cols = head.shape
        if head.is_transpose:
            key = -key
        if key<0:
            k = min(rows, cols) + key
        else:
            k = min(rows, cols) - key
        if isinstance(value, MatrixBase):
            m, n = value.head.shape
            assert n==1,`m,n`
            assert m==k,`m,k,rows,cols`
            if key<0:
                for (i,j),x in value.data.items():
                    data[i-key,i] = x
            else:
                for (i,j),x in value.data.items():
                    data[i,i+key] = x
        else:
            if key < 0:
                for i in range(k):
                    data[i-key, i] = value
            else:
                for i in range(k):
                    data[i, i+key] = value

    def __setitem__(self, key, value):
        if not self.is_writable:
            raise TypeError('Matrix content is read-only')
        if isinstance(value, list):
            self[key] = Matrix(value)
            return

        head, data = self.pair
        if head.is_diagonal:
            return self._set_diagonal(key, value)
        tkey = type(key)
        if tkey is int or tkey is slice:
            if head.is_diagonal:
                raise NotImplementedError(`head, str(head)`)
            self[key, :] = value
            return
        i, j = key
        ti, tj = type(i), type(j)
        if ti is int and i<0:
            i += head.rows
        if tj is int and j<0:
            j += head.cols
        if ti is int and tj is int:
            if head.is_transpose:
                key = j, i
            elif head.is_diagonal:
                raise NotImplementedError(`head, str(head)`)
            else:
                key = i, j
            if value:
                data[key] = value
            else:
                try:
                    del data[key]
                except KeyError:
                    pass
            return
        if ti is slice and tj is int:
            row_indices = [(i0,k) for k,i0 in enumerate(xrange(*i.indices(head.rows)))]
            col_indices = [(j,0)]
        elif ti is int and tj is slice:
            row_indices = [(i,0)]
            col_indices = [(j0,k) for k,j0 in enumerate(xrange(*j.indices(head.cols)))]
        elif ti is slice and tj is slice:
            row_indices = [(i0,k) for k,i0 in enumerate(xrange(*i.indices(head.rows)))]
            col_indices = [(j0,k) for k,j0 in enumerate(xrange(*j.indices(head.cols)))]
        elif is_integer(i) and tj in (int, slice):
            
            return self.__setitem__((int(i),j), value)
        elif is_integer(j) and ti in (int, slice):
            return self.__setitem__((i,int(j)), value)
        else:
            raise IndexError('tuple index must contain int or slice but got %r' % ((ti, tj),))
        if isinstance(value, MatrixBase):
            m, n = value.head.shape
            assert len(row_indices)==m,`len(row_indices),m`
            assert len(col_indices)==n,`len(col_indices),n`
            if head.is_transpose:
                for i,ki in row_indices:
                    for j,kj in col_indices:
                        v = value[ki, kj]
                        if v:
                            data[j, i] = v
                        else:
                            try:
                                del data[j, i]
                            except KeyError:
                                pass
            elif head.is_diagonal:
                raise NotImplementedError(`head, str (head)`)
            else:
                for i,ki in row_indices:
                    for j,kj in col_indices:
                        v = value[ki, kj]
                        if v:
                            data[i, j] = v
                        else:
                            try:
                                del data[i, j]
                            except KeyError:
                                pass
        else:
            if value:
                if head.is_transpose:
                    for i,ki in row_indices:
                        for j,kj in col_indices:
                            data[j, i] = value
                elif head.is_diagonal:
                    raise NotImplementedError(`head`)
                else:
                    for i,ki in row_indices:
                        for j,kj in col_indices:
                            data[i, j] = value
            else:
                if head.is_transpose:
                    for i,ki in row_indices:
                        for j,kj in col_indices:
                            try:
                                del data[j, i]
                            except KeyError:
                                pass
                elif head.is_diagonal:
                    raise NotImplementedError(`head`)
                else:
                    for i,ki in row_indices:
                        for j,kj in col_indices:
                            try:
                                del data[i, j]
                            except KeyError:
                                pass
    
    def copy(self):
        """ Return a copy of a matrix.
        """
        head, data = self.pair
        return MatrixDict(head, dict(data))


    def __array__(self):
        """ Return matrix as numpy array.
        """
        import numpy
        head, data = self.pair
        storage = head.storage
        shape = head.shape
        r = numpy.zeros(shape,dtype=object)
        if head.is_transpose:
            for (j, i),e in data.items():
                r[i][j] = e
        elif head.is_diagonal:
            raise NotImplementedError(`head`)
        else:
            for (i,j),e in data.items():
                r[i][j] = e
        return r

    def __neg__(self):
        head, data = self.pair
        return MatrixDict(head, dict([(index,-element) for index, element in data.items()]))

    def inv_l(self):
        """ Return inverse of lower triangular matrix L.
        """

        M,N = self.head.shape
        linv = Matrix(M,N)
        assert M==N,`M,N`
        for j in xrange(N):
            linv[j,j] = div(1, self[j, j])
            for i in xrange(j+1,N):
                r = 0
                for k in range(j,i):
                    r = r - self[i, k] * linv[k, j]
                linv[i, j] = div(r, self[i, i])
        return linv

    def inv(self):
        """ Return inverse of a square matrix.
        """
        head, data = self.pair
        m, n = head.shape
        assert m==n,`m,n`
        if head.is_transpose:
            d = {}
            for (i,j),x in data.items():
                d[j,i] = x
        elif head.is_diagonal:
            raise NotImplementedError(`head`)
        else:
            d = dict(data)
        for i in range(m):
            d[i,i+m] = 1
        a = MatrixDict(MATRIX(m, 2*m, MATRIX_DICT), d)
        # XXX: catch singular matrices
        b = a.gauss_jordan_elimination(overwrite=True)
        return b[:,m:]

    def solve(self, rhs):
        """ Solve a system of linear equations A * x = rhs

        where A is square matrix. Alternative usage::

          A // rhs -> x
        
        For example::

          Matrix([[1,2], [3,4]]) // [1,2] -> Matrix([[0],[1/2]])
        
        See also
        --------
        solve_null

        """
        t = type(rhs)
        if t is tuple or t is list:
            rhs = Matrix(rhs)
        head, data = self.pair
        m, n = head.shape
        assert m==n,`m,n`
        if head.is_transpose:
            d = {}
            for (i,j), x in data.items():
                d[j,i] = x
        elif head.is_diagonal:
            raise NotImplementedError(`head`)
        else:
            d = dict(data)
        rhead, rdata = rhs.pair
        p, q = rhead.shape
        assert p==m,`p,m`
        if rhead.is_transpose:
            for (j,i),x in rdata.items():
                d[i,j+n] = x
        else:
            for (i,j),x in rdata.items():
                d[i,j+n] = x
        a = MatrixDict(MATRIX(m, n+q, MATRIX_DICT), d)
        b = a.gauss_jordan_elimination(overwrite=True)
        return b[:,m:]

    __floordiv__ = solve

    def normal(self):
        d = {}
        head, data = self.pair
        for ij, x in data.iteritems():
            if hasattr(x, 'normal'):
                d[ij] = x.normal()
            else:
                d[ij] = x
        return type(self)(head, d)

    def expand(self):
        d = {}
        head, data = self.pair
        for ij, x in data.iteritems():
            if hasattr(x, 'expand'):
                d[ij] = x.expand()
            else:
                d[ij] = x
        return type(self)(head, d)

    def subs (self, subexpr, newexpr = None):
        data = {}
        for key, value in self.data.items():
            if isinstance (value, classes.Expr):
                value = value.subs(subexpr, newexpr)
            if value:
                data[key] = value
        return type(self)(self.head, data)

    @property
    def rank(self):
        """ Return row-rank of a matrix.
        """
        return self.gauss_jordan_elimination().rows

    @property
    def nullity (self):
        return self.cols - self.rank

    def solve_null(self, labels = None, check=False):
        """ Solve a system of homogeneous equations: A * x = 0

        where A is m x n matrix, x is n vector.

        Parameters
        ----------
        labels : {list, tuple}
          A list of column symbols.
        check : bool
          When True the check the definitions and raise a runtime
          exception on invalid results.

        Returns
        -------
        xd : dict
          A dictionary of null solutions. xd keys are labels of
          columns and xd values are rows of A null space basis
          matrix. That is, let us define::

            ker = Matrix([xd[l].tolist()[0] for l in labels])

          then::

            A * ker = 0

          or equivalently, if

            x = ker * Matrix(labels)

          then

            A * x = 0.

        dependent, independent : list
          Lists of dependent and independent column labels, respectively.

          If
            x_dep = [sum(dot(xd[l], independent)) for l in labels]
          then
            A * x_dep = 0

        See also
        --------
        solve
        """
        m, n = self.head.shape
        if labels is None:
            labels = range(n)
        gj, (dep, indep) = self.gauss_jordan_elimination(labels = labels)
        xd = {}
        rank = gj.rows
        nullity = n - rank
        for i, label in enumerate (dep):
            xd[label] = -gj[i, rank:]
        for i, label in enumerate (indep):
            xd[label] = Matrix(1,nullity, {(0,i):1})
        if check:
            ker = Matrix([xd[l].tolist () [0] for l in labels])
            null = self * ker
            if not null.is_zero:
                raise RuntimeError('%s.solve_null: A*ker=0 test failed' % (self.__class__.__name__))
            c = Matrix([classes.Calculus('X%s'%(label)) for label in indep])
            x = ker * c
            null = self*x
            if not null.is_zero:
                raise RuntimeError('%s.solve_null: A*(ker*indep)=0 test failed' % (self.__class__.__name__))

        return xd, dep, indep

    def gram_schmidt_process(self):
        """ Orthogonalized matrix columns using Gram-Schmidt process.

        Returns
        -------
        U : Matrix

        Notes
        -----
        U columns are orthogonal (but not orthonormal). U * self
        produces upper triangular matrix.
        """
        m, n = self.head.shape
        U = Matrix(m,n)
        k1 = 0
        l = []
        for k in range(n):
            ak = self[:,k]
            uk = ak
            flag = True
            for j in range(k1):
                uj = U[:,j]
                cdenom = (uj.T*uj)[0,0]
                c = div((uj.T*ak)[0,0], cdenom)
                uk = uk - c * uj
            if not uk.is_zero:
                U[:,k1] = uk
                k1 += 1
                l.append(k)
        U = U.resize(m, k1)
        return U, l

def switch_dep_indep_variables (xd, dep, indep, dep_label, indep_label):
    m = rank = len (dep)
    nullity = len (indep)
    n = m + nullity
    u = Matrix(m,n)
    for i in range(m):
        u[i,i] = 1
    for i,label in enumerate(dep):
        u[i,m:] = -xd[label]
    j1 = dep.index (dep_label)
    j2 = m + indep.index (indep_label)
    dep, indep = dep[:], indep[:]
    dep[j1] = indep_label
    indep[j2-m] = dep_label
    u.swap_cols (j1, j2)
    gj = u[:].gauss_jordan_elimination(overwrite=True)
    new_xd = {}
    for i in range (m):
        label = dep[i]
        new_xd[label] = -gj[i, rank:]
    for i in range (nullity):
        label = indep[i]
        new_xd[label] = Matrix(1, nullity, {(0,i):1})
    return new_xd, dep, indep

from .matrix_operations import MATRIX_DICT_iadd, MATRIX_DICT_imul
from .linalg import (MATRIX_DICT_swap_rows, MATRIX_DICT_swap_cols,
                     MATRIX_DICT_lu, MATRIX_DICT_crop,
                     MATRIX_DICT_gauss_jordan_elimination,
                     MATRIX_DICT_trace,
                     MATRIX_DICT_get_gauss_jordan_elimination_operations,
                     MATRIX_DICT_apply_row_operations)
from .linalg_determinant import MATRIX_DICT_determinant
from .linalg_lp import MATRIX_DICT_LP_solve

MatrixDict.__iadd__ = MATRIX_DICT_iadd
MatrixDict.__imul__ = MATRIX_DICT_imul
MatrixDict.swap_rows = MATRIX_DICT_swap_rows
MatrixDict.swap_cols = MATRIX_DICT_swap_cols
MatrixDict.crop = MATRIX_DICT_crop
MatrixDict.lu = MATRIX_DICT_lu
MatrixDict.gauss_jordan_elimination = MATRIX_DICT_gauss_jordan_elimination
MatrixDict.det = MATRIX_DICT_determinant
MatrixDict.trace = MATRIX_DICT_trace
MatrixDict.get_gauss_jordan_elimination_operations = MATRIX_DICT_get_gauss_jordan_elimination_operations
MatrixDict.apply_row_operations = MATRIX_DICT_apply_row_operations
MatrixDict.LP_solve = MATRIX_DICT_LP_solve
