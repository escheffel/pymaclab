# Author: Pearu Peterson
# Created: March 2008

import sys

from ..utils import MATRIX, MATRIX_DICT, MATRIX_DICT_T
from ..arithmetic.numbers import div
from .algebra import MatrixDict, Matrix

from ..core import init_module
init_module.import_lowlevel_operations()

def MATRIX_DICT_get_gauss_jordan_elimination_operations(self, overwrite=False, leading_cols = None, trailing_cols = None,
                                                        leading_row_selection = None, leading_column_selection = None,
                                                        verbose = False):
    """ Compute the row operations of Gauss-Jordan elimination of a m x n matrix A.

    Parameters
    ----------
    overwrite : bool
      When True then discard the content of matrix A.

    leading_cols : {None, list}
      A list of column indices that defines the order of leading
      columns used for Gauss-Jordan elimination.  By default,
      leading_cols = range(n).

    trailing_cols : {None, list}
      A list of column indices that are selected as leading columns
      when the leading_cols list has been consumed.

    leading_row_selection : {None, 'sparsest first'}
      A rule to select the leading row index. By default, the first
      valid (unused and non-zero entity at leading column) row is
      choosen as a leading row.

    leading_column_selection : {None, 'sparsest first'}
      A rule to select the leading column index. By default, the order
      of leading columns is defined by the leading_cols list.

    Returns
    -------
    B : MatrixDict
      m x n matrix

    row_operations : list
    
      A list of row operations. A row operation is represented as
      tuple: 

        ((k,j), ) denotes swapping of the k-th and j-th row;
        (k, c) denotes multiplication of k-th row by scalar c;
        (k, c, j) denotes adding j-th row multiplied by c to k-th row.

    leading_rows : list
      A list of leading row indices.

    leading_cols : list
      A list of leading column indices. leading_cols list is
      equivalent to a list of dependent variables.

    zero_rows : list
      A list of row indices that correspond to rows that Gauss-Jordan
      elimination zerod out or were zero, that is, B[i,:]==0 for i in
      zero_rows. The rest of the B rows are linearly independent.

    Notes
    -----
    
    The used algorithm does not swap rows nor columns. One can use
    leading_rows information to turn the matrix B to row-echelon
    form. However, then also the data in row_operations must be
    updated accordingly.

    See also
    --------
    apply_row_operations, gauss_jordan_elimination
    """
    head, data = self.pair
    m, n = head.shape
    k = min(m,n)
    if overwrite and self.is_writable:
        new_data = data
    else:
        new_data = dict(data)
    if leading_cols is None:
        leading_cols = range(n)
    if trailing_cols is None:
        trailing_cols = []
    if head.is_transpose:
        ops, leading_rows, leading_cols, zero_rows = \
            get_gauss_jordan_elimination_operations_MATRIX_T(m, n, new_data, leading_cols, trailing_cols, 
                                                             leading_row_selection, leading_column_selection,
                                                             verbose)
    elif head.is_diagonal:
        raise NotImplementedError(`head, head.is_diagonal`)
    else:
        ops, leading_rows, leading_cols, zero_rows = \
            get_gauss_jordan_elimination_operations_MATRIX(m, n, new_data, leading_cols, trailing_cols, 
                                                           leading_row_selection, leading_column_selection, 
                                                           verbose)
    if new_data is data:
        B = self
    else:
        B = MatrixDict(head, new_data)
    return B, ops, leading_rows, leading_cols, zero_rows

def MATRIX_DICT_apply_row_operations(self, row_operations, overwrite=False):
    """ Apply row operations to a m x n matrix A.

    Parameters
    ----------
    row_operations : list
    
      A list of row operations. A row operation is represented as tuple:

      ((k,j), ) denotes swapping of the k-th and j-th row;
      (k, c) denotes multiplication of k-th row by scalar c;
      (k, c, j) denotes adding j-th row multiplied by c to k-th row.

    overwrite : bool
      When True then discard the content of matrix A.

    Results
    -------
    B : MatrixDict
      m x n matrix

    See also
    --------
    get_gauss_jordan_elimination_operations, gauss_jordan_elimination
    
    """
    head, data = self.pair
    m, n = head.shape
    if overwrite and self.is_writable:
        new_data = data
    else:
        new_data = dict(data)    
    if head.is_transpose:
        apply_row_operations_MATRIX_T(m, n, new_data, row_operations)
    elif head.is_diagonal:
        raise NotImplementedError(`head, head.is_diagonal`)
    else:
        apply_row_operations_MATRIX(m, n, new_data, row_operations)
    if new_data is data:
        return self
    return MatrixDict(head, new_data)

def MATRIX_DICT_gauss_jordan_elimination(self, swap_columns=False, overwrite=False, labels = None, return_pivot_info = False):
    """ Perform Gauss-Jordan elimination of a m x n matrix A.

    Parameters
    ----------

    swap_columns : bool

      When True then ensure that row echelon form has maximum number
      of nonzero diagonal elements by swapping columns.

    overwrite : bool
      When True then discard the content of matrix A.

    labels : {None,list}
      A list of column labels.

    return_pivot_info : bool
      When True, return pivot information, see below.

    Returns
    -------
    B : MatrixDict
      m1 x n matrix that is in row echelon form. m1 is the number of non-zero rows.

    pivot_table : list
      A n-list of column indices. pivot_table is returned only when
      swap_columns is True and labels is not specified.

    (dep, indep) : tuple
      A 2-tuple of column labels corresponding to dependent and
      independent variables. The 2-tuple is returned only when labels
      is specified.

    (row_pivot_table, column_pivot_table) : tuple
      A 2-tuple of lists corresponding to row and column pivot tables.
      The 2-tuple is returned only when return_pivot_info is True.

    See also
    --------
    get_gauss_jordan_elimination_operations, apply_row_operations
    """
    head, data = self.pair
    m, n = head.shape
    if overwrite and self.is_writable:
        new_data = data
    else:
        new_data = dict(data)
    if labels:
        swap_columns = True
    if head.is_transpose:
        m, row_pivot_table, pivot_table = gauss_jordan_elimination_MATRIX_T(m, n, new_data, swap_columns = swap_columns)
        new_head = MATRIX(m, n, MATRIX_DICT_T)
    elif head.is_diagonal:
        raise NotImplementedError(`head, head.is_diagonal`)
    else:
        m, row_pivot_table, pivot_table = gauss_jordan_elimination_MATRIX(m, n, new_data, swap_columns = swap_columns)
        new_head = MATRIX(m, n, MATRIX_DICT)
    if data is new_data:
        B = self
    else:
        B = MatrixDict(new_head, new_data)

    if labels:
        dep = [labels[pivot_table[i]] for i in range (B.rows)]
        indep = [labels[pivot_table[i]] for i in range (B.rows, B.cols)]
        if return_pivot_info:
            return B, (dep, indep), (row_pivot_table, pivot_table)
        return B, (dep, indep)
    if swap_columns:
        if return_pivot_info:
            return B, (row_pivot_table, pivot_table)
        return B, pivot_table
    if return_pivot_info:
        return B, (row_pivot_table, pivot_table)
    return B

def MATRIX_DICT_lu(self, overwrite=False):
    """ Perform LU factorization of a m x n matrix A.

    Parameters
    ----------

      overwrite : bool
        When True then discard the content of matrix A.

    Returns
    -------

      P : MatrixDict
      L : MatrixDict
      U : MatrixDict
        LU decomposition matrices of A.

    Notes
    -----
        
      P - m x m permuation matrix
      L - m x k lower triangular or trapezoidal matrix with unit-diagonal
      U - k x n upper triangular or trapezoidal matrix
      k = min(m,n)
      
      A = P * L * U
    """
    head, data = self.pair
    m, n = head.shape
    k = min(m, n)
    ldata = {}
    if overwrite and self.is_writable:
        udata = data
    else:
        udata = dict(data)
    L = MatrixDict(MATRIX(m, k, MATRIX_DICT), ldata)
    if head.is_transpose:
        U = MatrixDict(MATRIX(k, n, MATRIX_DICT_T), udata)
        pivot_table = lu_MATRIX_T(m, n, k, ldata, udata)
    elif head.is_diagonal:
        raise NotImplementedError(`head`)
    else:
        U = MatrixDict(MATRIX(k, n, MATRIX_DICT), udata)
        pivot_table = lu_MATRIX(m, n, k, ldata, udata)
    P = Matrix(pivot_table, permutation=True).T
    return P, L, U

def get_gauss_jordan_elimination_operations_MATRIX(m, n, data, leading_cols0, trailing_cols, 
                                                   leading_row_selection,
                                                   leading_column_selection, 
                                                   verbose):
    data_get = data.get
    data_has = data.has_key
    rows, cols = get_rc_maps(data)
    leading_rows = []
    leading_cols = []
    zero_rows = sorted(set(range(m)).difference(rows))
    ops = []
    free_cols = [i for i in leading_cols0+trailing_cols if i<n]
    n0 = len(free_cols)
    select_sparsest_leading_row = leading_row_selection=='sparsest first'
    select_sparsest_leading_column = leading_column_selection=='sparsest first'
    for i0 in range(n0):
        if verbose:
            sys.stdout.write('\rget_gauss_jordan_elimination_operations: %5.2f%%' % (100.0*i0/n0))
            sys.stdout.flush()

        if select_sparsest_leading_column:
            lmin, i = m, n
            lmin0, i0 = m, n
            for i1 in free_cols:
                l = len (cols[i1])
                if i1 in trailing_cols:
                    if l < lmin0:
                        lmin0, i0 = l, i1
                else:
                    if l < lmin:
                        lmin, i = l, i1
            if i==n:
                i = i0
            if i==n:
                break
            free_cols.remove(i)
        else:
            i = free_cols[i0]

        if select_sparsest_leading_row:
            lmin, k = n, m
            for k1 in cols[i]:
                if k1 in leading_rows:
                    continue
                l = len(rows[k1])
                if l < lmin:
                    lmin, k = l, k1
        else:
            k = 0
            while not data_has((k, i)) or k in leading_rows:
                k += 1
                if k==m:
                    break

        a_ki = data_get((k, i))

        if a_ki is None:
            continue

        leading_rows.append(k)
        leading_cols.append(i)

        krow = rows[k]
        for j in range(m):
            if j==k:
                continue
            u_ji = data_get((j,i))
            if u_ji is None:
                continue
            c = div(-u_ji, a_ki)
            ops.append((j, c, k)) # add c times k-th row to j-th row
            jrow = rows[j]
            jrow_add = jrow.add
            jrow_remove = jrow.remove
            for p in krow:
                pcol = cols[p]
                b = data_get((j,p))
                if b is None:
                    data[j,p] = data[k,p] * c
                    jrow_add(p)
                    pcol.add(j)
                else:
                    b = b + data[k,p] * c
                    if b:
                        data[j,p] = b
                    else:
                        del data[j,p]
                        jrow_remove (p)
                        pcol.remove(j)
                if not pcol:
                    del cols[p]
            if not jrow:
                zero_rows.append(j)
                del rows[j]


        if a_ki==1:
            continue
        ia_ki = div(1,a_ki)
        ops.append((k, ia_ki)) # multiply k-th row with 1/a_ki
        for p in krow:
            data[k,p] *= ia_ki
    if verbose:
        sys.stdout.write ('\n')
    return ops, leading_rows, leading_cols, zero_rows

def get_gauss_jordan_elimination_operations_MATRIX_T(m, n, data, leading_cols0, trailing_cols, 
                                                     leading_row_selection, leading_column_selection, 
                                                     verbose):
    data_get = data.get
    data_has = data.has_key
    cols, rows = get_rc_maps(data)
    leading_rows = []
    leading_cols = []
    zero_rows = sorted(set(range(m)).difference(rows))
    ops = []
    free_cols = [i for i in leading_cols0 if i<n]
    select_sparsest_leading_row = leading_row_selection=='sparsest first'
    select_sparsest_leading_column = leading_column_selection=='sparsest first'
    n0 = len(free_cols)
    for i0 in range(n0):
        if verbose:
            sys.stdout.write('\rget_gauss_jordan_elimination_operations[T]: %5.2f%%' % (100.0*i0/n0))
            sys.stdout.flush()

        if select_sparsest_leading_column:
            lmin, i = m, n
            lmin0, i0 = m, n
            for i1 in free_cols:
                l = len (cols[i1])
                if i1 in trailing_cols:
                    if l < lmin0:
                        lmin0, i0 = l, i1
                else:
                    if l < lmin:
                        lmin, i = l, i1
            if i==n:
                i = i0
            free_cols.remove(i)
        else:
            i = free_cols[i0]

        if select_sparsest_leading_row:
            lmin, k = n, m
            for k1 in cols[i]:
                if k1 in leading_rows:
                    continue
                l = len(rows[k1])
                if l < lmin:
                    lmin, k = l, k1
        else:
            k = 0
            while not data_has((i, k)) or k in leading_rows:
                k += 1
                if k==m:
                    break

        a_ki = data_get((i, k))
        if a_ki is None:
            continue

        leading_rows.append(k)
        leading_cols.append(i)

        krow = rows[k]
        for j in range(m):
            if j==k:
                continue
            u_ji = data_get((i,j))
            if u_ji is None:
                continue
            c = div(u_ji, a_ki)
            ops.append((j, -c, k)) # add (-c) times k-th row to j-th row
            jrow = rows[j]
            jrow_add = jrow.add
            jrow_remove = jrow.remove
            for p in krow:
                pcol = cols[p]
                u_kp_c = data[p,k] * c
                jp = p,j
                b = data_get(jp)
                if b is None:
                    data[jp] = -u_kp_c
                    jrow_add(p)
                    pcol.add(j)
                else:
                    if u_kp_c==b:
                        del data[jp]
                        jrow_remove(p)
                        pcol.remove(j)
                    else:
                        data[jp] = b - u_kp_c
                if not pcol:
                    del cols[p]
            if not jrow:
                zero_rows.append(j)
                del rows[j]
        if a_ki==1:
            continue
        ia_ki = div(1,a_ki)
        ops.append((k, ia_ki)) # multiply k-th row with 1/a_ki
        for p in krow:
            data[p,k] *= ia_ki
    return ops, leading_rows, leading_cols, zero_rows


def apply_row_operations_MATRIX(m, n, data, ops):
    for op in ops:
        if len (op)==1:
            k, j = op[0]
            swap_rows_MATRIX(data, k, j)            
        elif len(op)==2:
            # multiply k-th row with c
            k, c = op
            for i, j in data:
                if i==k:
                    data[i,j] *= c
        elif len(op)==3:
            # add c times k-th row to j-th row
            p, c, k = op
            assert p!=k,`op`
            for i, j in data.keys():
                if i==k:
                    ck = data[k,j] * c
                    pj = p,j
                    a_pj = data.get(pj)
                    if a_pj is None:
                        data[pj] = ck
                    else:
                        a_pj += ck
                        if a_pj:
                            data[pj] = a_pj
                        else:
                            del data[pj]
        else:
            raise ValueError('unknown row operation: %r' % (ops,))

def apply_row_operations_MATRIX_T(m, n, data, ops):
    for op in ops:
        if len (op)==1:
            k, j = op[0]
            swap_rows_MATRIX_T(data, k, j)
        elif len(op)==2:
            # multiply k-th row with c
            k, c = op
            for j,i in data:
                if i==k:
                    data[j,i] *= c
        elif len(op)==3:
            # add c times k-th row to j-th row
            p, c, k = op
            assert p!=k,`op`
            for j,i in data.keys():
                if i==k:
                    ck = data[j,k] * c
                    pj = j,p
                    a_pj = data.get(pj)
                    if a_pj is None:
                        data[pj] = ck
                    else:
                        a_pj += ck
                        if a_pj:
                            data[pj] = a_pj
                        else:
                            del data[pj]
        else:
            raise ValueError('unknown row operation: %r' % (ops,))


def gauss_jordan_elimination_MATRIX(m, n, data, swap_columns = False):
    data_get = data.get
    data_has = data.has_key
    rows = get_rc_map(data)
    row_pivot_table = range(m)
    if swap_columns:
        pivot_table = range(n)
    else:
        pivot_table = None
    jpiv = 0    
    for i in xrange(m):
        ipiv = i

        while not data_has((ipiv, jpiv)):
            ipiv += 1
            if ipiv==m:
                ipiv = i
                jpiv += 1
                if jpiv==n:
                    break

        a_ii = data_get((ipiv, jpiv))
        if a_ii is None:
            break

        if swap_columns:
            if jpiv != i:
                swap_cols_MATRIX(data, jpiv, i)
                pivot_table[i], pivot_table[jpiv] = pivot_table[jpiv], pivot_table[i]
                rows = get_rc_map(data)
                jpiv = i

        if i!=ipiv:
            swap_rows_MATRIX(data, i, ipiv, rows)
            row_pivot_table[i], row_pivot_table[ipiv] = row_pivot_table[ipiv], row_pivot_table[i]
            swap_rc_map(rows, i, ipiv)
            ipiv = i
        irow = rows[i]
        for j in range(m):
            if j==i:
                continue
            u_ji = data_get((j,jpiv))
            if u_ji is None:
                continue
            c = div(u_ji, a_ii)
            jrow = rows[j]
            jrow_add = jrow.add
            jrow_remove = jrow.remove
            for p in irow:
                if p < jpiv: 
                    continue
                u_ip_c = data[i,p] * c
                jp = j,p
                b = data_get(jp)
                if b is None:
                    data[jp] = -u_ip_c
                    jrow_add(p)
                else:
                    if u_ip_c==b:
                        del data[jp]
                        jrow_remove(p)
                    else:
                        data[jp] = b - u_ip_c
            if not jrow:
                del rows[j]
        data[i,jpiv] = 1
        
        for p in irow:
            if p <= jpiv: 
                continue
            ip = i,p
            data[ip] = div(data[ip], a_ii)
    if rows:
        return max(rows)+1,row_pivot_table, pivot_table
    return 0, row_pivot_table, pivot_table

def gauss_jordan_elimination_MATRIX_T(m, n, data, swap_columns = False):
    data_get = data.get
    data_has = data.has_key
    rows = get_rc_map_T(data)
    row_pivot_table = range(m)
    if swap_columns:
        pivot_table = range(n)
    else:
        pivot_table = None
    jpiv = 0
    for i in xrange(m):
        ipiv = i

        while not data_has((jpiv, ipiv)):
            ipiv += 1
            if ipiv==m:
                ipiv = i
                jpiv += 1
                if jpiv==n:
                    break
        a_ii = data_get((jpiv, ipiv))
        if a_ii is None:
            break

        if swap_columns:
            if jpiv != i:
                swap_cols_MATRIX_T(data, jpiv, i)
                pivot_table[i], pivot_table[jpiv] = pivot_table[jpiv], pivot_table[i]
                rows = get_rc_map_T(data)
                jpiv = i

        if i!=ipiv:
            swap_rows_MATRIX_T(data, i, ipiv, rows)
            row_pivot_table[i], row_pivot_table[ipiv] = row_pivot_table[ipiv], row_pivot_table[i]
            swap_rc_map(rows, i, ipiv)
            ipiv = i
        irow = rows[i]
        for k in range(m):
            if k==ipiv:
                continue
            u_ji = data_get((jpiv,k))
            if u_ji is None:
                continue
            c = div(u_ji, a_ii)
            jrow = rows[k]
            jrow_add = jrow.add
            jrow_remove = jrow.remove
            for p in irow:
                if p < jpiv: 
                    continue
                u_ip_c = data[p,i] * c
                kp = p,k
                b = data_get(kp)
                if b is None:
                    data[kp] = -u_ip_c
                    jrow_add(p)
                else:
                    if u_ip_c==b:
                        del data[kp]
                        jrow_remove(p)
                    else:
                        data[kp] = b - u_ip_c
            if not jrow:
                del rows[k]
        data[jpiv,i] = 1
        
        for p in irow:
            if p <= jpiv: 
                continue
            ip = p,i
            data[ip] = div(data[ip], a_ii)

    if rows:
        return max(rows)+1, row_pivot_table, pivot_table
    return 0, row_pivot_table, pivot_table

def get_rc_map(data):
    """ Return mapping between row indices and column indices with existing values.
    """
    rows = {}
    for i, j in data:
        s = rows.get (i)
        if s is None:
            s = rows[i] = set ()
        s.add(j)
    return rows

def get_rc_map_T(data):
    rows = {}
    for j, i in data:
        s = rows.get (i)
        if s is None:
            s = rows[i] = set ()
        s.add(j)
    return rows

def get_rc_maps(data):

    rows = {}
    cols = {}
    for i, j in data:
        s = rows.get(i)
        if s is None:
            s = rows[i] = set ()
        s.add(j)
        s = cols.get(j)
        if s is None:
            s = cols[j] = set()
        s.add(i)
    return rows, cols

def swap_rc_map(rows, i, j):
    ri = rows.get(i)
    rj = rows.get(j)
    if rj is not None: del rows[j]
    if ri is not None: del rows[i]
    if rj is not None:
        rows[i] = rj
    if ri is not None:
        rows[j] = ri

def lu_MATRIX(m, n, k, ldata, udata):
    if m>n:
        for i in xrange(n,m):
            udata[(i,i)] = 1
    pivot_table = range(m)
    udata_get = udata.get
    udata_has = udata.has_key

    urows = get_rc_map(udata)
    lrows = get_rc_map(ldata)
    for i in xrange(m-1):
        ncols, j = n, None
        for j1 in range(i,m):
            if udata_has((j1,i)):
                l = len(urows[j1])
                if l <= ncols:
                    ncols, j = l, j1
        if j is None:
            continue
        a_ii = udata[j,i]
        if i!=j:
            pivot_table[i], pivot_table[j] = pivot_table[j], pivot_table[i]
            swap_rows_MATRIX(udata, i, j, urows)
            swap_rows_MATRIX(ldata, i, j, lrows)
            swap_rc_map(urows, i, j)
            swap_rc_map(lrows, i, j)
        irow = urows[i]
        for j in range (i+1,m):
            u_ji = udata_get((j,i))
            if u_ji is None:
                continue
            c = div(u_ji, a_ii)
            jrow = urows[j]
            jrow_add = jrow.add
            jrow_remove = jrow.remove
            for p in irow:
                if p < i: continue
                u_ip_c = udata[i,p] * c
                jp = j,p
                b = udata_get(jp)
                if b is None:
                    udata[jp] = -u_ip_c
                    jrow_add(p)
                else:
                    if b==u_ip_c:
                        del udata[jp]
                        jrow_remove(p)
                    else:
                        udata[jp] = b - u_ip_c
            if not jrow:
                del urows[j]
            if i<k and j<m:
                ldata[j, i] = c
                s = lrows.get(j)
                if s is None:
                    s = lrows[j] = set()
                s.add(i)
    for i in xrange(min(m, k)):
        ldata[i,i] = 1
    if m>n:
        crop_MATRIX(k, n, udata)
    return pivot_table

def lu_MATRIX_T(m, n, k, ldata, udata):
    if m>n:
        for i in xrange(n,m):
            udata[(i,i)] = 1    
    pivot_table = range(m)
    udata_get = udata.get
    udata_has = udata.has_key
    urows = get_rc_map_T(udata)
    lrows = get_rc_map(ldata)
    for i in xrange(m-1):
        ncols, j = n, None
        for j1 in range(i, m):
            if udata_has((i,j1)):
                l = len(urows[j1])
                if l <= ncols:
                    ncols, j = l, j1
        if j is None:
            continue
        a_ii = udata[i,j]
        if i!=j:
            pivot_table[i], pivot_table[j] = pivot_table[j], pivot_table[i]
            swap_rows_MATRIX_T(udata, i, j, urows)
            swap_rows_MATRIX(ldata, i, j, lrows)
            swap_rc_map(urows, i, j)
            swap_rc_map(lrows, i, j)
        irow = urows[i]
        for j in xrange(i+1,m):
            u_ji = udata_get((i,j))
            if u_ji is None:
                continue
            c = div(u_ji, a_ii)
            jrow = urows[j]
            jrow_add = jrow.add
            jrow_remove = jrow.remove
            for p in irow:
                if p < i: continue
                u_ip_c = udata[p,i] * c
                jp = p,j
                b = udata_get(jp)
                if b is None:
                    udata[jp] = -u_ip_c
                    jrow_add(p)
                else:
                    if b==u_ip_c:
                        del udata[jp]
                        jrow_remove(p)
                    else:
                        udata[jp] = b - u_ip_c
            if not jrow:
                del urows[j]
            if i<k and j<m:
                ldata[j, i] = c
                s = lrows.get(j)
                if s is None:
                    s = lrows[j] = set()
                s.add(i)
    for i in xrange(min(m, k)):
        ldata[i,i] = 1
    if m>n:
        crop_MATRIX_T(k, n, udata)
    return pivot_table

def MATRIX_DICT_crop(self):
    """ Remove matrix elements that are out of dimensions inplace and return the matrix.
    """
    if not self.is_writable:
        raise TypeError('Cannot crop read-only matrix inplace')
    head, data = self.pair
    m, n = head.shape
    if head.is_transpose:
        crop_MATRIX_T(m, n, data)
    elif head.is_diagonal:
        raise NotImplementedError(`head`)
    else:
        crop_MATRIX(m, n, data)
    return self

def crop_MATRIX(m, n, data):
    for (i,j) in data.keys():
        if 0<=i<m and 0<=j<n:
            continue
        del data[i,j]

def crop_MATRIX_T(m, n, data):
    for (j, i) in data.keys():
        if 0<=i<m and 0<=j<n:
            continue
        del data[j, i]

def swap_rows_MATRIX(data, i, j, rows=None):
    if i==j:
        return
    if rows is not None:
        d = {}
        data_pop = data.pop
        for k in rows.get(i,[]):
            d[j,k] = data_pop((i,k))
        for k in rows.get(j,[]):
            d[i,k] = data_pop((j,k))
        data.update(d)
        return

    row_i = []
    row_j = []
    for index, element in data.items():
        i0 = index[0]
        if i0==i:
            row_i.append((index[1], element))
            del data[index]
        elif i0==j:
            row_j.append((index[1], element))
            del data[index]
    for k, element in row_i:
        data[j, k] = element
    for k, element in row_j:
        data[i, k] = element

def swap_cols_MATRIX(data, i, j, cols=None):
    if i==j:
        return
    if cols is not None:
        d = {}
        data_pop = data.pop
        for k in cols.get(i,[]):
            d[k,j] = data_pop((k,i))
        for k in cols.get(j,[]):
            d[k,i] = data_pop((k,j))
        data.update(d)
        return
    column_i = []
    column_j = []
    for index, element in data.items():
        i0 = index[1]
        if i0==i:
            column_i.append((index[0], element))
            del data[index]
        elif i0==j:
            column_j.append((index[0], element))
            del data[index]
    for k, element in column_i:
        data[k, j] = element
    for k,element in column_j:
        data[k, i] = element

swap_rows_MATRIX_T = swap_cols_MATRIX
swap_cols_MATRIX_T = swap_rows_MATRIX


def MATRIX_DICT_swap_rows(self, i, j):
    if not self.is_writable:
        raise TypeError('Cannot swap rows of a read-only matrix')
    head, data = self.pair
    if head.is_transpose:
        swap_rows_MATRIX_T(data, i, j)
    elif head.is_diagonal:
        raise NotImplementedError(`head`)
    else:
        swap_rows_MATRIX(data, i, j)

def MATRIX_DICT_swap_cols(self, i, j):
    if not self.is_writable:
        raise TypeError('Cannot swap columns of a read-only matrix')
    head, data = self.pair
    if head.is_transpose:
        swap_cols_MATRIX_T(data, i, j)
    elif head.is_diagonal:
        raise NotImplementedError(`head`)
    else:
        swap_cols_MATRIX(data, i, j)

def MATRIX_DICT_trace(self):
    """ Return trace of a matrix.
    """
    head, data = self.pair
    m, n = head.shape
    s = 0
    if head.is_diagonal:
        dget = data.get
        for i in xrange(min(n,m)):
            s += dget((i, i), 0)
        return s
    if m != n:
        raise ValueError("matrix trace is only defined for square matrices but got %sx%s" % (m,n))
    sparse = len(data) < m
    if sparse:
        for (i, j), element in data.items():
            if i == j:
                s += element
    else:
        dget = data.get
        for i in xrange(n):
            s += dget((i, i), 0)
    return s

