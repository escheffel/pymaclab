# Author: Pearu Peterson
# Created: January 2012

import copy
from sympycore.arithmetic.numbers import div
from .linalg import get_rc_maps
from .algebra import Matrix

class LPError(Exception):
    """Generic Python-exception-derived object raised by sympycore.linalg LP related functions.
    """
    pass

def MATRIX_DICT_LP_solve(self, method='crisscross', overwrite=False):
    """ Solve LP problem ``max(c^T x: A*x <= b & x>=0)``.

    Given matrix must contain the problem data in dictionary form::

          [ 0  c^T ]
      D = [        ]
          [ b  -A  ]

    Parameters
    ----------
    method : {'crisscross'}
      Specify the method to be used for solving the LP problem.

    overwrite : bool
      When True then discard the content of given matrix. After
      computation the dictionary matrix will be in either feasible or
      inconsistent form.
      
    Returns
    -------
    xopt : Matrix
      Basic solution of the LP problem when the LP problem is feasible.
    vopt : number
      Optimal value of the LP problem: ``c^T xopt``


    When the LP problem is inconsistent, LPError exception will be
    returned.
    """
    head, data = self.pair
    m, n = head.shape
    if not (overwrite and self.is_writable):
        data = dict(data)
    if method=='crisscross':
        if head.is_transpose or head.is_diagonal:
            raise NotImplementedError('LP_solve_crisscross_MATRIX_T, LP_solve_crisscross_MATRIX_D')
        return LP_solve_crisscross_MATRIX(m, n, data)
    else:
        raise NotImplementedError(`method`)

def LP_solve_crisscross_MATRIX(m, n, data):
    """ Solve LP problem using crisscross method. See
    MATRIX_DICT_LP_solve for definition and user interface.

    The implementation of this function is based on the crisscross
    algorithm given in:

      http://www.ifor.math.ethz.ch/teaching/Courses/Fall_2011/intro_fall_11/Script2

    The algorithm has few modifications to achieve more optimal Python
    implementation. For instance, the code uses native indexing of the
    dictionary matrix instead of label indexing. Also, finding the
    minimum of labels is replaced with finding the maximum of indexes.

    The current implementation is compact, that is, all pivot
    operations are applied in-situ on an input matrix taking full
    advantage of possible sparsity of the input matrix.
    """
    N = range(1,n) # variables are indices
    B = range(-1,-m,-1) # slack variables are negative indices

    while True:
        Drows, Dcols = get_rc_maps(data)

        def rowmax(j, Dcols, B):
            k, index = None, None
            for i in Dcols.get(j,[]):
                if i and data[i, j] < 0:
                    label = B[i-1]
                    if label > k:
                        k = label
                        index = i
            return k, index
        def colmax(i, Drows, N):
            k, index = None, None
            for j in Drows.get(i,[]):
                if j and data[i, j] > 0:
                    label = N[j-1]
                    if label > k:
                        k = label
                        index = j
            return k, index

        k, kindex = max(rowmax(0, Dcols, B), colmax(0, Drows, N))
        if k is None:
            break

        if k in B:
            r, Br = k, kindex
            s, Ns = colmax (Br, Drows, N)
            if s is None:
                raise LPError("LP problem is inconsistent")

        else: # k in N
            s, Ns = k, kindex
            r, Br = rowmax(Ns, Dcols, B)
            if r is None:
                raise LPError("LP problem is dual inconsistent")

        B[Br-1] = s
        N[Ns-1] = r
        id_rs = div(1, data[Br, Ns])
        negid_rs = -id_rs
        newrows = copy.deepcopy(Drows)

        srows = Dcols[Ns]
        srows.discard(Br)
        rcols = Drows[Br]
        rcols.discard(Ns)

        for i in srows:
            cols = newrows[i]
            cols.update(rcols)
            cols.discard(Ns)

        for Bi in srows:
            for Nj in newrows[Bi]:
                d_rj = data.get((Br, Nj))
                if d_rj is not None:
                    d_ij = data.get((Bi,Nj), 0)
                    d_ij -= d_rj*data[Bi,Ns]*id_rs
                    if d_ij==0:
                        del data[Bi,Nj]
                    else:
                        data[Bi,Nj] = d_ij

        for Bi in srows:
            data[Bi,Ns] *= id_rs

        for Nj in rcols:
            data[Br,Nj] *= negid_rs

        data[Br,Ns] = id_rs

    optimal = Matrix(n-1,1)
    for i,k in enumerate(B):
        if k>0:
            optimal[k-1,0] = data.get((i+1,0), 0)
    optimal_value = data.get ((0,0),0)
    return optimal, optimal_value

