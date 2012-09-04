
#from __future__ import division

import copy
from sympycore import Logic, SymbolicEquality, Calculus, heads, Symbol
from sympycore.arithmetic.numbers import div
from . import Matrix

class Polyhedron(object):
    """ A basis class that provides a H-representation for general
    polyhedra.

    Consider the following polyhedron::

      P(A_I, b_I, A_L, b_L) = { x in R^n | A_I*x <= b_I, A_L*x == b_L }

    where ``A_I`` and ``A_L`` are matrices with ``|I|`` and ``|L|`` rows,
    respectively, and ``I``, ``L`` are partitions of ``{0,..,m-1}`` such that
    ``set(I + L)==set(range(m))``, and ``*`` denotes matrix dot product.

    The class ``Polyhedron`` stores inequalities and equalities in a
    homogenized form, that is, in a m times (n+1) matrix A where the
    first column corresponds to coefficients in ``-b_I`` and ``-b_L``,
    and the rest of the columns correspond to coefficents in ``A_I``
    and ``A_L``. With respect to coefficent storage, this class
    represents the following homogeneous polyhedron::

      P(A) = { x' in R^{n+1} | A[I] * x' <= 0, A[L] * x' == 0}
    
    where ``x'`` may contain special variables, see ``labels``
    attributes for details.

    Attributes
    ----------
    A : Matrix
      A constraints matrix.

    I, L : list
      Lists of row indices of ``A``. If ``i`` in ``I`` then the
      ``i``-th row corresponds to inequality relation and if ``i`` in
      ``L`` then to equality relation. See above for definitions.

    labels : list
      List of column labels of matrix A.  Some labels have special
      meaning, as defined below:

        ``"1"`` --- the name corresponds to homogenization variable
          and the corresponding column in A corresponds to coefficents
          in ``-b_I`` and ``-b_L``. To relate P(A) to P(A_I, b_I, A_L,
          b_L), ``"1"`` value should be taken to 1.

        ``"SLACK#"`` --- the name corresponds to slack variable. The
          slack variable is used in a LP problem where it is assumed
          to be non-negative. Here ``"#"`` denotes integers used to
          index slack variables corresponding to inequality rows.

        ``"MAX"`` or ``"MIN"`` --- the labels correspond to LP
          objective function values when used in an equality relation.
          Minimization LP problem is automatically converted to
          maximization LP problem.

    objective_row, objective_column : {None, int}
      Row and column indices of the objective coefficients and ``"MAX"``
      column name, respectively.

    rhs_column : {None, int}
      Column index corresponding to RHS coefficients.
    """

    def __init__(self, *exprs):
        self.A = Matrix(0,0)
        self.I = []
        self.L = []
        self.labels = []
        self.objective_row = None
        self.objective_column = None
        self.rhs_column = None
        for expr in exprs:
            self.add(expr)

    def copy(self):
        cpy = Polyhedron()
        cpy.A = self.A.copy()
        for n in ['I', 'L', 'labels']:
            getattr(cpy, n).extend(getattr(self, n))
        for n in ['objective_row', 'objective_column']:
            setattr(cpy, n, getattr(self, n))
        return cpy
    
    @property
    def A_L(self):
        """ Return a matrix of equalities.
        """
        return self.A[tuple(self.L),:]

    @property
    def A_I(self):
        """ Return a matrix of inequalities.
        """
        return self.A[tuple(self.I),:]    

    def add(self, expr):
        """ Add extra constraint to polyhedron.

        Patameters
        ----------
        expr : {str, Logic}

          Symbolic expression of a constraint. The expression must be
          an equality or inequality relation of arbitrary linear
          combination of variables.
        """
        if isinstance(expr, str):
            if not expr or expr.strip ().startswith ('#'):
                return
            if '=' in expr or '<' in expr or '>' in expr:
                expr = Logic(expr)
            else:
                expr = Logic(expr+'==0')
        elif isinstance (expr, Calculus):
            expr = Logic (heads.EQ, (expr, Calculus(0)))
        if not isinstance(expr, Logic):
            raise NotImplementedError (`expr`) # expected Logic instance

        op, (lhs, rhs) = expr.pair
        i = self.A.rows

        if op==heads.EQ:
            expr = lhs - rhs
            self.L.append(i)
        elif op in [heads.LT, heads.LE]:
            op = heads.LE
            expr = lhs - rhs
            self.I.append(i)
        elif op in [heads.GT, heads.GE]:
            op = heads.LE
            expr = rhs - lhs
            self.I.append(i)
        else:
            raise NotImplementedError(`op`) # expected ==, <, <=, >, >=

        expr = expr.to(heads.TERM_COEFF_DICT)

        negate_row = False
        if expr.head is heads.NUMBER:
            try:
                j = self.labels.index('1')
            except ValueError:
                j = len(self.labels)
                self.labels.append(Symbol('1'))
            if self.rhs_column is None:
                self.rhs_column = j
            else:
                assert self.rhs_column == j,`self.rhs_column,j`
            self.A[i,j] = expr.data
        elif expr.head is heads.TERM_COEFF:
            term, coeff = expr.data
            if term=='MIN':
                negate_row = True
                coeff = -coeff
                term = Symbol('MAX')
            try:
                j = self.labels.index(term)
            except ValueError:
                j = len(self.labels)
                self.labels.append(term)
            if term=='MAX':
                assert self.objective_row is None, `self.objective_row, term`
                self.objective_row = i
                self.objective_column = j
            self.A[i,j] = coeff
        elif expr.head is heads.TERM_COEFF_DICT:
            for term, coeff in expr.data.iteritems():
                if term.head == heads.NUMBER:
                    assert term.data==1,`term.pair` # expected Calculus('1')
                    try:
                        j = self.labels.index('1')
                    except ValueError:
                        j = len(self.labels)
                        self.labels.append(Symbol('1'))
                    if self.rhs_column is None:
                        self.rhs_column = j
                    else:
                        assert self.rhs_column == j,`self.rhs_column,j`
                    self.A[i,j] = coeff
                else:
                    if term=='MIN':
                        negate_row = True
                        coeff = -coeff
                        term = Symbol('MAX')
                    try:
                        j = self.labels.index(term)
                    except ValueError:
                        j = len(self.labels)
                        self.labels.append(term)
                    if term=='MAX':
                        assert self.objective_row is None, `self.objective_row, term`
                        self.objective_row = i
                        self.objective_column = j
                    self.A[i,j] = coeff
        else:
            try:
                j = self.labels.index(expr)
            except ValueError:
                j = len(self.labels)
                self.labels.append(expr)
            self.A[i,j] = 1            

        self.A = self.A.resize(i+1, len (self.labels))

    def get_LP(self):
        """Construct LP dictionary matrix from polyhedra constraints.

        Parameters
        ----------
        None

        Returns
        -------
        variables : list
          List of variables.
        D : Matrix
          LP dictionary matrix.

        See also
        --------
        MatrixDict.LP_solve
        """
        assert len (self.L)==1,`self.L, self.objective_row`
        colindices = []

        variables = []
        for name, i in sorted (zip(map (str, self.labels), range (len (self.labels)))):
            if name=='1':
                pass
            elif name=='MAX':
                pass
            else:
                colindices.append(i)
                variables.append(name)
        colindices = tuple(colindices)
        rowindices = tuple(self.I)                
        D = Matrix(self.A.rows, len (colindices)+1)
        D[1:,1:] = -self.A[rowindices,colindices]
        D[0,1:] = self.A[self.objective_row, colindices] / (-self.A[self.objective_row, self.objective_column])
        if self.rhs_column is not None:
            D[:,0] = -self.A[:,self.rhs_column]

        return variables, D

    def show(self):
        print 'A:',', '.join(map(str, self.labels))
        print self.A
        print 'L:',self.L
        print 'I:',self.I
