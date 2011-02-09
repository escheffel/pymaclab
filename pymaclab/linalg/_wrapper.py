"""
Python drivers for LAPACK's dgges and dtgsen.
"""
from __future__ import with_statement
from qz import dgges, zgges
from ordqz import dtgsen#, ztgsen
from numpy import zeros, asanyarray, iscomplex, diag, errstate, matrix
from numpy import any as npany
from numpy.linalg import LinAlgError

def qz(A, B, mode='complex'):
    """
    QZ decompostion for generalized eigenvalues of a pair of matrices.

    Parameters
    ----------
    A : array-like
        
    B : array-like
    
    mode : str {'real','complex'}

    Returns
    -------
    AA, BB, Q.T, and Z

    Notes
    -----
    If either A or B is a matrix, matrices are returned
    
    """
    A = asanyarray(A)
    B = asanyarray(B)
    if isinstance(A,matrix) or isinstance(B, matrix):
        usematrix = True
    else:
        usematrix = False
    #NOTE: the first arguments to both are the selector functions but
    #sort = 'N' is hard coded in qz.pyf so both are ignored, but the
    #arguments are checked
    n = A.shape[1]
    if mode == 'real':
#        sort = lambda x,y,z: (x+y)/z <= 1.
        (AA,BB,sdim,alphar,alphai,
            beta, Q, Z, work, info) = dgges(lambda x,y,z: None,A ,B)
    if mode == 'complex':
        rwork = zeros(8*n)
        (AA,BB,sdim,alpha,beta,
            Q, Z, work, info) = zgges(lambda x,y: None,A ,B,rwork)
        # note that it will mostly be given real arrays, so it should
        # return real arrays, check that this is the case and return the reals
        if not npany(iscomplex(AA)):
            AA = AA.real
        if not npany(iscomplex(BB)):
            BB = BB.real
        if not npany(iscomplex(Q)):
            Q = Q.real
        if not npany(iscomplex(Z)):
            Z = Z.real
    if info < 0:
        raise ValueError("Incorrect value at pos. "+str(abs(i))+"in _gges")
    elif info > 0 and info <= n:
        raise LinAlgError("QZ iteration failed, but ALPHA(j) and BETA(j) \
should be correct for j = %s,...,%s" % (str(i),str(n))) # not zero-based
    elif info > n: # shouldn't get to this because sort='N'
        raise LinAlgError("Something other than QZ iteration failed. \
Return INFO = %s.  See _GGES.f for more information." % str(info))
    #return Q.T to follow matlab convention
    if usematrix:
        AA = matrix(AA)
        BB = matrix(BB)
        Q = matrix(Q)
        Z = matrix(Z)
    return AA,BB,Q.T,Z

def eigwithqz(A,B,mode='complex'):
    """
    Does QZ decomposition and also returns generalized eigenvalues.

    Returns
    -------
    eigvals, AA, BB, Q, Z

    Notes
    -----
    Just the eigenvalues can be recovered by scipy.linalg.eig(A,B)
    """
    AA, BB, Q, Z = qz(A,B,mode)
    with errstate(all='ignore'):
        eigvals = diag(AA)/diag(BB)
    return eigvals, AA, BB, Q, Z

def ordqz(AA,BB,Q,Z,select):
    """
    
    Parameters
    ----------
    AA
    BB
    Q
    Z
    select : array
        Boolean array that selects the eigenvalues.

    Notes
    -----
    Only the real valued dtgsen is included.
    """
    #NOTE: transpose Q to keep matlab convention
    (AAS,BBS,alphar,alphai,beta,QS,ZS,
        m,pl,pr,dif,info) = dtgsen(0,1,1,select,AA,BB,Q.T,Z)
    if info < 0:
        raise ValueError("Incorrect value at pos. "+str(abs(i))+"in dtgsen")
    elif info == 1:
        raise LinAlgError("Reordering failed because the transformed pair \
would be too far from generalized Schur form; the problem is very ill-conditioned")
    return AAS, BBS, QS.T, ZS

