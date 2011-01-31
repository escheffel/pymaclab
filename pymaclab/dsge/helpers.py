import datetime
import time
import copy as COP
import os as OPS
import sys as SYST
import re as RE
import string as STR
import scipy as S
from numpy import *
from numpy import matlib as MAT
from scipy import linalg as LIN
from numpy.ctypeslib import load_library, ndpointer
from ctypes import cdll, c_int, c_char, POINTER


def null(A, eps = 1e-15):
	"""
	Function to get the nullspace of a matrix.

	Usage:
	B = nullspace(A,eps)

	where A is some matrix, B is the returned
	corresponding nullspace and eps is the optional
	machine precision (smallest number).
	"""
	u, s, vh = LIN.svd(A)
	null_mask = (s <= eps)
	null_space = compress(null_mask, vh, axis=0)
	return null_space

def eig(A,B):
	"""
	To ensure matlab compatibility, we need to
	swap matrices A and B around !!
	"""
	(XX1,XX2) = LIN.eig(A,B)
	(XX3,XX4) = LIN.eig(B,A)
	return (mat(XX4),mat(XX1))

def sortrows(A, integer):
	"""
	Sort the rows of a matrix using the
	column given by the index number
	"""
	A_array = A.A
	A_array = A_array[A_array[:,integer].argsort(),]
	return mat(A_array)

def ctr(A):
	"""
	returns conjugate transpose
	"""
	return A.conjugate().transpose()

def eig2(A,B, lapackname='', lapackpath=''):
	'''Calculates generalized eigenvalues of pair (A,B).

    This should correspond to Matlab's lambda = eig(A,B),
    and also to some (the same?) scipy function.

    Eigenvalues will be of complex type, are unsorted, and are returned as 1d.
    '''
	AA,BB,dum1,dum2,VSL,VSR = zgges4numpy(A,B,
													  lapackname=lapackname,lapackpath=lapackpath)
	return diag(AA)/diag(BB)

# Time Stamp Stuff
def now_is():
	t = datetime.datetime.now()
	EpochSeconds = time.mktime(t.timetuple())
	now = datetime.datetime.fromtimestamp(EpochSeconds)
	return now.__str__()

### EVERYTHING BELOW IS COMMENTED OUT ###
### FUNCTIONALITY HAS BEEN REPLACED ###

#NOTE: numpy now includes linalg.matrix_rank
#def rank(A, eps=1e-15):
#	"""
#	Function to calcute the mathematical rank of matrix.
#	This was coded because NUMPY's rank functions give the
#	matrix rank of a matrix AND NOT the rank in a linear
#	algebra sense of returning the max. amount of linearly
#	independent rows or columns in a matrix.
#
#	Usage:
#	rank = rankm(A,eps)
#
#	where A is some matrix, rank is then returned integer
#	number for rank, and eps is optional machine precision
#	(smallest number).
#	"""
#	u, s, vh = LIN.svd(A)
#	tol = max(A.shape)*eps*abs(max(s))
#	rank = MAT.sum(s > tol)
#	return rank


## Sven Schreiber's Python functions
#'''
#QZ alias generalized Schur decomposition (complex or real) for Python/Numpy.
#
#You need to import the qz() function of this module, check out its docstring,
#especially what it says about the required lapack shared library. Run this 
#module for some quick tests of the setup.
#
#This is free but copyrighted software, distributed under the same license
#as Python 2.5, copyright Sven Schreiber.
#
#If you think a different license would make (more) sense, please say so
#on the Numpy mailing list (see scipy.org).
#''' 
#
#def setuplapack4xgges(A,B,lpname,lppath):
#	'''Loads the lapack shared lib and does some input checks.
#
#    The defaults for lapackname and location are platform-specific:
#        Win32: 'lapack' (due to scilab's lapack.dll)
#               'c:\\winnt\\system32\\'
#        Otherwise: 'liblapack' 
#                   '/usr/lib/'
#    '''
#	# some input checks
#	assert A.ndim == 2
#	assert A.shape == B.shape
#	assert A.shape[0] == A.shape[1]
#	# load the lapack shared library
#	if lpname == '':
#		if SYST.platform == 'win32':    lpname = 'lapack'
#		else:                          lpname = 'liblapack'  
#	if lppath == '':
#		if SYST.platform == 'win32': lppath = 'c:\\winnt\\system32\\'
#		else:                       lppath = '/usr/lib/' 
#	lapack = load_library(lpname, lppath)
#	return lapack
#
#def dgges4numpy(A,B, jobvsl='V', jobvsr='V', lapackname='', lapackpath=''):
#	'''wraps lapack function dgges, no sorting done'''
#	lapack = setuplapack4xgges(A,B,lapackname,lapackpath)
#	rows = A.shape[0]
#	# to determine matrix subclass
#	Aintype = type(A)
#
#	# actual inputs
#	A = asfortranarray(A, dtype=float64)
#	B = asfortranarray(B, dtype=float64)
#	# seems ok to pass strings directly, but the function expects only 1 char!
#	jobvsl = jobvsl[0]
#	jobvsr = jobvsr[0]
#
#	# dummy inputs
#	sort    = 'N'            # we don't want sorting
#	dummy   = 0     # 
#	info    = c_int(1)
#	lda     = c_int(rows)
#	ldb     = c_int(rows)
#	ldvsl   = c_int(rows)
#	ldvsr   = c_int(rows)
#	plwork  = 16*rows        # needed again later
#	lwork   = c_int(plwork)
#	n       = c_int(rows)
#	csdim   = c_int(0)    # because we don't sort
#
#	# auxiliary arrays
#	Alphar = asfortranarray(empty(rows), dtype=float64)
#	Alphai = asfortranarray(empty(rows), dtype=float64)
#	Beta  = asfortranarray(empty(rows), dtype=float64)
#	Vsl   = asfortranarray(empty([rows,rows]), dtype=float64)
#	Vsr   = asfortranarray(empty([rows,rows]), dtype=float64)
#	Work  = asfortranarray(empty(plwork), dtype=float64)
#	Rwork = asfortranarray(empty(8*rows), dtype=float64)
#
#	lapack.dgges_.argtypes = [  
#		POINTER(c_char),                                       # JOBVSL
#		POINTER(c_char),                                       # JOBVSR
#		POINTER(c_char),                                       # SORT
#		# for the dummy the POINTER thing didn't work, 
#		#  but plain c_int apparently does...
#		c_int,                                          # dummy SELCTG 
#		POINTER(c_int),                                        # N
#		ndpointer(dtype=float64, ndim=2, flags='FORTRAN'),  # A
#		POINTER(c_int),                                        # LDA
#		ndpointer(dtype=float64, ndim=2, flags='FORTRAN'),  # B
#		POINTER(c_int),                                        # LDB
#		POINTER(c_int),                                        # SDIM
#		ndpointer(dtype=float64, ndim=1, flags='FORTRAN'),  # ALPHAr
#		ndpointer(dtype=float64, ndim=1, flags='FORTRAN'),  # ALPHAi
#		ndpointer(dtype=float64, ndim=1, flags='FORTRAN'),  # BETA
#		ndpointer(dtype=float64, ndim=2, flags='FORTRAN'),  # VSL
#		POINTER(c_int),                                        # LDVSL
#		ndpointer(dtype=float64, ndim=2, flags='FORTRAN'),  # VSR
#		POINTER(c_int),                                        # LDVSR
#		ndpointer(dtype=float64, ndim=1, flags='FORTRAN'),  # WORK
#		POINTER(c_int),                                              # LWORK
#		# same as with SELCTG...
#		c_int,                                                 # dummy BWORK 
#		POINTER(c_int) ]                                       # INFO
#
#	lapack.dgges_(jobvsl,jobvsr,sort,dummy,n,A,lda,B,ldb,csdim,Alphar,Alphai,
#					  Beta,Vsl,ldvsl,Vsr,ldvsr,Work,lwork,dummy,info)
#
#	# preserve matrix subclass
#	if Aintype == type(mat(1)):
#		A=mat(A); B=mat(B); Vsl=mat(Vsl); Vsr=mat(Vsr)
#	if info.value == 0:
#		if   jobvsl=='V' and jobvsr=='V': return A,B,Alphar,Alphai,Beta,Vsl,Vsr
#		elif jobvsl=='V' and jobvsr=='N': return A,B,Alphar,Alphai,Beta,Vsl
#		elif jobvsl=='N' and jobvsr=='V': return A,B,Alphar,Alphai,Beta,Vsr
#		else:                             return A,B,Alphar,Alphai,Beta
#	elif info.value < 0:
#		raise ValueError, 'Illegal argument (' + str(abs(info.value)) + ')'
#	elif info.value <= rows: 
#		raise RuntimeError, 'QZ iteration failed'
#	elif info.value <= rows+3:
#		raise RuntimeError, 'something other than QZ iteration failed'
#	else: raise RuntimeError, 'INFO not updated by dgges, complete failure!?'
#
#
#def zgges4numpy(A,B, jobvsl='V', jobvsr='V', lapackname='', lapackpath=''):
#	'''Wraps lapack function zgges, no sorting done.
#
#    Returns complex arrays, use real_if_close() if needed/possible.
#    '''
#	lapack = setuplapack4xgges(A,B,lapackname,lapackpath)
#	rows = A.shape[0]
#	# determine matrix subclass
#	Aintype = type(A)
#
#	# actual inputs
#	# The COMPLEX*16 type in Fortran translates to numpy's complex128
#	A = asfortranarray(A, dtype=complex128)
#	B = asfortranarray(B, dtype=complex128)
#	# seems ok to pass strings directly, but the function expects only 1 char!
#	jobvsl = jobvsl[0]
#	jobvsr = jobvsr[0]
#
#	# dummy inputs
#	sort = 'N'         # we don't want sorting
#	dummy = 0           # a placeholder for what would be needed for sorting 
#	info = c_int(rows+4)  # >n+3 aren't used as error codes of zgges
#	lda = c_int(rows)
#	ldb = c_int(rows)
#	ldvsl = c_int(rows)
#	ldvsr = c_int(rows)
#	plwork = 16*rows        # needed again later
#	lwork = c_int(plwork)
#	n = c_int(rows)
#	sdim = c_int(0)    # because we don't sort
#
#	# auxiliary arrays
#	Alpha = asfortranarray(empty(rows), dtype=complex128)
#	Beta  = asfortranarray(empty(rows), dtype=complex128)
#	Vsl   = asfortranarray(empty([rows,rows]), dtype=complex128)
#	Vsr   = asfortranarray(empty([rows,rows]), dtype=complex128)
#	Work  = asfortranarray(empty(plwork), dtype=complex128)
#	Rwork = asfortranarray(empty(8*rows), dtype=float64)
#
#	lapack.zgges_.argtypes = [  
#		POINTER(c_char),                                         # JOBVSL
#		POINTER(c_char),                                         # JOBVSR
#		POINTER(c_char),                                         # SORT
#		c_int,                                             # dummy SELCTG 
#		POINTER(c_int),                                          # N
#		ndpointer(dtype=complex128, ndim=2, flags='FORTRAN'), # A
#		POINTER(c_int),                                          # LDA
#		ndpointer(dtype=complex128, ndim=2, flags='FORTRAN'), # B
#		POINTER(c_int),                                          # LDB
#		POINTER(c_int),                                          # SDIM
#		ndpointer(dtype=complex128, ndim=1, flags='FORTRAN'), # ALPHA
#		ndpointer(dtype=complex128, ndim=1, flags='FORTRAN'), # BETA
#		ndpointer(dtype=complex128, ndim=2, flags='FORTRAN'), # VSL
#		POINTER(c_int),                                          # LDVSL
#		ndpointer(dtype=complex128, ndim=2, flags='FORTRAN'), # VSR
#		POINTER(c_int),                                          # LDVSR
#		ndpointer(dtype=complex128, ndim=1, flags='FORTRAN'), # WORK
#		POINTER(c_int),                                          # LWORK
#		ndpointer(dtype=float64, ndim=1, flags='FORTRAN'),    # RWORK
#		c_int,                                             # dummy BWORK 
#		POINTER(c_int) ]                                         # INFO
#
#	lapack.zgges_(jobvsl,jobvsr,sort,dummy,n,A,lda,B,ldb,sdim,Alpha,
#					  Beta,Vsl,ldvsl,Vsr,ldvsr,Work,lwork,Rwork,dummy,info)
#
#	# preserve matrix subclass
#	if Aintype == type(mat(1)):
#		A=mat(A); B=mat(B); Vsl=mat(Vsl); Vsr=mat(Vsr)
#	# use .value for ctypes safety, although probably redundant
#	if info.value == 0:
#		if   jobvsl=='V' and jobvsr=='V': return A,B,Alpha,Beta,Vsl,Vsr
#		elif jobvsl=='V' and jobvsr=='N': return A,B,Alpha,Beta,Vsl
#		elif jobvsl=='N' and jobvsr=='V': return A,B,Alpha,Beta,Vsr
#		else:                             return A,B,Alpha,Beta
#	elif info.value < 0:
#		raise ValueError, 'Illegal argument (' + str(abs(info.value)) + ')'
#	elif info.value <= rows: 
#		raise RuntimeError, 'QZ iteration failed'
#	elif info.value <= rows+3:
#		raise RuntimeError, 'something other than QZ iteration failed'
#	else: raise RuntimeError, 'INFO not updated by zgges, complete failure!?'
#
#def qz(A,B, mode='complex', lapackname='', lapackpath=''):
#	'''Equivalent to Matlab's qz function [AA,BB,Q,Z] = qz(A,B).
#
#    Requires Lapack as a shared compiled library on the system (one that
#    contains the functions dgges for real and zgges for complex use -- on 
#    Windows the one shipped with Scilab works). The underlying defaults for 
#    lapackname and lapackpath are platform-specific:
#        Win32: 'lapack' (due to scilab's lapack.dll)
#               'c:\\winnt\\system32\\'
#        Otherwise: 'liblapack' 
#                   '/usr/lib/'
#
#    This function should exactly match Matlab's usage, unlike octave's qz 
#    function which returns the conjugate-transpose of one of the matrices. Thus
#    it holds that 
#        AA = Q*A*Z
#        BB = Q*B*Z,
#    where Q and Z are unitary (orthogonal if real).
#
#    If mode is 'complex', then:
#     returns complex-type arrays, 
#     AA and BB are upper triangular, 
#     and diag(AA)/diag(BB) are the generalized eigenvalues of (A,B).
#
#    If the real qz decomposition is explicitly requested --as in Matlab:  
#    qz(A,B,'real')-- then:
#     returns real-type arrays,
#     AA is only block upper triangular,
#     and calculating the eigenvalues is more complicated.
#
#    Other variants like [AA,BB,Q,Z,V,W] = qz(A,B) are not implemented, i.e.
#    no generalized eigenvectors are calculated.
#    '''
#	if mode == 'real':
#		AA,BB,dum1,dum2,dum3,VSL,VSR = dgges4numpy(A,B,
#																 lapackname=lapackname,lapackpath=lapackpath)
#		return AA, BB, VSL.T, VSR
#	elif mode == 'complex':
#		AA,BB,dum1,dum2,VSL,VSR = zgges4numpy(A,B,
#														  lapackname=lapackname,lapackpath=lapackpath)
#		return AA, BB, VSL.conj().T, VSR
#	else: raise ValueError, 'bogus choice for mode'


#def eigwithqz(A,B, lapackname='', lapackpath=''):
#	'''Does complex QZ decomp. and also returns the eigenvalues'''
#	AA, BB, Q, Z = qz(A,B,lapackname=lapackname,lapackpath=lapackpath)
#	evals = diag(AA)/diag(BB)
#	return evals,AA,BB,Q,Z
