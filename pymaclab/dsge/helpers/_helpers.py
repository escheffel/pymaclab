'''
.. module:: helpers
   :platform: Linux
   :synopsis: A module containing just a small number of useful functions, such as getting the current data and some special
              linear algebra routines compiled in externally.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''

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
