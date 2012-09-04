"""Low-level arithmetics support.

Arithmetic package provides:

  * low-level number types ``mpq``, ``mpf``, ``mpqc``, ``mpc``
  * base class ``Infinity`` for extended numbers
  * various low-level number theory functions: ``gcd``, ``multinomial_coefficients``, etc.

"""
__docformat__ = "restructuredtext"

from . import mpmath
from .numbers import mpq, mpf, mpqc, mpc, setdps, getdps
from .number_theory import gcd, lcm, multinomial_coefficients, f2q
from .infinity import Infinity
