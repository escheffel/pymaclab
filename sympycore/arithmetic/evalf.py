"""Provides low-level evalf function.
"""

__all__ = ['evalf']

import math
import cmath
import re
import mpmath

int_pattern = re.compile('\d+([.]\d+)?')

def quote_numbers(s, format):
    return int_pattern.sub(lambda m: (format % m.group()), s)

def replace_names(s, names):
    for name, new in names.items():
        s = s.replace(name, new)
    return s

def f_header(symbols):
    if isinstance(symbols, (tuple, list)):
        return "lambda %s: " % ",".join(map(str, symbols))
    if not symbols:
        return "lambda: "
    return "lambda %s: " % symbols

def convert_mpmath(expr):
    s = str(expr)
    s = quote_numbers(s, "mpf(%s)")
    s = replace_names(s, { 'I':'j', 'E':'e', 'oo':'inf','undefined':'nan',
                           'Sin':'sin','Cos':'cos', 'Exp':'exp'})
    return s

def compile_mpmath(symbols, expr):
    s = convert_mpmath(expr)
    f = eval(f_header(symbols) + s, vars(mpmath))
    return f

def evalf(expr, digits=15):
    s = convert_mpmath(expr)
    return eval(s, vars(mpmath))
