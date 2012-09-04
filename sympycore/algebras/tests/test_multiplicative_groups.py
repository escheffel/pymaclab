
from __future__ import with_statement

from sympycore import *
import sympycore.algebras.groups as groups

def test_additive_group_operations():
    Group = groups.MultiplicativeGroup
    x,y,z = map(Group, 'xyz')

    operands = [x, y, x*y, x*x, Group(1), 1, 'x', 'x*y']
    unary_operations = []
    binary_operations = ['*', '/']

    utils.test_operations(operands, expected_results_multiplicative_group+expected_results_numbers, unary_operations, binary_operations)

    operands = ([x, Group(1)], [1,2, -1, -2, 0])
    binary_operations = ['**']
    utils.test_operations(operands, expected_results_multiplicative_group+expected_results_numbers, unary_operations, binary_operations)

expected_results_multiplicative_group = '''

(x)/(x):1
(x)/(y):x/y
(x)/(x*y):x/y/x
(x)/(x**2):1/x
(x)/(1):x
(y)/(x):y/x
(y)/(y):1
(y)/(x*y):1/x
(y)/(x**2):y/x**2
(y)/(1):y
(x*y)/(x):x*y/x
(x*y)/(y):x
(x*y)/(x*y):1
(x*y)/(x**2):x*y/x**2
(x*y)/(1):x*y
(x**2)/(x):x
(x**2)/(y):x**2/y
(x**2)/(x*y):x**2/y/x
(x**2)/(x**2):1
(x**2)/(1):x**2
(1)/(x):1/x
(1)/(y):1/y
(1)/(x*y):1/y/x
(1)/(x**2):1/x**2


(x)**(1):x
(x)**(2):x**2
(x)**(-1):1/x
(x)**(-2):1/x**2
(x)**(0):1

(x)*(x):x**2
(x)*(y):x*y
(y)*(x):y*x
(y)*(y):y**2
(x)*(x*y):x**2*y
(y)*(x*y):y*x*y
(x*y)*(x):x*y*x
(x*y)*(y):x*y**2
(x*y)*(x*y):x*y*x*y
(x)*(x**2):x**3
(y)*(x**2):y*x**2
(x*y)*(x**2):x*y*x**2

(x**2)*(x):x**3
(x**2)*(y):x**2*y
(x**2)*(x*y):x**3*y
(x**2)*(x**2):x**4

(x)*(1):x
(y)*(1):y
(x*y)*(1):x*y
(x**2)*(1):x**2
(1)*(x):x
(1)*(y):y
(1)*(x*y):x*y
(1)*(x**2):x**2

'''

expected_results_numbers = '''
(1)*(1):1

(1)**(1):1
(1)**(2):1
(1)**(-1):1
(1)**(-2):1
(1)**(0):1
(1)/(1):1
'''
