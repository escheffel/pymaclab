
from __future__ import with_statement

from sympycore import *
import sympycore.algebras.groups as groups

def test_additive_abelian_group_inplace_operations():
    Group = groups.AdditiveAbelianGroup
    x,y,z = map(Group, 'xyz')
    s = s2 = x+y
    s_data_id = id(s.data)
    s += z
    assert s == x+y+z
    assert s_data_id==id(s.data)
    assert s == s2
    assert s is s2

    s *= 2
    assert s==2*x+2*y+2*z
    assert s_data_id==id(s.data)
    assert s == s2
    assert s is s2

    s -= z
    assert s==2*x+2*y+z
    assert s_data_id==id(s.data)
    assert s == s2
    assert s is s2

    s += x-y
    assert s==3*x+y+z
    assert s_data_id==id(s.data)
    assert s == s2
    assert s is s2

    s += 2
    assert s==3*x+y+z+2
    assert s_data_id==id(s.data)
    assert s == s2
    assert s is s2

def test_additive_abelian_group_operations():
    Group = groups.AdditiveAbelianGroup
    x,y,z = map(Group, 'xyz')

    operands = [x, y, -z, x+y, x+x, 0, Group(2), Group(1), 2, x-x, z-x, 'y + x', 'x']
    unary_operations = ['+', '-']
    binary_operations = ['+','-', '*']
    with core.UnevaluatedAddition(Group):
        operands.extend([x-x, x+x, 2*x, 2, Group(2), 0, x-z])

    utils.test_operations(operands, expected_results_additive_abelian_group+expected_results_numbers, unary_operations, binary_operations)    

def test_additive_abelian_group_unevaluated_operations():
    Group = groups.AdditiveAbelianGroup
    x,y,z = map(Group, 'xyz')
    operands = [x, y, x+x, x+y, x-x, x-z, '2*x', '0', 2]
    unary_operations = ['+', '-']
    binary_operations = ['+','-', '*']
    with core.UnevaluatedAddition(Group):
        utils.test_operations(operands, expected_results_additive_abelian_group_unevaluated+expected_results_numbers_unevaluated,
                              unary_operations, binary_operations)    

def test_additive_group_inplace_operations():
    Group = groups.AdditiveGroup
    x,y,z = map(Group, 'xyz')
    s = s2 = x + y
    s_data_id = id(s.data)
    s += z
    assert s == x+y+z
    assert s_data_id==id(s.data)
    assert s == s2
    assert s is s2
    
    s *= 2
    assert s == x+y+z+x+y+z
    assert s is not s2 # because creating new list and discarding old
                       # is more efficient that creating a new list
                       # and then copying it into the old list

    s_data_id = id(s.data)
    s2 = s
    s -= z
    assert s == x+y+z+x+y
    assert s_data_id==id(s.data)
    assert s is s2

    s -= x + y
    assert s==x+y+z
    assert s_data_id==id(s.data)
    assert s is s2
    
def test_additive_group_unevaluated_operations():
    Group = groups.AdditiveGroup
    x,y,z = map(Group, 'xyz')
    operands = [x, y, -z, x+y, x+x+x, x-x]
    unary_operations = ['+', '-']
    binary_operations = ['+','-', '*']
    with core.UnevaluatedAddition(Group):
        utils.test_operations(operands, expected_unevaluated_results_additive_group, unary_operations, binary_operations)

def test_additive_group_operations():
    Group = groups.AdditiveGroup
    x,y,z = map(Group, 'xyz')

    operands = [x, y, -z, x+x, x+y, y+x, x-x, 0, 2, Group(1), Group(2), 'x','x + y']
    unary_operations = ['+', '-']
    binary_operations = ['+','-', '*']

    utils.test_operations(operands, expected_results_additive_group+expected_results_numbers, unary_operations, binary_operations)
    operands = [x+y+x, 3, -4, x+y+z, x-y]
    utils.test_operations(operands, expected_results_additive_group+expected_results_numbers, unary_operations, binary_operations)

expected_results_additive_abelian_group_unevaluated = '''
+(x):x
-(x):-x
(x)+(x):x + x
(x)-(x):x + -x
(x)*(x):unsupported

(x)+(y):x + y
(x)-(y):x + -y
(x)*(y):unsupported
+(y):y
-(y):-y
(y)+(x):y + x
(y)-(x):y + -x
(y)*(x):unsupported
(y)+(y):y + y
(y)-(y):y + -y
(y)*(y):unsupported

(x)+(2*x):x + 2*x
(x)-(2*x):x + -(2*x)
(x)*(2*x):unsupported
(y)+(2*x):y + 2*x
(y)-(2*x):y + -(2*x)
(y)*(2*x):unsupported
+(2*x):2*x
-(2*x):-(2*x)

(2*x)+(x):2*x + x
(2*x)-(x):2*x + -x
(2*x)*(x):unsupported
(2*x)+(y):2*x + y
(2*x)-(y):2*x + -y
(2*x)*(y):unsupported
(2*x)+(2*x):2*x + 2*x
(2*x)-(2*x):2*x + -(2*x)
(2*x)*(2*x):unsupported

(x)+(y + x):x + y + x
(x)*(y + x):unsupported
(y)+(y + x):y + y + x
(y)*(y + x):unsupported
(2*x)+(y + x):2*x + y + x
(2*x)*(y + x):unsupported
+(y + x):y + x

(x)-(y + x):x + -(y + x)
(y)-(y + x):y + -(y + x)
(2*x)-(y + x):2*x + -(y + x)
-(y + x):-(y + x)
(y + x)+(x):y + x + x
(y + x)-(x):y + x + -x
(y + x)*(x):unsupported
(y + x)+(y):y + x + y
(y + x)-(y):y + x + -y
(y + x)*(y):unsupported
(y + x)+(2*x):y + x + 2*x
(y + x)-(2*x):y + x + -(2*x); y + x + -(x*2)
(y + x)*(2*x):unsupported
(y + x)+(y + x):y + x + y + x
(y + x)-(y + x):y + x + -(y + x)
(y + x)*(y + x):unsupported

(x)+(0):x + 0
(x)-(0):x + -0
(y)+(0):y + 0
(y)-(0):y + -0
(2*x)+(0):2*x + 0
(2*x)-(0):2*x + -0
(y + x)+(0):y + x + 0
(y + x)-(0):y + x + -0
(x)*(0):x*0
(y)*(0):y*0
(2*x)*(0):2*x*0
(0)*(x):0*x
(0)*(y):0*y
(0)*(2*x):0*2*x
(0)*(y + x):0*(y + x)
(y + x)*(0):(y + x)*0
(0)+(x):0 + x
(0)-(x):0 + -x
(0)+(y):0 + y
(0)-(y):0 + -y
(0)+(2*x):0 + 2*x
(0)-(2*x):0 + -(2*x)
(0)+(y + x):0 + y + x
(0)-(y + x):0 + -(y + x)

(x)+(x + -z):x + x + -z
(x)-(x + -z):x + -(x + -z)
(x)*(x + -z):unsupported
(y)+(x + -z):y + x + -z
(y)-(x + -z):y + -(x + -z)
(y)*(x + -z):unsupported
(2*x)+(x + -z):2*x + x + -z
(2*x)-(x + -z):2*x + -(x + -z)
(2*x)*(x + -z):unsupported
(y + x)+(x + -z):y + x + x + -z
(y + x)-(x + -z):y + x + -(x + -z)
(y + x)*(x + -z):unsupported
(0)+(x + -z):0 + x + -z
(0)-(x + -z):0 + -(x + -z)
(0)*(x + -z):0*(x + -z)
+(x + -z):x + -z
-(x + -z):-(x + -z)
(x + -z)+(x):x + -z + x
(x + -z)-(x):x + -z + -x
(x + -z)*(x):unsupported
(x + -z)+(y):x + -z + y
(x + -z)-(y):x + -z + -y
(x + -z)*(y):unsupported
(x + -z)+(2*x):x + -z + 2*x
(x + -z)-(2*x):x + -z + -(2*x)
(x + -z)*(2*x):unsupported
(x + -z)+(y + x):x + -z + y + x
(x + -z)-(y + x):x + -z + -(y + x)
(x + -z)*(y + x):unsupported
(x + -z)+(0):x + -z + 0
(x + -z)-(0):x + -z + -0
(x + -z)*(0):(x + -z)*0
(x + -z)+(x + -z):x + -z + x + -z
(x + -z)-(x + -z):x + -z + -(x + -z)
(x + -z)*(x + -z):unsupported

(x)+(2):x + 2
(x)-(2):x + -2
(x)*(2):x*2
(y)+(2):y + 2
(y)-(2):y + -2
(y)*(2):y*2
(2*x)+(2):2*x + 2
(2*x)-(2):2*x + -2
(2*x)*(2):2*x*2
(y + x)+(2):y + x + 2
(y + x)-(2):y + x + -2
(y + x)*(2):(y + x)*2
(x + -z)+(2):x + -z + 2
(x + -z)-(2):x + -z + -2
(x + -z)*(2):(x + -z)*2
(2)+(x):2 + x
(2)-(x):2 + -x
(2)*(x):2*x
(2)+(y):2 + y
(2)-(y):2 + -y
(2)*(y):2*y
(2)+(2*x):2 + 2*x
(2)-(2*x):2 + -(2*x)
(2)*(2*x):2*2*x
(2)+(y + x):2 + y + x
(2)-(y + x):2 + -(y + x)
(2)*(y + x):2*(y + x)

(2)+(x + -z):2 + x + -z
(2)-(x + -z):2 + -(x + -z)
(2)*(x + -z):2*(x + -z)

'''

expected_results_additive_abelian_group = '''
+(x):x
-(x):-x
(x)+(x):2*x
(x)-(x):0
(x)*(x):unsupported
(x)+(y):y + x
(x)-(y):-y + x
(x)*(y):unsupported
+(y):y
-(y):-y
(y)+(x):y + x
(y)-(x):y - x
(y)*(x):unsupported
(y)+(y):2*y
(y)-(y):0
(y)*(y):unsupported
(x)+(-z):x - z
(x)-(-z):x + z
(x)*(-z):unsupported
(y)+(-z):y - z
(y)-(-z):y + z
(y)*(-z):unsupported
+(-z):-z
-(-z):z
(-z)+(x):x - z

(-z)-(x):-x - z
(-z)*(x):unsupported
(-z)+(y):y - z
(-z)-(y):-y - z
(-z)*(y):unsupported
(-z)+(-z):-2*z
(-z)-(-z):0
(-z)*(-z):unsupported

(x)+(y + x):y + 2*x
(x)-(y + x):-y
(x)*(y + x):unsupported
(y)+(y + x):2*y + x
(y)-(y + x):-x
(y)*(y + x):unsupported

(-z)+(y + x):y + x - z
(-z)-(y + x):-y - x - z
(-z)*(y + x):unsupported
+(y + x):y + x
-(y + x):-y - x
(y + x)+(x):y + 2*x
(y + x)-(x):y
(y + x)*(x):unsupported
(y + x)+(y):2*y + x
(y + x)-(y):x
(y + x)*(y):unsupported
(y + x)+(-z):y + x - z
(y + x)-(-z):y + x + z
(y + x)*(-z):unsupported
(y + x)+(y + x):2*y + 2*x
(y + x)-(y + x):0
(y + x)*(y + x):unsupported

(x)+(2*x):3*x
(x)-(2*x):-x
(x)*(2*x):unsupported
(y)+(2*x):y + 2*x
(y)-(2*x):y - 2*x
(y)*(2*x):unsupported
(-z)+(2*x):2*x - z
(-z)-(2*x):-2*x - z
(-z)*(2*x):unsupported
(y + x)+(2*x):y + 3*x
(y + x)-(2*x):y - x
(y + x)*(2*x):unsupported
+(2*x):2*x
-(2*x):-2*x
(2*x)+(x):3*x
(2*x)-(x):x

(2*x)*(x):unsupported
(2*x)+(y):y + 2*x
(2*x)-(y):-y + 2*x
(2*x)*(y):unsupported
(2*x)+(-z):2*x - z
(2*x)-(-z):2*x + z
(2*x)*(-z):unsupported
(2*x)+(y + x):y + 3*x
(2*x)-(y + x):-y + x
(2*x)*(y + x):unsupported
(2*x)+(2*x):4*x
(2*x)-(2*x):0
(2*x)*(2*x):unsupported

(x)+(0):x
(x)-(0):x
(x)*(0):0
(y)+(0):y
(y)-(0):y
(y)*(0):0
(-z)+(0):-z
(-z)-(0):-z
(-z)*(0):0

(y + x)+(0):y + x
(y + x)-(0):y + x
(y + x)*(0):0
(2*x)+(0):2*x
(2*x)-(0):2*x
(2*x)*(0):0
(0)+(x):x
(0)-(x):-x
(0)*(x):0
(0)+(y):y
(0)-(y):-y
(0)*(y):0
(0)+(-z):-z
(0)-(-z):z
(0)*(-z):0
(0)+(y + x):y + x
(0)-(y + x):-y - x
(0)*(y + x):0
(0)+(2*x):2*x
(0)-(2*x):-2*x
(0)*(2*x):0

(x)+(2):2 + x
(x)-(2):-2 + x
(x)*(2):2*x
(y)+(2):y + 2
(y)-(2):y - 2
(y)*(2):2*y
(-z)+(2):2 - z
(-z)-(2):-2 - z

(y + x)+(2):y + x + 2
(y + x)-(2):y + x - 2
(y + x)*(2):2*y + 2*x
(2*x)+(2):2*x + 2
(2*x)-(2):2*x - 2


(2)+(x):2 + x
(2)-(x):-x + 2
(2)*(x):2*x
(2)+(y):y + 2
(2)-(y):-y + 2
(2)*(y):2*y
(2)+(-z):2 - z
(2)-(-z):2 + z

(2)*(-z):-2*z
(2)+(y + x):y + x + 2
(2)-(y + x):-y - x + 2
(2)*(y + x):2*y + 2*x
(2)+(2*x):2*x + 2
(2)-(2*x):-2*x + 2
(2)*(2*x):4*x
(2)+(0):2
(2)-(0):2
(2)*(0):0
(2)+(2):4
(2)-(2):0
(2)*(2):4

(-z)*(2):-2*z
(2*x)*(2):4*x

(x)+(1):1 + x
(x)-(1):-1 + x
(x)*(1):x
(y)+(1):y + 1
(y)-(1):y - 1
(y)*(1):y
(-z)+(1):1 - z
(-z)-(1):-1 - z
(-z)*(1):-z
(y + x)+(1):y + x + 1
(y + x)-(1):y + x - 1
(y + x)*(1):y + x
(2*x)+(1):2*x + 1
(2*x)-(1):2*x - 1
(2*x)*(1):2*x
(0)+(1):1
(0)-(1):-1
(0)*(1):0
(2)+(1):3
(2)-(1):1
(2)*(1):2
+(1):1

-(1):-1
(1)+(x):1 + x
(1)-(x):-x + 1
(1)*(x):x
(1)+(y):y + 1
(1)-(y):-y + 1
(1)*(y):y
(1)+(-z):1 - z
(1)-(-z):1 + z
(1)*(-z):-z
(1)+(y + x):y + x + 1
(1)-(y + x):-y - x + 1
(1)*(y + x):y + x
(1)+(2*x):2*x + 1
(1)-(2*x):-2*x + 1
(1)*(2*x):2*x
(1)+(0):1
(1)-(0):1
(1)*(0):0
(1)+(2):3
(1)-(2):-1
(1)*(2):2
(1)+(1):2
(1)-(1):0
(1)*(1):1

(x)+(-x + z):z
(x)-(-x + z):2*x - z
(x)*(-x + z):unsupported
(y)+(-x + z):y - x + z
(y)-(-x + z):y + x - z
(y)*(-x + z):unsupported
(-z)+(-x + z):-x
(-z)-(-x + z):x - 2*z
(-z)*(-x + z):unsupported
(y + x)+(-x + z):y + z
(y + x)-(-x + z):y + 2*x - z
(y + x)*(-x + z):unsupported
(2*x)+(-x + z):x + z
(2*x)-(-x + z):3*x - z
(2*x)*(-x + z):unsupported
(0)+(-x + z):-x + z
(0)-(-x + z):x - z
(0)*(-x + z):0


(2)+(-x + z):-x + z + 2
(2)-(-x + z):x - z + 2
(2)*(-x + z):-2*x + 2*z
(1)+(-x + z):-x + z + 1
(1)-(-x + z):x - z + 1
(1)*(-x + z):-x + z
(2)+(-x + z):-x + z + 2
(2)-(-x + z):x - z + 2
(2)*(-x + z):-2*x + 2*z
+(-x + z):-x + z
-(-x + z):x - z
(-x + z)+(x):z
(-x + z)-(x):-2*x + z
(-x + z)*(x):unsupported
(-x + z)+(y):y - x + z
(-x + z)-(y):-y - x + z
(-x + z)*(y):unsupported
(-x + z)+(-z):-x
(-x + z)-(-z):-x + 2*z
(-x + z)*(-z):unsupported
(-x + z)+(y + x):y + z
(-x + z)-(y + x):-y - 2*x + z
(-x + z)*(y + x):unsupported
(-x + z)+(2*x):x + z
(-x + z)-(2*x):-3*x + z

(-x + z)*(2*x):unsupported
(-x + z)+(0):-x + z
(-x + z)-(0):-x + z
(-x + z)*(0):0
(-x + z)+(2):-x + z + 2
(-x + z)-(2):-x + z - 2
(-x + z)*(2):-2*x + 2*z
(-x + z)+(1):-x + z + 1
(-x + z)-(1):-x + z - 1
(-x + z)*(1):-x + z
(-x + z)+(2):-x + z + 2
(-x + z)-(2):-x + z - 2
(-x + z)*(2):-2*x + 2*z
(-x + z)+(0):-x + z
(-x + z)-(0):-x + z
(-x + z)*(0):0
(-x + z)+(-x + z):-2*x + 2*z
(-x + z)-(-x + z):0
(-x + z)*(-x + z):unsupported

(x)+(x - x):x
(x)-(x - x):x
(x)*(x - x):unsupported
(y)+(x - x):y
(y)-(x - x):y
(y)*(x - x):unsupported
(-z)+(x - x):-z
(-z)-(x - x):-z
(-z)*(x - x):unsupported
(y + x)+(x - x):y + x
(y + x)-(x - x):y + x
(y + x)*(x - x):unsupported
(2*x)+(x - x):2*x
(2*x)-(x - x):2*x
(2*x)*(x - x):unsupported
(0)+(x - x):0
(0)-(x - x):0
(0)*(x - x):0
(2)+(x - x):2
(2)-(x - x):2
(2)*(x - x):0
(1)+(x - x):1
(1)-(x - x):1
(1)*(x - x):0
(2)+(x - x):2
(2)-(x - x):2
(2)*(x - x):0
(0)+(x - x):0
(0)-(x - x):0
(0)*(x - x):0

(-x + z)+(x - x):-x + z
(-x + z)-(x - x):-x + z
(-x + z)*(x - x):unsupported
+(x - x):0
-(x - x):0
(x - x)+(x):x
(x - x)-(x):-x
(x - x)*(x):unsupported
(x - x)+(y):y
(x - x)-(y):-y
(x - x)*(y):unsupported
(x - x)+(-z):-z
(x - x)-(-z):z
(x - x)*(-z):unsupported
(x - x)+(y + x):y + x
(x - x)-(y + x):-y - x

(x - x)*(y + x):unsupported
(x - x)+(2*x):2*x
(x - x)-(2*x):-2*x
(x - x)*(2*x):unsupported
(x - x)+(0):0
(x - x)-(0):0
(x - x)*(0):0
(x - x)+(2):2
(x - x)-(2):-2
(x - x)*(2):0
(x - x)+(1):1
(x - x)-(1):-1
(x - x)*(1):0
(x - x)+(2):2
(x - x)-(2):-2
(x - x)*(2):0
(x - x)+(0):0
(x - x)-(0):0
(x - x)*(0):0
(x - x)+(-x + z):-x + z
(x - x)-(-x + z):x - z
(x - x)*(-x + z):unsupported
(x - x)*(y + x):unsupported
(x - x)+(x - x):0
(x - x)-(x - x):0
(x - x)*(x - x):unsupported

(x)+(x + x):3*x
(x)-(x + x):-x
(x)*(x + x):unsupported
(y)+(x + x):y + 2*x
(y)-(x + x):y - 2*x
(y)*(x + x):unsupported
(-z)+(x + x):2*x - z
(-z)-(x + x):-2*x - z
(-z)*(x + x):unsupported
(y + x)+(x + x):y + 3*x
(y + x)-(x + x):y - x
(y + x)*(x + x):unsupported
(2*x)+(x + x):4*x
(2*x)-(x + x):0
(2*x)*(x + x):unsupported
(0)+(x + x):2*x
(0)-(x + x):-2*x
(0)*(x + x):0
(2)+(x + x):2*x + 2
(2)-(x + x):-2*x + 2
(2)*(x + x):4*x
(1)+(x + x):2*x + 1
(1)-(x + x):-2*x + 1

(1)*(x + x):2*x
(-x + z)+(x + x):x + z
(-x + z)-(x + x):-3*x + z
(-x + z)*(x + x):unsupported
(x - x)+(x + x):2*x
(x - x)-(x + x):-2*x
(x - x)*(x + x):unsupported
+(x + x):2*x
-(x + x):-2*x
(x + x)+(x):3*x
(x + x)-(x):x
(x + x)*(x):unsupported
(x + x)+(y):y + 2*x
(x + x)-(y):-y + 2*x
(x + x)*(y):unsupported
(x + x)+(-z):2*x - z
(x + x)-(-z):2*x + z
(x + x)*(-z):unsupported
(x + x)+(y + x):y + 3*x
(x + x)-(y + x):-y + x
(x + x)*(y + x):unsupported
(x + x)+(2*x):4*x
(x + x)-(2*x):0
(x + x)*(2*x):unsupported

(x + x)+(0):2*x
(x + x)-(0):2*x
(x + x)*(0):0
(x + x)+(2):2*x + 2
(x + x)-(2):2*x - 2
(x + x)*(2):4*x
(x + x)+(1):2*x + 1
(x + x)-(1):2*x - 1
(x + x)*(1):2*x
(x + x)+(2):2*x + 2
(x + x)-(2):2*x - 2
(x + x)*(2):4*x
(x + x)+(0):2*x
(x + x)-(0):2*x
(x + x)*(0):0
(x + x)+(-x + z):x + z
(x + x)-(-x + z):3*x - z
(x + x)*(-x + z):unsupported
(x + x)+(x - x):2*x
(x + x)-(x - x):2*x
(x + x)*(x - x):unsupported
(x + x)+(x + x):4*x
(x + x)-(x + x):0
(x + x)*(x + x):unsupported

(x)+(x - z):2*x - z
(x)-(x - z):z
(x)*(x - z):unsupported
(y)+(x - z):y + x - z
(y)-(x - z):y - x + z
(y)*(x - z):unsupported
(-z)+(x - z):x - 2*z
(-z)-(x - z):-x
(-z)*(x - z):unsupported
(y + x)+(x - z):y + 2*x - z
(y + x)-(x - z):y + z
(y + x)*(x - z):unsupported
(2*x)+(x - z):3*x - z
(2*x)-(x - z):x + z

(2*x)*(x - z):unsupported
(0)+(x - z):x - z
(0)-(x - z):-x + z
(0)*(x - z):0
(2)+(x - z):x - z + 2
(2)-(x - z):-x + z + 2
(2)*(x - z):2*x - 2*z
(1)+(x - z):x - z + 1
(1)-(x - z):-x + z + 1
(1)*(x - z):x - z
(2)+(x - z):x - z + 2
(2)-(x - z):-x + z + 2
(2)*(x - z):2*x - 2*z
(0)+(x - z):x - z
(0)-(x - z):-x + z
(0)*(x - z):0
(-x + z)+(x - z):0
(-x + z)-(x - z):-2*x + 2*z
(-x + z)*(x - z):unsupported
(x - x)+(x - z):x - z
(x - x)-(x - z):-x + z
(x - x)*(x - z):unsupported
(x + x)+(x - z):3*x - z
(x + x)-(x - z):x + z
(x + x)*(x - z):unsupported
(2*x)*(x - z):unsupported

+(x - z):x - z
-(x - z):-x + z
(x - z)+(x):2*x - z
(x - z)-(x):-z
(x - z)*(x):unsupported
(x - z)+(y):y + x - z
(x - z)-(y):-y + x - z
(x - z)*(y):unsupported
(x - z)+(-z):x - 2*z
(x - z)-(-z):x
(x - z)*(-z):unsupported
(x - z)+(y + x):y + 2*x - z
(x - z)-(y + x):-y - z
(x - z)*(y + x):unsupported
(x - z)+(2*x):3*x - z
(x - z)-(2*x):-x - z
(x - z)*(2*x):unsupported
(x - z)+(0):x - z
(x - z)-(0):x - z
(x - z)*(0):0
(x - z)+(2):x - z + 2
(x - z)-(2):x - z - 2
(x - z)*(2):2*x - 2*z
(x - z)+(1):x - z + 1
(x - z)-(1):x - z - 1
(x - z)*(1):x - z
(x - z)+(2):x - z + 2
(x - z)-(2):x - z - 2
(x - z)*(2):2*x - 2*z
(x - z)+(0):x - z
(x - z)-(0):x - z
(x - z)*(0):0

(x - z)+(-x + z):0
(x - z)-(-x + z):2*x - 2*z
(x - z)*(-x + z):unsupported
(x - z)+(x - x):x - z
(x - z)-(x - x):x - z
(x - z)*(x - x):unsupported
(x - z)+(x + x):3*x - z
(x - z)-(x + x):-x - z
(x - z)*(x + x):unsupported
(x - z)+(x - z):2*x - 2*z
(x - z)-(x - z):0
(x - z)*(x - z):unsupported
'''

expected_unevaluated_results_additive_group = '''
+(x):x
-(x):-x
(x)+(x):x + x
(x)-(x):x + -x
(x)*(x):unsupported
(x)+(y):x + y
(x)-(y):x + -y
(x)*(y):unsupported
+(y):y
-(y):-y
(y)+(x):y + x
(y)-(x):y + -x
(y)*(x):unsupported
(y)+(y):y + y
(y)-(y):y + -y
(y)*(y):unsupported

(x)+(-z):x + -z
(x)-(-z):x + --z
(x)*(-z):unsupported
(y)-(-z):y + --z
(y)*(-z):unsupported
+(-z):-z
-(-z):--z
(-z)+(x):-z + x
(-z)*(x):unsupported
(-z)+(y):-z + y
(-z)*(y):unsupported
(-z)-(-z):-z + --z
(-z)*(-z):unsupported

(y)+(-z):y + -z
(-z)-(x):-z + -x
(-z)-(y):-z + -y
(-z)+(-z):-z + -z

(x)+(x + y):x + x + y
(x)-(x + y):x + -(x + y)
(x)*(x + y):unsupported
(y)+(x + y):y + x + y
(y)-(x + y):y + -(x + y)
(y)*(x + y):unsupported
(-z)+(x + y):-z + x + y
(-z)-(x + y):-z + -(x + y)
(-z)*(x + y):unsupported
+(x + y):x + y
-(x + y):-(x + y)
(x + y)+(x):x + y + x
(x + y)-(x):x + y + -x
(x + y)*(x):unsupported
(x + y)+(y):x + y + y
(x + y)-(y):x + y + -y
(x + y)*(y):unsupported
(x + y)+(-z):x + y + -z
(x + y)-(-z):x + y + --z
(x + y)*(-z):unsupported
(x + y)+(x + y):x + y + x + y
(x + y)-(x + y):x + y + -(x + y)
(x + y)*(x + y):unsupported

(x)+(3*x):x + 3*x
(x)-(3*x):x + -(3*x)
(x)*(3*x):unsupported
(y)+(3*x):y + 3*x
(y)-(3*x):y + -(3*x)
(y)*(3*x):unsupported
(-z)+(3*x):-z + 3*x
(-z)-(3*x):-z + -(3*x)
(-z)*(3*x):unsupported
(x + y)+(3*x):x + y + 3*x
(x + y)-(3*x):x + y + -(3*x);x + y + -(x*3)
(x + y)*(3*x):unsupported
+(3*x):3*x
-(3*x):-(3*x)
(3*x)+(x):3*x + x
(3*x)-(x):3*x + -x
(3*x)*(x):unsupported
(3*x)+(y):3*x + y
(3*x)-(y):3*x + -y
(3*x)*(y):unsupported
(3*x)+(-z):3*x + -z
(3*x)-(-z):3*x + --z
(3*x)*(-z):unsupported
(3*x)+(x + y):3*x + x + y
(3*x)-(x + y):3*x + -(x + y)
(3*x)*(x + y):unsupported
(3*x)+(3*x):3*x + 3*x
(3*x)-(3*x):3*x + -(3*x)
(3*x)*(3*x):unsupported

(x)+(0):x + 0
(x)-(0):x + -0
(x)*(0):x*0
(y)+(0):y + 0
(y)-(0):y + -0
(y)*(0):y*0
(-z)+(0):-z + 0
(-z)-(0):-z + -0
(-z)*(0):-z*0
(x + y)+(0):x + y + 0
(x + y)-(0):x + y + -0
(x + y)*(0):(x + y)*0
(3*x)+(0):3*x + 0
(3*x)-(0):3*x + -0
(3*x)*(0):3*x*0
+(0):0
-(0):-0
(0)+(x):0 + x
(0)-(x):0 + -x

(0)+(y):0 + y
(0)-(y):0 + -y
(0)+(-z):0 + -z
(0)-(-z):0 + --z
(0)+(x + y):0 + x + y
(0)-(x + y):0 + -(x + y)
(0)+(3*x):0 + 3*x
(0)-(3*x):0 + -(3*x)
(0)+(0):0 + 0
(0)-(0):0 + -0
(0)*(0):0*0

(0)*(x):0*x
(0)*(y):0*y
(0)*(-z):0*-z
(0)*(x + y):0*(x + y)
(0)*(3*x):0*3*x

'''
    
expected_results_additive_group = '''
+(x):x
-(x):-x
(x)+(x):2*x
(x)-(x):0
(x)+(y):x + y
(x)-(y):x - y
+(y):y
-(y):-y
(y)+(x):y + x
(y)-(x):y - x
(y)+(y):2*y
(y)-(y):0
(x)+(-z):x - z
(x)-(-z):x + z
(y)+(-z):y - z
(y)-(-z):y + z
+(-z):-z
-(-z):z

(-z)+(x):-z + x
(-z)-(x):-z - x
(-z)+(y):-z + y
(-z)-(y):-z - y
(-z)+(-z):-2*z
(-z)-(-z):0

(x)+(2*x):3*x
(x)-(2*x):-x
(y)+(2*x):y + 2*x
(y)-(2*x):y - 2*x
(-z)+(2*x):-z + 2*x
(-z)-(2*x):-z - 2*x
+(2*x):2*x
-(2*x):-2*x
(2*x)+(x):3*x
(2*x)-(x):x
(2*x)+(y):2*x + y
(2*x)-(y):2*x - y
(2*x)+(-z):2*x - z
(2*x)-(-z):2*x + z
(2*x)+(2*x):4*x
(2*x)-(2*x):0

(x)+(x + y):2*x + y
(x)-(x + y):x - y - x
(y)+(x + y):y + x + y
(y)-(x + y):-x
(-z)+(x + y):-z + x + y
(-z)-(x + y):-z - y - x
(2*x)+(x + y):3*x + y
(2*x)-(x + y):2*x - y - x
+(x + y):x + y
-(x + y):-y - x
(x + y)+(x):x + y + x
(x + y)-(x):x + y - x
(x + y)+(y):x + 2*y
(x + y)-(y):x
(x + y)+(-z):x + y - z
(x + y)-(-z):x + y + z
(x + y)+(2*x):x + y + 2*x
(x + y)-(2*x):x + y - 2*x
(x + y)+(x + y):x + y + x + y
(x + y)-(x + y):0

(x)+(y + x):x + y + x
(x)-(y + x):-y
(y)+(y + x):2*y + x
(y)-(y + x):y - x - y
(-z)+(y + x):-z + y + x
(-z)-(y + x):-z - x - y
(2*x)+(y + x):2*x + y + x
(2*x)-(y + x):x - y
(x + y)+(y + x):x + 2*y + x
(x + y)-(y + x):x + y - x - y
+(y + x):y + x
-(y + x):-x - y
(y + x)+(x):y + 2*x
(y + x)-(x):y
(y + x)+(y):y + x + y
(y + x)-(y):y + x - y
(y + x)+(-z):y + x - z
(y + x)-(-z):y + x + z
(y + x)+(2*x):y + 3*x
(y + x)-(2*x):y - x
(y + x)+(x + y):y + 2*x + y
(y + x)-(x + y):y + x - y - x
(y + x)+(y + x):y + x + y + x
(y + x)-(y + x):0

(x)+(0):x
(x)-(0):x
(y)+(0):y
(y)-(0):y
(-z)+(0):-z
(-z)-(0):-z
(2*x)+(0):2*x
(2*x)-(0):2*x
(x + y)+(0):x + y
(x + y)-(0):x + y
(y + x)+(0):y + x
(y + x)-(0):y + x
+(0):0
-(0):0
(0)+(x):x
(0)-(x):-x
(0)+(y):y
(0)-(y):-y
(0)+(-z):-z
(0)-(-z):z
(0)+(2*x):2*x
(0)-(2*x):-2*x
(0)+(x + y):x + y
(0)-(x + y):-y - x
(0)+(y + x):y + x
(0)-(y + x):-x - y
(0)+(0):0
(0)-(0):0

(x)*(x):unsupported
(x)*(y):unsupported
(x)*(-z):unsupported
(x)*(2*x):unsupported
(x)*(x + y):unsupported
(x)*(y + x):unsupported

(x)*(0):0
(x)*(0):0
(y)*(x):unsupported
(y)*(y):unsupported
(y)*(-z):unsupported
(y)*(2*x):unsupported
(y)*(x + y):unsupported
(y)*(y + x):unsupported
(y)*(0):0
(y)*(0):0
(-z)*(x):unsupported
(-z)*(y):unsupported
(-z)*(-z):unsupported
(-z)*(2*x):unsupported
(-z)*(x + y):unsupported
(-z)*(y + x):unsupported

(-z)*(0):0
(-z)*(0):0
(2*x)*(x):unsupported
(2*x)*(y):unsupported
(2*x)*(-z):unsupported
(2*x)*(2*x):unsupported
(2*x)*(x + y):unsupported
(2*x)*(y + x):unsupported
(2*x)*(0):0
(2*x)*(0):0
(x + y)*(x):unsupported
(x + y)*(y):unsupported
(x + y)*(-z):unsupported
(x + y)*(2*x):unsupported
(x + y)*(x + y):unsupported
(x + y)*(y + x):unsupported
(x + y)*(0):0
(x + y)*(0):0
(y + x)*(x):unsupported
(y + x)*(y):unsupported
(y + x)*(-z):unsupported
(y + x)*(2*x):unsupported
(y + x)*(x + y):unsupported
(y + x)*(y + x):unsupported
(y + x)*(0):0
(y + x)*(0):0
(0)*(x):0
(0)*(y):0
(0)*(-z):0
(0)*(2*x):0
(0)*(x + y):0
(0)*(y + x):0
(0)*(0):0
(0)*(0):0

(x)+(2):x + 2
(x)-(2):x - 2
(x)*(2):2*x
(y)+(2):y + 2
(y)-(2):y - 2
(y)*(2):2*y
(-z)+(2):-z + 2
(-z)-(2):-z - 2
(0)+(2):2
(0)-(2):-2
(0)*(2):0
(2)+(x):2 + x
(2)-(x):2 - x
(2)*(x):2*x
(2)+(y):2 + y
(2)-(y):2 - y
(2)*(y):2*y
(2)+(-z):2 - z
(2)-(-z):2 + z
(2)+(2*x):2 + 2*x
(2)-(2*x):2 - 2*x
(2)*(2*x):4*x
(2)+(x + y):2 + x + y
(2)-(x + y):2 - y - x
(2)+(y + x):2 + y + x
(2)-(y + x):2 - x - y
(2)+(0):2
(2)-(0):2
(2)*(0):0
(2*x)+(2):2*x + 2
(2*x)-(2):2*x - 2
(2*x)*(2):4*x
(x + y)+(2):x + y + 2
(x + y)-(2):x + y - 2
(y + x)+(2):y + x + 2
(y + x)-(2):y + x - 2
(-z)*(2):-2*z
(2)*(-z):-2*z
(x + y)*(2):x + y + x + y
(y + x)*(2):y + x + y + x
(2)*(x + y):x + y + x + y
(2)*(y + x):y + x + y + x


+(x + y + x):x + y + x
-(x + y + x):-x - y - x
(x + y + x)+(x + y + x):x + y + 2*x + y + x
(x + y + x)-(x + y + x):0
(x + y + x)*(x + y + x):unsupported
(x + y + x)+(3):x + y + x + 3
(x + y + x)-(3):x + y + x - 3
(x + y + x)*(3):x + y + 2*x + y + 2*x + y + x
(3)+(x + y + x):3 + x + y + x
(3)-(x + y + x):3 - x - y - x
(3)*(x + y + x):x + y + 2*x + y + 2*x + y + x

(x + y + x)+(-4):x + y + x - 4
(x + y + x)-(-4):x + y + x + 4
(x + y + x)*(-4):-x - y - 2*x - y - 2*x - y - 2*x - y - x
(-4)+(x + y + x):-4 + x + y + x
(-4)-(x + y + x):-4 - x - y - x
(-4)*(x + y + x):-x - y - 2*x - y - 2*x - y - 2*x - y - x

(x + y + x)+(x + y + z):x + y + 2*x + y + z
(x + y + x)-(x + y + z):x + y + x - z - y - x
(x + y + x)*(x + y + z):unsupported
(3)+(x + y + z):3 + x + y + z
(3)-(x + y + z):3 - z - y - x
(3)*(x + y + z):x + y + z + x + y + z + x + y + z
(-4)+(x + y + z):-4 + x + y + z
(-4)-(x + y + z):-4 - z - y - x
(-4)*(x + y + z):-z - y - x - z - y - x - z - y - x - z - y - x
+(x + y + z):x + y + z
-(x + y + z):-z - y - x
(x + y + z)+(x + y + x):x + y + z + x + y + x
(x + y + z)-(x + y + x):x + y + z - x - y - x
(x + y + z)*(x + y + x):unsupported
(x + y + z)+(3):x + y + z + 3
(x + y + z)-(3):x + y + z - 3
(x + y + z)*(3):x + y + z + x + y + z + x + y + z
(x + y + z)+(-4):x + y + z - 4
(x + y + z)-(-4):x + y + z + 4
(x + y + z)*(-4):-z - y - x - z - y - x - z - y - x - z - y - x
(x + y + z)+(x + y + z):x + y + z + x + y + z
(x + y + z)-(x + y + z):0
(x + y + z)*(x + y + z):unsupported


(x + y + x)+(x - y):x + y + 2*x - y
(x + y + x)-(x - y):x + y + x + y - x
(x + y + x)*(x - y):unsupported
(3)+(x - y):3 + x - y
(3)-(x - y):3 + y - x
(3)*(x - y):x - y + x - y + x - y
(-4)+(x - y):-4 + x - y
(-4)-(x - y):-4 + y - x
(-4)*(x - y):y - x + y - x + y - x + y - x
(x + y + z)+(x - y):x + y + z + x - y
(x + y + z)-(x - y):x + y + z + y - x
(x + y + z)*(x - y):unsupported
+(x - y):x - y
-(x - y):y - x
(x - y)+(x + y + x):x - y + x + y + x
(x - y)-(x + y + x):x - y - x - y - x
(x - y)*(x + y + x):unsupported
(x - y)+(3):x - y + 3
(x - y)-(3):x - y - 3
(x - y)*(3):x - y + x - y + x - y
(x - y)+(-4):x - y - 4
(x - y)-(-4):x - y + 4
(x - y)*(-4):y - x + y - x + y - x + y - x
(x - y)+(x + y + z):x - y + x + y + z
(x - y)-(x + y + z):x - y - z - y - x
(x - y)*(x + y + z):unsupported
(x - y)+(x - y):x - y + x - y
(x - y)-(x - y):0
(x - y)*(x - y):unsupported

(x)+(1):x + 1
(x)-(1):x - 1
(x)*(1):x
(y)+(1):y + 1
(y)-(1):y - 1
(y)*(1):y
(-z)+(1):-z + 1
(-z)-(1):-z - 1
(-z)*(1):-z
(2*x)+(1):2*x + 1
(2*x)-(1):2*x - 1
(2*x)*(1):2*x
(x + y)+(1):x + y + 1
(x + y)-(1):x + y - 1
(x + y)*(1):x + y
(y + x)+(1):y + x + 1
(y + x)-(1):y + x - 1
(y + x)*(1):y + x

(0)+(1):1
(0)-(1):-1
(0)*(1):0
(0)+(1):1
(0)-(1):-1
(0)*(1):0
(2)+(1):3
(2)-(1):1
(2)*(1):2
+(1):1
-(1):-1
(1)+(x):1 + x
(1)-(x):1 - x
(1)*(x):x
(1)+(y):1 + y
(1)-(y):1 - y
(1)*(y):y

(1)+(-z):1 - z
(1)-(-z):1 + z
(1)*(-z):-z
(1)+(2*x):1 + 2*x
(1)-(2*x):1 - 2*x
(1)*(2*x):2*x
(1)+(x + y):1 + x + y
(1)-(x + y):1 - y - x
(1)*(x + y):x + y
(1)+(y + x):1 + y + x
(1)-(y + x):1 - x - y
(1)*(y + x):y + x


'''


expected_results_numbers = '''
(0)+(0):0
(0)-(0):0
(0)*(0):0
+(0):0
-(0):0
(0)+(0):0
(0)-(0):0
(0)*(0):0
(0)+(0):0
(0)-(0):0
(0)*(0):0
(1)+(0):1
(1)-(0):1
(1)*(0):0
(1)+(0):1
(1)-(0):1
(1)*(0):0
(1)+(2):3
(1)-(2):-1
(1)*(2):2
(1)+(1):2
(1)-(1):0
(1)*(1):1

(2)+(2):4
(2)-(2):0
(2)*(2):4
+(2):2
-(2):-2
(2)+(2):4
(2)-(2):0
(2)*(2):4
(2)+(2):4
(2)-(2):0
(2)*(2):4
(0)+(2):2
(0)-(2):-2
(0)*(2):0
+(2):2
-(2):-2
'''

expected_results_numbers_unevaluated = '''
+(0):0
-(0):-0
(0)+(0):0 + 0
(0)-(0):0 + -0
(0)*(0):0*0
(0)+(2):0 + 2
(0)-(2):0 + -2
(0)*(2):0*2

(2)+(0):2 + 0
(2)-(0):2 + -0
(2)*(0):2*0

'''
