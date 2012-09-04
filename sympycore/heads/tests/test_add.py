
from sympycore import *

def test_combine_add_list():
    combine_add_list = heads.ADD.combine_add_list
    A = Ring
    x,y,z = map(A, 'xyz')
    assert combine_add_list(A, [x]) == [x]
    assert combine_add_list(A, [x,-x]) == []
    assert combine_add_list(A, [x,y]) == [x,y]
    assert combine_add_list(A, [x,x]) == [2*x]
    assert combine_add_list(A, [x,y,x]) == [x,y,x]
    assert combine_add_list(A, [x,y,-y]) == [x]
    assert combine_add_list(A, [x,y,-y,x]) == [2*x]
    assert combine_add_list(A, [-x,y,-y,x]) == []
    assert combine_add_list(A, [x,3*x]) == [4*x]
    assert combine_add_list(A, [2*x,3*x]) == [5*x]
    assert combine_add_list(A, [2*x,3*x,x]) == [6*x]
