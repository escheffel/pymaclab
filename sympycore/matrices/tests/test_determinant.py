
from sympycore import *

def test_det1():
    assert Matrix([2]).det()==2

def test_det2():
    assert Matrix([[1,2],[3,4]]).det()==-2
    assert Matrix([[1,2],[2,4]]).det()==0

def test_det3():
    assert Matrix([[1,2,3],[4,5,6],[7,8,9]]).det()==0
    assert Matrix([[1,2,3],[4,5,6],[7,8,10]]).det()==-3

def test_det4():
    assert Matrix([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]).det()==0
    assert Matrix([[10,2,3,4],[5,60,7,8],[9,10,110,12],[13,14,15,16]]).det()==561168
