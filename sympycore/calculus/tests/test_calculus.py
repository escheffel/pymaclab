
from sympycore.calculus import *

def test_factorial():
    assert isinstance(factorial(2), Calculus)
    assert factorial(2)==2
    assert factorial(3)==6
    assert factorial(4)==24
    assert factorial(1)==1
    assert factorial(0)==1
    
    
