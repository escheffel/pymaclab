
from sympycore import *

def test_eye():
    assert eye(3).tolist()==[[1,0,0],[0,1,0],[0,0,1]]
    assert eye(3,2).tolist()==[[1,0],[0,1],[0,0]]
    assert eye(2,3).tolist()==[[1,0,0],[0,1,0]]

    assert eye(3,k=1).tolist()==[[0,1,0],[0,0,1],[0,0,0]]
    assert eye(3,k=2).tolist()==[[0,0,1],[0,0,0],[0,0,0]]
    assert eye(3,k=3).tolist()==[[0,0,0],[0,0,0],[0,0,0]]
    assert eye(3,k=-1).tolist()==[[0,0,0],[1,0,0],[0,1,0]]

def test_concatenate():

    m = Matrix([[1,2],[3,4]])
    assert concatenate(m,m).tolist()==[[1,2,1,2],
                                       [3,4,3,4]]
    assert concatenate(m,m,[5,6],[[7,8]]).tolist()==[[1,2,1,2,5,7,8],
                                                     [3,4,3,4,6,0,0]]

    assert concatenate(m,m,axis=1).tolist()==[[1,2],[3,4],[1,2],[3,4]]
    assert concatenate(m,m,[5,6],[[7,8]],axis=1).tolist()==[[1,2],
                                                            [3,4],
                                                            [1,2],
                                                            [3,4],
                                                            [5,0],
                                                            [6,0],
                                                            [7,8]]

    assert concatenate(m,m,diagonal=True).tolist()==[[1,2,0,0],
                                                     [3,4,0,0],
                                                     [0,0,1,2],
                                                     [0,0,3,4]]
    assert concatenate(m,m,[5,6],[[7,8]], diagonal=True).tolist()\
           ==[[1,2,0,0,0,0,0],
              [3,4,0,0,0,0,0],
              [0,0,1,2,0,0,0],
              [0,0,3,4,0,0,0],
              [0,0,0,0,5,0,0],
              [0,0,0,0,6,0,0],
              [0,0,0,0,0,7,8]]
