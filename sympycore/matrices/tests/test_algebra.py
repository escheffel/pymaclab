
from sympycore import *

def test_Matrix1():
    a = Matrix(2)
    assert a.rows==2
    assert a.cols==2
    assert a.tolist()==[[0,0], [0,0]]

    a = Matrix([1,2])
    assert a.rows==2
    assert a.cols==1
    assert a.tolist()==[[1], [2]]

    a = Matrix([1,2], diagonal=True)
    assert a.rows==2
    assert a.cols==2
    assert a.tolist()==[[1,0], [0,2]]

    a = Matrix([0,1], permutation=True)
    assert a.rows==2
    assert a.cols==2
    assert a.tolist()==[[1,0], [0,1]]

    a = Matrix([1,0], permutation=True)
    assert a.rows==2
    assert a.cols==2
    assert a.tolist()==[[0,1], [1,0]]

    a = Matrix([1,[2]])
    assert a.rows==2
    assert a.cols==1
    assert a.tolist()==[[1], [2]]

    a = Matrix([[1,2], [3,4]])
    assert a.rows==2
    assert a.cols==2
    assert a.tolist()==[[1,2], [3,4]]

    a = Matrix([[1,2,3], [4,5,6]])
    assert a.rows==2
    assert a.cols==3
    assert a.tolist()==[[1,2,3], [4,5,6]]

def test_Matrix_properties():
    a = Matrix([[1,0,0,0],[1,1,0,0],[1,1,1,0],[1,1,1,1]])
    assert a.is_square
    assert a.is_lower
    assert not a.is_upper
    assert a.T.is_upper
    assert not a.is_row_echelon_form
    a[0,2] = 1
    assert not a.is_lower
    assert not a.is_upper
    assert not a.T.is_upper
    assert Matrix([[1,0],[0,1]]).is_orthogonal
    assert Matrix([[1,0],[0,-1]]).is_orthogonal
    assert Matrix([[0,0,0,1],[0,0,1,0],[1,0,0,0],[0,1,0,0]]).is_orthogonal
    assert not Matrix([[1,0],[2,-1]]).is_orthogonal
    assert not Matrix([[1,0],[0,1],[0,1]]).is_orthogonal

    assert Matrix([[1,0],[0,1]]).is_row_echelon_form
    assert Matrix([[1,0],[0,1]]).is_row_canonical_form
    assert not Matrix([[1,0],[1,1]]).is_row_echelon_form
    assert not Matrix([[1,0],[1,1]]).is_row_canonical_form
    assert Matrix([[1,1],[0,1]]).is_row_echelon_form
    assert not Matrix([[1,1],[0,1]]).is_row_canonical_form
    assert Matrix([[1,1,0],[0,0,1],[0,0,0]]).is_row_echelon_form
    a = Matrix([[1,1,0],[0,0,1],[0,0,0]])
    assert a.is_row_canonical_form
    assert not Matrix([[1,1,0],[0,0,0],[0,0,1]]).is_row_echelon_form
    assert Matrix([[0,1,0],[0,0,1],[0,0,0]]).is_row_echelon_form
    assert Matrix([[0,1,0],[0,0,1],[0,0,0]]).is_row_canonical_form
    assert Matrix([[0,1,0,1],[0,0,1,1],[0,0,0,2]]).is_row_echelon_form
    assert not Matrix([[0,1,0,1],[0,0,1,1],[0,0,0,2]]).is_row_canonical_form
    assert not Matrix([[0,1,0,1],[0,0,1,1],[0,0,0,1]]).is_row_canonical_form
    assert Matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1]]).is_row_canonical_form

    a = Matrix([[1,2,0,0,3,4,0,5,0],[0,0,1,0,2,3,0,4,0],[0,0,0,1,2,3,0,4,0],[0.0,0,0,0,0,0,1,2,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,0]])
    assert a.is_row_echelon_form
    assert a.is_row_canonical_form

def test_Matrix2():
    a = Matrix(2, 3)
    assert a.rows==2
    assert a.cols==3
    assert a.tolist()==[[0,0,0], [0,0,0]]

def test_Matrix3():
    a = Matrix(2, 3, {(0,0):1, (1,2):2})
    assert a.rows==2
    assert a.cols==3
    assert a.tolist()==[[1,0,0], [0,0,2]]

    a = Matrix(2, 3, [1,2,3,4,5,6])
    assert a.rows==2
    assert a.cols==3
    assert a.tolist()==[[1,2,3],[4,5,6]]

def test_T():
    a = Matrix([[1,2], [3,4]]).T
    assert a.tolist()==[[1,3],[2,4]]

def test_get_items():
    a = Matrix([[1,2], [3,4]])
    assert a[0,0]==1
    assert a[0,1]==2
    assert a[1,0]==3
    assert a[1,1]==4
    assert a[0].tolist()==[[1,2]]
    assert a[1].tolist()==[[3,4]]
    assert a.T[0].tolist()==[[1,3]]
    assert a.T[1].tolist()==[[2,4]]
    assert a.T[0].T.tolist()==[[1],[3]]
    assert a.T[1].T.tolist()==[[2],[4]]

def test_get_row_slice():
    a = Matrix([[1,2], [3,4]])
    assert a[:].tolist()==[[1,2], [3,4]]
    assert a[:2].tolist()==[[1,2], [3,4]]
    assert a[:1].tolist()==[[1,2]]
    assert a[1:2].tolist()==[[3,4]]
    assert a[1:].tolist()==[[3,4]]
    assert a[:-1].tolist()==[[1,2]]

def test_get_column_slice():
    a = Matrix([[1,2], [3,4]])
    assert a[:,0].tolist()==[[1],[3]]
    assert a[:,1].tolist()==[[2],[4]]
    assert a[:,:].tolist()==[[1,2],[3,4]]

def test_get_slice():
    a = Matrix([[1,2,3], [4,5,6], [7,8,9]])
    assert a[:2,:2].tolist()==[[1,2],[4,5]]
    assert a[1:,1:].tolist()==[[5,6],[8,9]]
    assert a.T[:2,:2].tolist()==[[1,4],[2,5]]
    assert a[:2,1:].tolist()==[[2,3],[5,6]]

def test_set_item():
    a = Matrix(3,3)
    a[0,0] = 1
    a[0,1] = 2
    a[0,2] = 3
    assert a.tolist()==[[1,2,3],[0,0,0],[0,0,0]]

    a = Matrix(3,3)
    a[0,:] = [[1,2,3]]
    a[2,:] = [[7,8,9]]
    assert a.tolist()==[[1,2,3],[0,0,0],[7,8,9]]

    a = Matrix(3,3)
    a.T[0,:] = [[1,2,3]]
    a.T[2,:] = [[7,8,9]]
    assert a.tolist()==[[1,0,7],[2,0,8],[3,0,9]]

    a = Matrix(3,3)
    a[:,0] = [[1],2,3]
    a[:,2] = [[7],8,9]
    assert a.tolist()==[[1,0,7],[2,0,8],[3,0,9]]

    a = Matrix(3,3)
    b = Matrix([[1,2,3],[7,8,9]])
    a[::2] = b
    assert a.tolist()==[[1,2,3],[0,0,0],[7,8,9]]
    a[:2,:] = 0*b
    assert a.tolist()==[[0,0,0],[0,0,0],[7,8,9]]
    a.T[0,:] = [[1,0,0]]
    assert a.tolist()==[[1,0,0],[0,0,0],[0,8,9]]
    
    a = Matrix([[1,2],[3,4]])
    a.T[1,0] = 22
    a.T[0,1] = 0
    assert a.tolist()==[[1,22],[0,4]]
    a.T[0,1] = 0
    assert a.tolist()==[[1,22],[0,4]]

    a = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    a[:2,:2] = 11
    assert a.tolist()==[[11,11,3],[11,11,6],[7,8,9]]
    a[1:,:] = 0
    assert a.tolist()==[[11,11,3],[0,0,0],[0,0,0]]
    a[:,:1] = 0
    assert a.tolist()==[[0,11,3],[0,0,0],[0,0,0]]
    a.T[::2,1] = 22
    assert a.tolist()==[[0,11,3],[22,0,22],[0,0,0]]
    a.T[:,0] = 0
    assert a.tolist()==[[0,0,0],[22,0,22],[0,0,0]]
    a.T[:1,:] = 0
    assert a.tolist()==[[0,0,0],[0,0,22],[0,0,0]]

def test_get_diagonal():
    a = Matrix([[1,2,3],
                [4,5,6],
                [7,8,9]])
    assert a.D[0].tolist()==[[1],
                             [5],
                             [9]]
    assert a.D[1].tolist()==[[2],
                             [6],
                             ]
    assert a.D[-1].tolist()==[[4],
                             [8],
                             ]

def test_set_diagonal():
    a = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    a.D[0] = 0
    assert a.tolist()==[[0,2,3],
                        [4,0,6],
                        [7,8,0]]
    a.D[1] = 1
    assert a.tolist()==[[0,1,3],
                        [4,0,1],
                        [7,8,0]]
    a.D[-1] = 2
    assert a.tolist()==[[0,1,3],
                        [2,0,1],
                        [7,2,0]]

    a.T.D[1] = 3
    assert a.tolist()==[[0,1,3],
                        [3,0,1],
                        [7,3,0]]

    a.D[0] = [1,2,3]
    assert a.tolist()==[[1,1,3],
                        [3,2,1],
                        [7,3,3]]

    a.D[1] = [7,8]
    assert a.tolist()==[[1,7,3],
                        [3,2,8],
                        [7,3,3]]
    
def test_iadd():
    a = a2 = Matrix([[1,2], [3,4]])
    a += 1
    assert a.tolist()==[[2,3],[4,5]]
    assert a.data is a2.data

    a = a2 = Matrix([[1,2], [3,4]])
    hash(a)
    a += 1
    assert a.tolist()==[[2,3],[4,5]]
    assert a2.tolist()==[[1,2],[3,4]]
    assert a.data is not a2.data

    a = a2 = Matrix([[1,2], [3,4]])
    b = Matrix([[1,-2], [-3,4]])
    a += b
    assert a.tolist()==[[2,0],[0,8]]
    assert a.data is a2.data

    a = a2 = Matrix([[1,2], [3,4]])
    b = Matrix([[1,-2], [-3,4]])
    hash(a)
    a += b
    assert a.tolist()==[[2,0],[0,8]]
    assert a2.tolist()==[[1,2],[3,4]]
    assert a.data is not a2.data

    a = 1
    a2 = Matrix([[1,2], [3,4]])
    a += a2
    assert a.tolist()==[[2,3],[4,5]]

    a += a.T
    assert a.tolist()==[[4,7],[7,10]]

    r = 2
    r += a
    assert r==2+a

def test_add():
    a = Matrix([[1,2], [3,4]])
    assert (a+1).tolist()==[[2,3],[4,5]]
    assert (1+a).tolist()==[[2,3],[4,5]]

    b = Matrix([[1,-2], [-3,4]])
    assert (a+b).tolist()==[[2,0],[0,8]]

def test_isub():
    a = a2 = Matrix([[1,2], [3,4]])
    a -= 1
    assert a.tolist()==[[0,1],[2,3]]
    assert a.data is a2.data

    b = Matrix([[1,-2], [-3,4]])
    a -= b
    assert a.tolist()==[[-1,3],[5,-1]]
    assert a.data is a2.data

    a -= a.T
    assert a.tolist()==[[0,-2],[2,0]]
    assert a.data is a2.data

    r = 2
    r -= a
    assert r==2-a
    

def test_sub():
    a = Matrix([[1,2],[3,4]])
    assert (1-a).tolist()==[[0,-1],[-2,-3]]
    assert (a-1).tolist()==[[0,1],[2,3]]
    b = Matrix([[3,4],[1,2]])
    assert (a-b).tolist()==[[-2,-2],[2,2]]

def test_posneg():
    a = Matrix([[1,2], [3,4]])
    assert (+a).tolist() == [[1,2],[3,4]]
    assert (-a).tolist() == [[-1,-2],[-3,-4]]

def test_imul():
    a = a2 = Matrix([[1,2], [3,4]])
    a *= 2
    assert a.tolist() == [[2,4],[6,8]]
    assert a.data is a2.data

    a = a2 = Matrix([[1,2], [3,4]])
    hash(a)
    a *= 2
    assert a.tolist() == [[2,4],[6,8]]
    assert a2.tolist() == [[1,2],[3,4]]
    assert a.data is not a2.data

    a = a2 = Matrix([[1,2], [3,4]])
    a *= a
    assert a.tolist()==[[7,10],[15,22]]
    assert a2.tolist() == [[1,2],[3,4]]

    a = a2 = Matrix([[1,2], [3,4]])
    a *= a.A
    assert a.tolist()==[[1,4],[9,16]]
    assert a.data is a2.data

    a = a2 = Matrix([[1,2], [3,4]])
    hash(a)
    a *= a.A
    assert a.tolist()==[[1,4],[9,16]]
    assert a2.tolist()==[[1,2],[3,4]]
    assert a.data is not a2.data

    r = 2
    r *= a
    assert r==2*a

def test_mul():
    a = Matrix([[1,2], [3,4]])
    assert (a*a).tolist() == [[7,10],[15,22]],`(a*a).tolist()`
    assert (a.A*a).tolist() == [[1,4],[9,16]]
    assert (a*a.A).tolist() == [[1,4],[9,16]]
    assert (a*a.T).tolist() == [[5,11],[11,25]]
    assert (a.T*a).tolist() == [[10,14],[14,20]]
    assert (a.T*a.T).tolist() == [[7,15],[10,22]]
    assert (a.A*a.T).tolist() == [[1,6],[6,16]]
    assert (a*a.T.A).tolist() == [[1,6],[6,16]]

    assert (a*2).tolist() == [[2,4],[6,8]]
    assert (2*a).tolist() == [[2,4],[6,8]]

def test_div():
    a = a2 = Matrix([[1,2], [3,4]])*2
    assert (a/2).tolist() == [[1,2],[3,4]]
    assert (a/4).tolist() == [[mpq((1,2)),1],[mpq((3,2)),2]]

    a /= 2
    assert a.tolist() == [[1,2],[3,4]]
    assert a.data is a2.data

def test_rdiv():
    a = Matrix([[1,2],[3,4]])
    assert (1/a) == a.inv()
    assert (a/a).is_identity

    r = 1
    r /= a
    assert r==a.inv()
    assert (1/a.A).tolist()==[[1,mpq((1,2))],[mpq((1,3)),mpq((1,4))]]
    
def test_pow():
    a = Matrix([[1,2],[3,4]])
    assert a ** 0 == eye(2)
    assert a ** 1 == a
    assert a ** 2 == a*a
    assert a ** 3 == a*a*a
    assert a ** 4 == a*a*a*a
    assert a ** 5 == a*a*a*a*a
    assert a ** (-1) == a.inv()
    assert a ** (-2) == a.inv()*a.inv()
    assert a ** (-3) == a.inv()**3

    assert (a.A**1).tolist()==[[1,2],[3,4]]
    assert (a.A**2).tolist()==[[1,4],[9,16]]
    assert (a.A**-1).tolist()==[[1,mpq((1,2))],[mpq((1,3)),mpq((1,4))]]
    
def test_views():
    a = Matrix([[1,2], [3,4]])
    assert not a.head.is_array
    assert not a.M.head.is_array
    assert a.A.head.is_array
    b = a.A
    assert b.A is b

def test_trace():
    assert Matrix([[1,3],[6,9]]).trace() == 10
    assert Matrix([[1,3],[6,9]]).D.trace() == 10
    assert Matrix([[1,2,3],[4,5,6],[7,8,9]]).trace() == 15
    b = Matrix(10000, 10000)
    assert b.trace() == 0
    b[100,100] = 3
    b[1000,1000] = 4
    assert b.trace() == 7


def test_solve():
    assert Matrix([[1,2],[3,4]]).solve([1,2]).tolist()==[[0],[mpq((1,2))]]
    assert Matrix([[1,2],[3,4]]).solve([2,1]).tolist()==[[-3],[mpq((5,2))]]

    while 1:
        m = Matrix(3,3,random=True)
        if m.det():
            break

    b = Matrix(3,1,random=True)
    assert m * (m//b) == b
    
    b = Matrix(3,2,random=True)
    assert m * (m//b) == b

    b = Matrix(3,15,random=True)
    assert m * (m//b) == b

    while 1:
        m = Matrix(5,5,random=True)
        if m.det():
            break

    b = Matrix(5,1,random=True)
    assert m * (m//b) == b
    
    b = Matrix(5,2,random=True)
    assert m * (m//b) == b

    b = Matrix(5,15,random=True)
    assert m * (m//b) == b

def test_resize():

    m = Matrix([[1,2],[3,4]])
    assert m.resize(3,3).tolist()==[[1,2,0],[3,4,0],[0,0,0]]

    m = Matrix([[1,2],[3,4]])
    assert m.resize(2,3).tolist()==[[1,2,0],[3,4,0]]

    m = Matrix([[1,2],[3,4]])
    assert m.resize(2,1).tolist()==[[1],[3]]

    m = Matrix([[1,2],[3,4]])
    assert m.resize(2,1).resize(2,2).tolist()==[[1,2],[3,4]]

    m = Matrix([[1,2],[3,4]])
    assert m.resize(2,1).crop().resize(2,2).tolist()==[[1,0],[3,0]]
    assert m.tolist()==[[1,0],[3,0]]

def test_inv():
    m = Matrix([[1,2],[3,4]])
    assert (m*m.inv()).is_identity
    assert (m.inv()*m).is_identity

def test_lu_issue63():
    d = {(7, 3): -1, (6, 6): -1, (5, 6): 1, (2, 8): -1, (0, 3): 1, (1, 0): -1, (1, 2): -1, (4, 9): 1, (2, 9): 1, (1, 5): 1, (3, 0): 1, (7, 10): 1, (0, 4): -1, (4, 10): -1, (2, 6): 1, (5, 0): -1, (2, 10): 1, (3, 9): -1, (0, 5): -1, (6, 4): 1, (6, 1): -1, (5, 7): 1, (2, 4): 1}
    a = Matrix (8, 11, d)
    p,l,u = a.lu ()
    assert p*l*u == a
    assert (p*l*u - a).data == {}

def test_solve_null():
    x = ['x1', 'x2', 'x3', 'x4', 'x5', 'x6']
    a = Matrix ([[2,3,5],[-4,2,3]])
    xd,dep,indep = a.solve_null(x)
    assert xd['x1'] == Matrix([mpq((-1, 16))]),`xd`
    assert xd['x2'] == Matrix([mpq((-13, 8))]),`xd`
    assert xd['x3'] == Matrix([1]),`xd`
    assert set(dep) == set(['x1', 'x2']),`dep`
    assert indep == ['x3'],`dep`

    ker = Matrix([xd[s] for s in x[:3]])
    assert (a*ker).is_zero,`a*ker`

    a = Matrix(2,3,random=True)
    xd,dep,indep = a.solve_null(x)
    ker = Matrix([xd[s] for s in x[:3]])
    assert (a*ker).is_zero,`a*ker`

    a = Matrix ([[1,0,-3,0,2,-8],[0,1,5,0,-1,4],[0,0,0,1,7,-9],[0,0,0,0,0,0]])
    xd,dep,indep = a.solve_null(x)
    ker = Matrix([xd[s] for s in x[:6]])
    assert (a*ker).is_zero,`a*ker`

    a = Matrix(4,6,random=True)
    xd,dep,indep = a.solve_null(x)

    ker = Matrix([xd[s] for s in x[:6]])
    assert (a*ker).is_zero,`a*ker`

    a = Matrix([[1,0,0,0,-1,-1,-1,0,0,0,0],
                [0,1,0,0,1,-1,-1,0,0,0,0],
                [0,0,-1,0,0,1,0,1,0,0,0],
                [0,0,0,-1,0,0,1,-1,0,0,0],
                [0,0,0,0,0,0,0,0,1,-1,0],
                [0,0,0,0,0,0,0,0,-1,0,1]])
    xd,dep,indep = a.solve_null()
    ker = Matrix([xd[s] for s in sorted (dep+indep)])
    assert (a*ker).is_zero,`a*ker`

def test_gauss_jordan_elimination():
    a = Matrix([[1,0,0,0,-1,-1,-1,0,0,0,0],
                [0,1,0,0,1,-1,-1,0,0,0,0],
                [0,0,-1,0,0,1,0,1,0,0,0],
                [0,0,0,-1,0,0,1,-1,0,0,0],
                [0,0,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,1,-1,0],
                [0,0,0,0,0,0,0,0,-1,0,1]])
    r = Matrix([[1,0,0,0,-1,-1,-1,0,0,0,0],
                [0,1,0,0,1,-1,-1,0,0,0,0],
                [0,0,1,0,0,-1,0,-1,0,0,0],
                [0,0,0,1,0,0,-1,1,0,0,0],
                [0,0,0,0,0,0,0,0,1,0,-1],
                [0,0,0,0,0,0,0,0,0,1,-1]
                ])
    b = a.gauss_jordan_elimination()
    assert (b==r)

    a = a.T[:].T
    b = a.gauss_jordan_elimination()
    assert (b-r).is_zero,'\n%s==\n%s' % (b,r)

def test_gauss_jordan_elimination_swap_columns():
    a = Matrix([[1,0,0,0,-1,-1,-1,0,0,0,0],
                [0,1,0,0,1,-1,-1,0,0,0,0],
                [0,0,-1,0,0,1,0,1,0,0,0],
                [0,0,0,0,0,0,0,0,0,0,0],
                [0,0,0,-1,0,0,1,-1,0,0,0],
                [0,0,0,0,0,0,0,0,1,-1,0],
                [0,0,0,0,0,0,0,0,-1,0,1],
                ])
    r = Matrix([[1,0,0,0,-1,-1,-1,0,0,0,0],
                [0,1,0,0,1,-1,-1,0,0,0,0],
                [0,0,1,0,0,-1,0,-1,0,0,0],
                [0,0,0,1,0,0,-1,1,0,0,0],
                [0,0,0,0,0,0,0,0,1,0,-1],
                [0,0,0,0,0,0,0,0,0,1,-1]
                ])
    b,p = a.gauss_jordan_elimination(swap_columns=True)

    for k,j in enumerate (p):
        assert (b[:,k]==r[:,j])

    a = a.T[:].T
    b,p = a.gauss_jordan_elimination(swap_columns=True)
    for k,j in enumerate (p):
        assert (b[:,k]==r[:,j])

def test_get_gauss_jordan_elimination_operations():
    a = Matrix([[1,0,0,0,-1,-1,-1,0,0,0,0],
                [0,1,0,0,1,-1,-1,0,0,0,0],
                [0,0,-1,0,0,1,0,1,0,0,0],
                [0,0,-1,0,0,1,0,1,0,0,0],
                [0,0,0,-1,0,0,1,-1,0,0,0],
                [0,0,0,0,0,0,0,0,1,-1,0],
                [0,0,0,0,0,0,0,0,-1,0,1],
                [0,0,0,0,0,0,0,0,0,0,0],
                ]).T[:].T
    ab = Matrix ([[1, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 1, -1, -1, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 1, 0, 0, -1, 1, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1],
                  [0,0,0,0,0,0,0,0,0,0,0],
                  ])
    b, ops, r, c, z = a.get_gauss_jordan_elimination_operations()
    assert r==[0, 1, 2, 4, 5, 6],`r`
    assert z==[7, 3],`z`
    assert c==[0, 1, 2, 3, 8, 9]
    if not b[:]==ab:
        print
        print b
        print ab
    assert b[:]==ab
    a1 = a.apply_row_operations(ops)
    assert a1[:]==ab

def test_general_solver():
    a = Matrix([[1,0,0,0,-1,-1,-1,0,0,0,0],
                [0,1,0,0,1,-1,-1,0,0,0,0],
                [0,0,-1,0,0,1,0,1,0,0,0],
                [0,0,-1,0,0,1,0,1,0,0,0],
                [0,0,0,-1,0,0,1,-1,0,0,0],
                [0,0,0,0,0,0,0,0,1,-1,0],
                [0,0,0,0,0,0,0,0,-1,0,1],
                ])
    cols = ['A_D','A_B','AB_C','B_D','C_DE','Eout','Ain','Dout','Cin']
    a = Matrix([[-1,-1,-1,0 ,0 ,0, 1, 0,0],    # A
                [0 ,1 ,-1,-1,0 ,0, 0, 0,0],    # B
                [0 ,0 ,1 ,0 ,-1,0, 0, 0,1],    # C
                [1 ,0 ,0 ,1 ,1 ,0, 0,-1,0],    # D
                [0 ,0 ,0 ,0 ,1 ,-1,0, 0,0],    # E
                [0 ,0 ,0 ,0 ,0 ,0,-1, 0,0],    # Asrc
                [0 ,0 ,0 ,0 ,0 ,0, 0, 1,0],    # Dt
                [0 ,0 ,0 ,0 ,0 ,0, 0, 0,-1],    # Csrc
                ])
    print
    print 'a'
    print a

    for i,row in enumerate(a):
        break
        if len(row.data)<=1:
            a[i] = 0

    print "a'"
    print a
    
    b, ops, r, c, z = a.get_gauss_jordan_elimination_operations(leading_cols=range(5),
                                                                leading_column_selection='sparsest first'
                                                                )
    print b
    print r
    print c
    print z
