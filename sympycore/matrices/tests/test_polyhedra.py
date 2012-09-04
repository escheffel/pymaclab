
import time
from sympycore import Polyhedron, Matrix, mpq


def test_LP_solve_Chateau_ETH_Production():
    constraints = '''\
3*x1+4*x2+2*x3==MAX
2*x1<=4
x1+2*x3<=8
3*x2+x3<=6
'''
    p = Polyhedron(*constraints.split('\n'))
    #p.show()
    names, D = p.get_LP()
    Dwork = D[:]
    xopt, vopt = Dwork.LP_solve(overwrite=True)
    assert (D[0,1:] * xopt)[0,0]==vopt==16,`D[0,1:] * xopt,vopt`

def test_LP_solve_cyclic():
    constraints = '''
x1-2*x2+x3==MAX
2*x1-x2+x3<=0
3*x1+x2+x3<=0
-5*x1+3*x2-2*x3<=0
'''
    p = Polyhedron(*constraints.split('\n'))
    #p.show()
    names, D = p.get_LP()
    Dwork = D[:]
    xopt, vopt = Dwork.LP_solve(overwrite=True)
    assert (D[0,1:] * xopt)[0,0]==vopt==0,`D[0,1:] * xopt,vopt`

def test_LP_solve_Karmarkar_example ():
    p = Polyhedron('x1+x2==MAX')
    for i in range (11):
        p1 = mpq((i,10))*1
        p.add('2*%s*x1+x2<=%s+1' % (p1,p1**2))
    names, D = p.get_LP()
    Dwork = D[:]
    xopt, vopt = Dwork.LP_solve(overwrite=True)
    assert (D[0,1:] * xopt)[0,0]==vopt==mpq ((5,4)),`D[0,1:] * xopt,vopt`

def test_Av98a():
    constraints = '''\
#x1+x2+x3==1
1+x1>=0
1+x2>=0
1-x1>=0
1-x2>=0
1-x1+x3>=0
1-x2+x3>=0
1+x1+x3>=0
1+x2+x3>=0'''
    p = Polyhedron(*constraints.split('\n'))
    #p.show()

    return

    vertices = [(1,1,0),(-1,1,0),(1,-1,0),(-1,-1,0),(0,0,-1)]
    rays = [(0,0,1)]

    print 'A_L'
    print p.A_L

    xd, dep, indep = p.A_L.solve_null(labels=p.names)
    Aker = Matrix([xd[l].tolist()[0] for l in p.names])
    print 'Aker:'
    print Aker

    print 'A * Aker:'
    print p.A * Aker

    xd, dep, indep = p.A_I[(0,1),:].solve_null(labels=p.names)
    Aker = Matrix([xd[l].tolist()[0] for l in p.names])
    print 'Aker:'
    print Aker
    print 'A * Aker:'
    print p.A * Aker

    print p.A[(0,1,2,4),:].I*2
    print p.A[(0,1,2,4),:]*p.A[(0,1,2,4),:].I
