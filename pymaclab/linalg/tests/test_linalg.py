import os
from scipy.io import loadmat
from numpy.testing import assert_almost_equal
from pymaclab.linalg import qz, ordqz

current_path = os.path.dirname(os.path.abspath(__file__))

def test_qz():
    res = loadmat(current_path+'/results/qz.mat')
    A,B = res['A'],res['B']
    AA,BB,Q,Z = qz(A,B)
    assert_almost_equal(AA,res['AA'],8)
    assert_almost_equal(BB,res['BB'],8)
    assert_almost_equal(Q,res['Q'],8)
    assert_almost_equal(Z,res['Z'],8)

def test_qz_real():
    res = loadmat(current_path+'/results/qzreal.mat')
    A,B = res['A'],res['B']
    AA,BB,Q,Z = qz(A,B,mode='real')
    assert_almost_equal(AA,res['AA'],8)
    assert_almost_equal(BB,res['BB'],8)
    assert_almost_equal(Q,res['Q'],8)
    assert_almost_equal(Z,res['Z'],8)

def test_ordqz():
    res = loadmat(current_path+'/results/ordqz.mat')
    AA,BB,Q,Z,select = res['AA'],res['BB'],res['Q'],res['Z'],res['selector']
    AAS,BBS,QS,ZS = ordqz(AA,BB,Q,Z,select)
    assert_almost_equal(AAS,res['AAS'],8)
    assert_almost_equal(BBS,res['BBS'],8)
    assert_almost_equal(QS,res['QS'],8)
    assert_almost_equal(ZS,res['ZS'],8)

if __name__ == "__main__":
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb'], exit=False)
