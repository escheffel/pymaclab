'''
.. module:: modsolvers
   :platform: Linux
   :synopsis: The modsolvers module is an integral part of PyMacLab and is called extensively from the DSGEmodel class in module
              macrolab. It contains all of the dynamic solver methods most of which make use of the DSGE models' derivatives, i.e.
              the computed (numerical) Jacobian and Hessian.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''

from numpy import matlib as MAT
from numpy import linalg as LIN
from numpy.linalg import matrix_rank
from .. import helpers as HLP
from scipy.linalg import qz
import numpy as np
import pylab as P
from matplotlib import pyplot as PLT
import copy as COP
from pymaclab.filters import hpfilter
from pymaclab.filters import bkfilter
from pymaclab.filters import cffilter
from isolab import isolab
try:
    import mlabraw
except:
    pass


class MODsolvers(object):
    def __init__(self):
        pass


class PyUhlig(MODsolvers):
    def __init__(self,intup):
        self.ntup = intup[0]
        self.nendo = self.ntup[0]
        self.ncon = self.ntup[1]
        self.nexo = self.ntup[2]

        self.eqindx = intup[1]
        self.vreg = intup[2]
        self.llsys_list = intup[3]
        self.diffli1 = intup[4]
        self.diffli2 = intup[5]
        diffli1 = self.diffli1
        diffli2 = self.diffli2

        deteq = []
        for x in self.eqindx['det']:
            deteq.append(self.llsys_list[x])
        self.deteq = deteq
        expeq = []
        for x in self.eqindx['exp']:
            expeq.append(self.llsys_list[x])
        self.expeq = expeq
        erreq = []
        for x in self.eqindx['err']:
            erreq.append(self.llsys_list[x])
        self.erreq = erreq

        detdif1 = []
        detdif2 = []
        for x in self.eqindx['det']:
            detdif1.append(diffli1[x])
            detdif2.append(diffli2[x])
        expdif1 = []
        expdif2 = []
        for x in self.eqindx['exp']:
            expdif1.append(diffli1[x])
            expdif2.append(diffli2[x])
        errdif1 = []
        errdif2 = []
        for x in self.eqindx['err']:
            errdif1.append(diffli1[x])
            errdif2.append(diffli2[x])
        self.detdif1 = detdif1
        self.detdif2 = detdif2
        self.expdif1 = expdif1
        self.expdif2 = expdif2
        self.errdif1 = errdif1
        self.errdif2 = errdif2

        self.mkmats()

    def mkmats(self):
        deteq = self.deteq
        expeq = self.expeq
        erreq = self.erreq
        nendo = self.nendo
        ncon = self.ncon
        nexo = self.nexo        
        vreg = self.vreg
        detdif1 = self.detdif1
        detdif2 = self.detdif2
        expdif1 = self.expdif1
        expdif2 = self.expdif2
        errdif1 = self.errdif1
        errdif2 = self.errdif2

        # Create argument list for matrix creation
        argli = [['AA',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','0')],
             ['BB',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','-1')],
             ['CC',(detdif2,detdif1),deteq,(len(deteq),ncon),(None,'con','0')],
             ['DD',(detdif2,detdif1),deteq,(len(deteq),nexo),(None,'exo','0')],
             ['FF',(expdif2,expdif1),expeq,(len(expeq),nendo),('0','endo','1')],
             ['GG',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','0')],
             ['HH',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','-1')],
             ['JJ',(expdif2,expdif1),expeq,(len(expeq),ncon),('0','con','1')],
             ['KK',(expdif2,expdif1),expeq,(len(expeq),ncon),(None,'con','0')],
             ['LL',(expdif2,expdif1),expeq,(len(expeq),nexo),('0','exo','1')],
             ['MM',(expdif2,expdif1),expeq,(len(expeq),nexo),(None,'exo','0')],
             ['NN',(errdif2,errdif1),erreq,(len(erreq),nexo),(None,'exo','0')]]

        # Create all matrices in a loop over argli
        for argx in argli:
            XX = MAT.zeros(argx[3])
            dic1 = dict([[x,'0'] for x in range(argx[3][1])])
            sXX = dict([[x,COP.deepcopy(dic1)] for x in range(len(argx[2]))])
            for y,x in enumerate(argx[2]):
                if vreg(argx[4],x,False,'max'):
                    cont=[[z[0],z[1][1]] for z in vreg(argx[4],x,True,'min')]
                    for z in cont:
                        XX[y,z[1]] = argx[1][0][y][z[0]]
                        sXX[y][z[1]] = argx[1][1][y][z[0]]
            exec('self.'+argx[0]+' = XX')
            exec('self.'+'s'+argx[0]+' = sXX')

    def solve(self,sol_method='do_QZ',Xi_select='all'):
        self.output = {}
        self.Xi_select = Xi_select
        self.sol_method = sol_method
        if self.sol_method == 'do_PP':
            self.do_PP(self.Xi_select)
        elif self.sol_method == 'do_QZ':
            self.do_QZ(self.Xi_select,Tol=1.0000e-06)
        self.do_QRS()

    def do_QZ(self,Xi_select='all',Tol=1.0000e-06):
        # Make uhlig matrices locally available for computations
        AA = self.AA
        BB = self.BB
        CC = self.CC
        DD = self.DD
        FF = self.FF
        GG = self.GG
        HH = self.HH
        JJ = self.JJ
        KK = self.KK
        LL = self.LL
        MM = self.MM
        NN = self.NN
        q_expectational_equ = np.shape(FF)[0]
        m_states = np.shape(FF)[1]
        l_equ = np.shape(CC)[0]
        n_endog = np.shape(CC)[1]
        k_exog = min(np.shape(NN))

#        if HLP.rank(CC) < n_endog:
        if matrix_rank(CC) < n_endog:
            raise uhlerr, 'uhlerror: Rank(CC) needs to be at least\n\
                  equal to the number of endogenous variables'
        if l_equ == 0:
            CC_plus = LIN.pinv(CC)
            CC_0 = (HLP.null(CC.T)).T
            Psi_mat = FF
            Gamma_mat = -GG
            Theta_mat = -HH
            Xi_mat = \
                   np.concatenate((Gamma_mat,Theta_mat,MAT.eye(m_states),\
                          MAT.zeros((m_states,m_states))),1)
            Delta_mat = \
                  np.concatenate((Psi_mat,MAT.zeros((m_states,m_states)),\
                         MAT.zeros((m_states,m_states)),\
                         MAT.zeros((m_states,m_states)),MAT.eye(m_states),1))


        else:
            CC_plus = LIN.pinv(CC)
            CC_0 = HLP.null(CC.T)
            if l_equ > n_endog:
                Psi_mat = \
                    np.concatenate((MAT.zeros((l_equ-n_endog,m_states)),FF-JJ*CC_plus*AA),1)
            elif l_equ == n_endog:
                Psi_mat = np.concatenate((FF-JJ*CC_plus*AA),1)
            Gamma_mat = np.concatenate((CC_0*AA, JJ*CC_plus*BB-GG+KK*CC_plus*AA))
            Theta_mat = np.concatenate((CC_0*BB, KK*CC_plus*BB-HH))
            Xi_mat = \
                   np.concatenate((\
                       np.concatenate((Gamma_mat,Theta_mat),1),\
                       np.concatenate((MAT.identity(m_states),MAT.zeros((m_states,m_states))),1)\
                   ))
            Psi_mat = Psi_mat.reshape(m_states,m_states)
            Delta_mat = \
                  np.concatenate((\
                      np.concatenate((Psi_mat,MAT.zeros((m_states,m_states))),1),\
                      np.concatenate((MAT.zeros((m_states,m_states)),MAT.identity(m_states)),1)\
                  ))

#        (Delta_up,Xi_up,UUU,VVV)=\
#         HLP.qz(Delta_mat,Xi_mat,\
#            mode='complex',\
#            lapackname=lapackname,\
#            lapackpath=lapackpath)
        (Delta_up,Xi_up,UUU,VVV) = qz(Delta_mat,Xi_mat)
        

        d_Delta_up = MAT.diag(Delta_up)
        d_Xi_up = MAT.diag(Xi_up)
        Xi_eigval = MAT.zeros(N.shape(Delta_up))
        for i1 in range(0,N.shape(Delta_up)[0],1):
            Xi_eigval[i1,i1] = d_Xi_up[i1]/d_Delta_up[i1]
        d_Xi_eigval = np.diag(Xi_eigval)
        mat_tmp = MAT.zeros((N.shape(Xi_eigval)[0],3))
        i1=0
        for x in d_Xi_eigval:
            mat_tmp[i1,0] = d_Xi_eigval[i1]
            mat_tmp[i1,1] = abs(d_Xi_eigval)[i1]
            mat_tmp[i1,2] = i1
            i1=i1+1

        Xi_sortval = HLP.sortrows(mat_tmp,1)[:,0]
        # Need to do an argsort() on index column to turn to integers (Booleans) !!
        Xi_sortindex = HLP.sortrows(mat_tmp,1)[:,2].argsort(axis=0)
        Xi_sortabs = HLP.sortrows(mat_tmp,1)[:,1]

        # Root selection branch with Xi_select
        if Xi_select == 'all':
            Xi_select = np.arange(0,m_states,1)

        stake = max(abs(Xi_sortval[Xi_select])) + Tol
        (Delta_up,Xi_up,UUU,VVV) = self.qzdiv(stake,Delta_up,Xi_up,UUU,VVV)
        Lamda_mat = np.diag(Xi_sortval[Xi_select])
        trVVV = VVV.conjugate().T
        VVV_2_1 = trVVV[m_states:2*m_states,0:m_states]
        VVV_2_2 = trVVV[m_states:2*m_states,m_states:2*m_states]
        UUU_2_1 = UUU[m_states:2*m_states,0:m_states]

        PP = -(VVV_2_1.I*VVV_2_2)
        PP = np.real(PP)

        self.PP = PP
        self.CC_plus = CC_plus

    def do_PP(self,Xi_select='all'):
        # Make uhlig matrices locally available for computations
        AA = self.AA
        BB = self.BB
        CC = self.CC
        DD = self.DD
        FF = self.FF
        GG = self.GG
        HH = self.HH
        JJ = self.JJ
        KK = self.KK
        LL = self.LL
        MM = self.MM
        NN = self.NN
        q_expectational_equ = np.shape(FF)[0]
        m_states = np.shape(FF)[1]
        l_equ = np.shape(CC)[0]
        n_endog = np.shape(CC)[1]
        k_exog = min(N.shape(NN))

#        if HLP.rank(CC) < n_endog:
        if matrix_rank(CC) < n_endog:
            raise uhlerr, 'uhlerror: Rank(CC) needs to be at least\n\
                  equal to the number of endogenous variables'
        if l_equ == 0:
            CC_plus = LIN.pinv(CC)
            CC_0 = (HLP.null(CC.T)).T
            Psi_mat = FF
            Gamma_mat = -GG
            Theta_mat = -HH
            Xi_mat = \
                   np.concatenate((Gamma_mat,Theta_mat,MAT.eye(m_states),\
                          MAT.zeros((m_states,m_states))),1)
            Delta_mat = \
                  np.concatenate((Psi_mat,MAT.zeros((m_states,m_states)),\
                         MAT.zeros((m_states,m_states)),\
                         MAT.zeros((m_states,m_states)),MAT.eye(m_states),1))


        else:
            CC_plus = LIN.pinv(CC)
            CC_0 = HLP.null(CC.T)
            if l_equ > n_endog:
                Psi_mat = \
                    np.concatenate((MAT.zeros((l_equ-n_endog,m_states)),FF-JJ*CC_plus*AA),1)
            elif l_equ == n_endog:
                Psi_mat = np.concatenate((FF-JJ*CC_plus*AA),1)
            Gamma_mat = np.concatenate((CC_0*AA, JJ*CC_plus*BB-GG+KK*CC_plus*AA))
            Theta_mat = np.concatenate((CC_0*BB, KK*CC_plus*BB-HH))
            Xi_mat = \
                   np.concatenate((\
                       np.concatenate((Gamma_mat,Theta_mat),1),\
                       np.concatenate((MAT.eye(m_states),MAT.zeros((m_states,m_states))),1)\
                   ))
            Delta_mat = \
                  np.concatenate((\
                      np.concatenate((Psi_mat,MAT.zeros((m_states,m_states))),1),\
                      np.concatenate((MAT.zeros((m_states,m_states)),MAT.eye(m_states)),1)\
                  ))

        (Xi_eigvec,Xi_eigval) = HLP.eig(Xi_mat,Delta_mat)
        tmp_mat = MAT.eye(N.rank(Xi_eigvec))
        for i1 in range(0,N.rank(Xi_eigvec),1):
            tmp_mat[i1,i1] = float(Xi_eigval[0,i1])
        Xi_eigval = tmp_mat


#        if HLP.rank(Xi_eigvec) < m_states:
        if matrix_rank(Xi_eigvec) < m_states:
            raise uhlerr, 'uhlerror: Xi_mat is not diagonalizable!\n\
                  Cannot solve for PP. Maybe you should try the Schur-Decomposition\n\
                  method instead, use do_QZ()!!'


        d_Xi_eigval = np.diag(Xi_eigval)
        mat_tmp = MAT.zeros((N.rank(Xi_eigval),3))
        i1=0
        for x in d_Xi_eigval:
            mat_tmp[i1,0] = d_Xi_eigval[i1]
            mat_tmp[i1,1] = abs(d_Xi_eigval[i1])
            mat_tmp[i1,2] = i1
            i1=i1+1
        Xi_sortval = HLP.sortrows(mat_tmp,1)[:,0]
        # Need to do an argsort() on index column to turn to integers (Booleans) !!
        Xi_sortindex = HLP.sortrows(mat_tmp,1)[:,2].argsort(axis=0)
        Xi_sortabs = HLP.sortrows(mat_tmp,1)[:,1]
        Xi_sortvec = Xi_eigvec[0:2*m_states,:].take(Xi_sortindex.T.A, axis=1)

        # Root selection branch with Xi_select
        if Xi_select == 'all':
            Xi_select = np.arange(0,m_states,1)

        #Create Lambda_mat and Omega_mat
        Lambda_mat = MAT.zeros((len(Xi_select),len(Xi_select)))
        for i1 in range(0,len(Xi_select),1):
            Lambda_mat[i1,i1] = Xi_sortval[Xi_select][i1]
        Omega_mat = Xi_sortvec[m_states:2*m_states,:].take(Xi_select, axis=1)

#        if HLP.rank(Omega_mat) < m_states:
        if matrix_rank(Omega_mat) < m_states:
            raise uhlerr, 'uhlerror: Omega_mat is not invertible!!\n\
                  Therefore, solution for PP is not available'

        PP = Omega_mat*HLP.ctr((HLP.ctr(Omega_mat).I*HLP.ctr(Lambda_mat)))
        self.PP = PP
        self.CC_plus = CC_plus

    def do_QRS(self):
        # Make uhlig matrices locally available for computations
        AA = self.AA
        BB = self.BB
        CC = self.CC
        DD = self.DD
        FF = self.FF
        GG = self.GG
        HH = self.HH
        JJ = self.JJ
        KK = self.KK
        LL = self.LL
        MM = self.MM
        NN = self.NN
        PP = self.PP
        CC_plus = self.CC_plus
        (l_equ,m_states) = MAT.shape(AA)
        (l_equ,n_endog) = MAT.shape(CC)
        (l_equ,k_exog) = MAT.shape(DD)
        if l_equ == 0:
            RR = MAT.zeros((0,m_states))
            VV1 = MAT.kron(MAT.transpose(NN),FF)+MAT.kron(MAT.identity(k_exog),FF*PP+GG)
            VV2 = MAT.kron(MAT.transpose(NN),JJ)+MAT.kron(MAT.identity(k_exog),KK)
            VV = MAT.hstack((VV1,VV2))
        else:
            RR = -CC_plus*(AA*PP+BB)
            VV1 = MAT.kron(MAT.identity(k_exog),AA)
            VV2 = MAT.kron(MAT.identity(k_exog),CC)
            VV3 = MAT.kron(MAT.transpose(NN),FF)+MAT.kron(MAT.identity(k_exog),FF*PP+JJ*RR+GG)
            VV4 = MAT.kron(MAT.transpose(NN),JJ)+MAT.kron(MAT.identity(k_exog),KK)
            VV = MAT.vstack((MAT.hstack((VV1,VV2)),MAT.hstack((VV3,VV4))))
        self.RR = RR

        LLNN_plus_MM = LL*NN+MM
        NON = MAT.hstack((DD.flatten(1),LLNN_plus_MM.flatten(1)))
        try:
            QQSS_vec = -(VV.I*MAT.transpose(NON))
        except MyErr:
            print 'Error: Matrix VV is not invertible!'
        QQ = QQSS_vec[0:m_states*k_exog,:].reshape(-1,m_states).transpose()
        SS = QQSS_vec[(m_states*k_exog):((m_states+n_endog)*k_exog),:].reshape(-1,n_endog).transpose()
        WW1 = MAT.hstack((MAT.identity(m_states),MAT.zeros((m_states,k_exog))))
        WW2 = MAT.hstack((RR*LIN.pinv(PP),SS-RR*LIN.pinv(PP)*QQ))
        WW3 = MAT.hstack((MAT.zeros((k_exog,m_states)),MAT.identity(k_exog)))
        WW = MAT.vstack((WW1,WW2,WW3))
        self.WW = WW
        self.QQ = QQ
        self.SS = SS
        del self.CC_plus

    def qzdiv(self,stake,A,B,Q,Z):
        n = np.shape(A)[0]
        root = np.hstack((abs(N.mat(N.diag(A)).T),abs(N.mat(N.diag(B)).T)))
        index_mat = (root[:,0]<1e-13).choose(root[:,0],1)
        index_mat = (index_mat>1e-13).choose(root[:,0],0)
        root[:,0] = root[:,0]-MAT.multiply(index_mat,(root[:,0]+root[:,1]))
        root[:,1] = root[:,1]/root[:,0]
        for i1 in range(n-1,-1,-1):
            m='none'
            for i2 in range(i1,-1,-1):
                if root[i2,1] > stake or root[i2,1] < -0.1:
                    m=i2
                    break
            if m == 'none':
                return (A,B,Q,Z)
            else:
                for i3 in range(m,i1,1):
                    (A,B,Q,Z) = self.qzswitch(i3,A,B,Q,Z)
                    tmp = COP.deepcopy(root[i3,1])
                    root[i3,1] = root[i3+1,1]
                    root[i3+1,1] = tmp

    def qzswitch(self,i,A,B,Q,Z):
        a = A[i,i]
        d = B[i,i]
        b = A[i,i+1]
        e = B[i,i+1]
        c = A[i+1,i+1]
        f = B[i+1,i+1]
        wz = np.mat(N.hstack(((c*e-f*b),(c*d-f*a).conjugate().T)))
        xy = np.mat(N.hstack(((b*d-e*a).conjugate().T,(c*d-f*a).conjugate().T)))
        n = np.mat(N.sqrt(wz*wz.conjugate().T))
        m = np.mat(N.sqrt(xy*xy.conjugate().T))
        if n.all() == 0:
            return
        else:
            wz = np.mat(wz/n)
            xy = np.mat(xy/m)
            wz = np.vstack((wz,N.hstack((-wz[:,1].T,wz[:,0].T))))
            xy = np.vstack((xy,N.hstack((-xy[:,1].T,xy[:,0].T))))
            A[i:i+2,:] = xy*A[i:i+2,:]
            B[i:i+2,:] = xy*B[i:i+2,:]
            A[:,i:i+2] = A[:,i:i+2]*wz
            B[:,i:i+2] = B[:,i:i+2]*wz
            Z[:,i:i+2] = Z[:,i:i+2]*wz
            Q[i:i+2,:] = xy*Q[i:i+2,:]
        return (A,B,Q,Z)

#----------------------------------------------------------------------------------------------------------------------
class MatUhlig:

    def __init__(self,intup):
        self.ntup = intup[0]
        self.nendo = self.ntup[0]
        self.ncon = self.ntup[1]
        self.nexo = self.ntup[2]

        self.eqindx = intup[1]
        self.vreg = intup[2]
        self.llsys_list = intup[3]
        self.diffli1 = intup[4]
        self.diffli2 = intup[5]
        self.sess1 = intup[6]
        self.vardic = intup[7]
        diffli1 = self.diffli1
        diffli2 = self.diffli2

        deteq = []
        for x in self.eqindx['det']:
            deteq.append(self.llsys_list[x])
        self.deteq = deteq
        expeq = []
        for x in self.eqindx['exp']:
            expeq.append(self.llsys_list[x])
        self.expeq = expeq
        erreq = []
        for x in self.eqindx['err']:
            erreq.append(self.llsys_list[x])
        self.erreq = erreq

        detdif1 = []
        detdif2 = []
        for x in self.eqindx['det']:
            detdif1.append(diffli1[x])
            detdif2.append(diffli2[x])
        expdif1 = []
        expdif2 = []
        for x in self.eqindx['exp']:
            expdif1.append(diffli1[x])
            expdif2.append(diffli2[x])
        errdif1 = []
        errdif2 = []
        for x in self.eqindx['err']:
            errdif1.append(diffli1[x])
            errdif2.append(diffli2[x])
        self.detdif1 = detdif1
        self.detdif2 = detdif2
        self.expdif1 = expdif1
        self.expdif2 = expdif2
        self.errdif1 = errdif1
        self.errdif2 = errdif2

        self.mkmats()

    def mkmats(self):
        deteq = self.deteq
        expeq = self.expeq
        erreq = self.erreq
        nendo = self.nendo
        ncon = self.ncon
        nexo = self.nexo        
        vreg = self.vreg
        detdif1 = self.detdif1
        detdif2 = self.detdif2
        expdif1 = self.expdif1
        expdif2 = self.expdif2
        errdif1 = self.errdif1
        errdif2 = self.errdif2

        # Create argument list for matrix creation
        argli = [['AA',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','0')],
             ['BB',(detdif2,detdif1),deteq,(len(deteq),nendo),(None,'endo','-1')],
             ['CC',(detdif2,detdif1),deteq,(len(deteq),ncon),(None,'con','0')],
             ['DD',(detdif2,detdif1),deteq,(len(deteq),nexo),(None,'exo','0')],
             ['FF',(expdif2,expdif1),expeq,(len(expeq),nendo),('0','endo','1')],
             ['GG',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','0')],
             ['HH',(expdif2,expdif1),expeq,(len(expeq),nendo),(None,'endo','-1')],
             ['JJ',(expdif2,expdif1),expeq,(len(expeq),ncon),('0','con','1')],
             ['KK',(expdif2,expdif1),expeq,(len(expeq),ncon),(None,'con','0')],
             ['LL',(expdif2,expdif1),expeq,(len(expeq),nexo),('0','exo','1')],
             ['MM',(expdif2,expdif1),expeq,(len(expeq),nexo),(None,'exo','0')],
             ['NN',(errdif2,errdif1),erreq,(len(erreq),nexo),(None,'exo','0')]]

        # Create all matrices in a loop over argli
        for argx in argli:
            XX = MAT.zeros(argx[3])
            dic1 = dict([[x,'0'] for x in range(argx[3][1])])
            sXX = dict([[x,COP.deepcopy(dic1)] for x in range(len(argx[2]))])
            for y,x in enumerate(argx[2]):
                if vreg(argx[4],x,False,'max'):
                    cont=[[z[0],z[1][1]] for z in vreg(argx[4],x,True,'min')]
                    for z in cont:
                        XX[y,z[1]] = argx[1][0][y][z[0]]
                        sXX[y][z[1]] = argx[1][1][y][z[0]]
            exec('self.'+argx[0]+' = XX')
            exec('self.'+'s'+argx[0]+' = sXX')

    def solve(self,sol_method='do_QZ',Xi_select='all'):
        # Create varnames with correct string length
        tmp_list =[x[1] for x in self.vardic['endo']['var']]\
             +[x[1] for x in self.vardic['con']['var']]\
             +[x[1] for x in self.vardic['exo']['var']]
        tmp_list2 = [[len(x),x] for x in tmp_list]
        tmp_list2.sort()
        max_len = tmp_list2[-1][0]
        i1=0
        for x in tmp_list:
            tmp_list[i1]=x+(max_len-len(x))*' '
            i1=i1+1
        varnstring = 'VARNAMES = ['
        for x in tmp_list:
            varnstring = varnstring + "'" +x[1]+ "'" + ','
        varnstring = varnstring[0:-1]+'];'

        # Start matlab session and calculate
        sess1 = self.sess1
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,varnstring)
        mlabraw.put(sess1,'AA',self.AA)
        mlabraw.put(sess1,'BB',self.BB)
        mlabraw.put(sess1,'CC',self.CC)
        mlabraw.put(sess1,'DD',self.DD)
        mlabraw.put(sess1,'FF',self.FF)
        mlabraw.put(sess1,'GG',self.GG)
        mlabraw.put(sess1,'HH',self.HH)
        mlabraw.put(sess1,'JJ',self.JJ)
        mlabraw.put(sess1,'KK',self.KK)
        mlabraw.put(sess1,'LL',self.LL)
        mlabraw.put(sess1,'MM',self.MM)
        mlabraw.put(sess1,'NN',self.NN)
        mlabraw.eval(sess1,'cd '+mlabpath)
        mlabraw.eval(sess1,'cd Toolkit41')
        message = ' '*70
        eval_list = ['message='+"'"+message+"'",
                 'warnings = [];',
                 'options;',
                 'solve;']
        try:
            for x in eval_list:
                mlabraw.eval(sess1,x)
            self.PP = np.matrix(mlabraw.get(sess1,'PP'))
            self.QQ = np.matrix(mlabraw.get(sess1,'QQ'))
            self.RR = np.matrix(mlabraw.get(sess1,'RR'))
            self.SS = np.matrix(mlabraw.get(sess1,'SS'))
            self.WW = np.matrix(mlabraw.get(sess1,'WW'))
        except MyErr:
            pass
        finally:
            return
#----------------------------------------------------------------------------------------------------------------------
class MatKlein:

    def __init__(self,intup):
        self.sess1 = intup[-1]
        self.uhlig = PyUhlig(intup[:-1])
        self.uhlig.mkmats()
        self.AA = self.uhlig.AA
        self.BB = self.uhlig.BB
        self.CC = self.uhlig.CC
        self.DD = self.uhlig.DD
        self.FF = self.uhlig.FF
        self.GG = self.uhlig.GG
        self.HH = self.uhlig.HH
        self.JJ = self.uhlig.JJ
        self.KK = self.uhlig.KK
        self.LL = self.uhlig.LL
        self.MM = self.uhlig.MM
        self.NN = self.uhlig.NN
        self.mkKleinMats()
        self.outdic = {}

    def mkKleinMats(self):
        # Make uhlig matrices locally available for computations
        AA = self.AA
        BB = self.BB
        CC = self.CC
        DD = self.DD
        FF = self.FF
        GG = self.GG
        HH = self.HH
        JJ = self.JJ
        KK = self.KK
        LL = self.LL
        MM = self.MM
        NN = self.NN

        # Determine size of states, endogenous
        exo_st = MAT.shape(NN)[1]
        endo_st = MAT.shape(BB)[1]
        endo_cn = MAT.shape(CC)[1]
        n_deteq = MAT.shape(AA)[0]
        n_expeq = MAT.shape(JJ)[0]
        tot_st = exo_st+endo_st
        self.tstates = tot_st
        tot_var = tot_st+endo_cn

        klein_A_rtwo = MAT.hstack((LL,GG))
        klein_A_rtwo = MAT.hstack((klein_A_rtwo,JJ))
        klein_A_rtwo_rows = MAT.shape(klein_A_rtwo)[0]
        klein_A_rtwo_cols = MAT.shape(klein_A_rtwo)[1]

        klein_A_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))
        klein_A_rone = MAT.hstack((MAT.identity(exo_st),klein_A_rone[:,exo_st:]))

        klein_A_rthree = MAT.hstack((MAT.zeros((n_deteq,exo_st)),AA))
        klein_A_rthree = MAT.hstack((klein_A_rthree,MAT.zeros((n_deteq,endo_cn))))

        klein_A = MAT.vstack((klein_A_rone,klein_A_rtwo))
        klein_A = MAT.vstack((klein_A,klein_A_rthree))

        klein_B_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))

        klein_B_rone = MAT.hstack((NN,klein_B_rone[:,exo_st:]))
        klein_B_rtwo = MAT.hstack((-MM,-HH))
        klein_B_rtwo = MAT.hstack((klein_B_rtwo,-KK))
        klein_B_rthree = MAT.hstack((-DD,-BB))
        klein_B_rthree = MAT.hstack((klein_B_rthree,-CC))

        klein_B = MAT.vstack((klein_B_rone,klein_B_rtwo))
        klein_B = MAT.vstack((klein_B,klein_B_rthree))  

        self.A = klein_A
        self.B = klein_B

    def solve(self):
        A = self.A
        B = self.B
        tstates = self.tstates
        sess1 = self.sess1
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,'cd '+mlabpath)
        mlabraw.eval(sess1,'cd Klein')
        mlabraw.put(sess1,'AA',A)
        mlabraw.put(sess1,'BB',B)
        mlabraw.put(sess1,'tstates',tstates)
        try:
            mlabraw.eval(sess1,'[F,P,Z11]=solab(AA,BB,tstates)')
            self.F = np.matrix(mlabraw.get(sess1,'F'))
            self.P = np.matrix(mlabraw.get(sess1,'P'))
            self.Z11 = np.matrix(mlabraw.get(sess1,'Z11'))
            self.outdic['P'] = self.P
            self.outdic['F'] = self.F
        except MyErr:
            pass
        finally:
            return
#----------------------------------------------------------------------------------------------------------------------
class MatKleinD:

    def __init__(self,intup):
        self.interm3 = intup[0]
        self.param = intup[1]
        self.sstate_list = intup[2]
        self.varvecs = intup[3]
        self.vardic = intup[4]
        self.vardic2 = intup[5]
        self.shockdic = intup[6]
        self.shockdic2 = intup[7]
        self.varns = intup[8]
        self.re_var = intup[9]
        self.re_var2 = intup[10]
        self.ishockdic = intup[11]
        self.ishockdic2 = intup[12]
        self.sess1 = intup[13]
        self.sstate = intup[14]
        self.mkeqs()
        self.mkvarl()
        self.mkeqs2()
        self.mkgrad()
        self.mkAB()
        self.solab()
        self.mkhess()
        self.solab2()

    def mkeqs(self):
        list_tmp = COP.deepcopy(self.interm3)
        reg2 = RE.compile('\*{2,2}')
        reva2 = self.re_var2
        i1=0
        for x in list_tmp:
            list_tmp[i1]=reg2.sub('^',x.split('=')[0]).strip()
            for y in self.subvars(tuple(self.varvecs['e'])):
                while y in list_tmp[i1]:
                    list_tmp[i1] = list_tmp[i1].replace(y,'1')
            i1=i1+1
        self.interm4 = list_tmp

    def mkvarl(self):
        varlist=[]
        nexost = len(self.varvecs['z'])
        nendost = len(self.varvecs['k'])
        nendocon = len(self.varvecs['y'])
        exost_f = makeForward(self.varvecs['z'],'1','0')
        for x in exost_f:
            varlist.append(x)
        for x in self.varvecs['k']:
            varlist.append(x)
        endocons_f = makeForward(self.varvecs['y'],'1','0')
        for x in endocons_f:
            varlist.append(x)
        exost_l = self.varvecs['z']
        for x in exost_l:
            varlist.append(x)
        endost_l = makeLags(self.varvecs['k'],'1')
        for x in endost_l:
            varlist.append(x)
        for x in self.varvecs['y']:
            varlist.append(x)
        varlist2 = list(self.subvars(tuple(varlist)))
        self.vlist = varlist2

        varlist=[]
        for x in self.varvecs['z']:
            varlist.append(x)
        for x in self.varvecs['k']:
            varlist.append(x)
        endocons_f = makeForward(self.varvecs['y'],'1','0')
        for x in endocons_f:
            varlist.append(x)
        exost_l = makeLags(self.varvecs['z'],'1')
        for x in exost_l:
            varlist.append(x)
        endost_l = makeLags(self.varvecs['k'],'1')
        for x in endost_l:
            varlist.append(x)
        for x in self.varvecs['y']:
            varlist.append(x)
        varlist2 = list(self.subvars(tuple(varlist)))
        self.vlist2 = varlist2

    def mkeqs2(self):
        nexost = len(self.varvecs['z'])
        nendost = len(self.varvecs['k'])
        nendocon = len(self.varvecs['y'])
        list_tmp1 = COP.deepcopy(self.interm4)
        eq_len = len(self.interm4)
        sub_list2 = self.vlist2
        sub_list = self.vlist
        reva2 = self.re_var2
        i1=0
        for x in list_tmp1:
            str_tmp1 = x[:]
            i2=0
            for y1,y2 in zip(sub_list,sub_list2):
                if i1 < eq_len-nexost:
                    while reva2.search(str_tmp1) != None:
                        ma1 = reva2.search(str_tmp1)
                        indx1 = sub_list.index(ma1.group(0))
                        str_tmp1 = str_tmp1.replace(ma1.group(0),'exp(x0('+str(indx1+1)+'))')[:]
                    i2 = i2 + 1
                else:
                    while reva2.search(str_tmp1) != None:
                        ma2 = reva2.search(str_tmp1)
                        indx2 = sub_list2.index(ma2.group(0))
                        str_tmp1 = str_tmp1.replace(ma2.group(0),'exp(x0('+str(indx2+1)+'))')[:]
                    i2 = i2 + 1
            list_tmp1[i1] = str_tmp1
            i1=i1+1

        self.interm5 = list_tmp1        

    def mkfunc(self):
        rege = RE.compile('\*{2,2}')
        inde = ' '*2
        eq_len=len(self.interm5)
        mfile=open(path_opts['MlabPath']['path']+'/Klein/'+'mfunc.m','w')
        mfile.write('function y = mfunc(x0)\n')
        mfile.write(inde+'%Parameter Values\n')
        for x in self.param.items():
            mfile.write(inde+x[0]+'='+str(x[1])+';\n')
        mfile.write('\n')
        mfile.write(inde+'%Steady State Values\n')
        for x in self.sstate_list:
            x[0]=rege.sub('^',x[0])
            x[1]=rege.sub('^',x[1])
            mfile.write(inde+x[0]+'='+str(x[1])+';\n')
        mfile.write('\n')
        mfile.write(inde+'y=zeros('+str(eq_len)+',1);\n')
        mfile.write('\n')
        for i1,i2 in zip(range(1,eq_len+1,1),self.interm5):
            mfile.write(inde+'y('+str(i1)+')='+i2+';'+'\n')
        mfile.close()

    def mkgrad(self):
        self.mkfunc()
        sess1 = self.sess1
        sstate = self.sstate
        ssdic = sstate[0]
        ssdic.update(sstate[1])
        ssdic2 = {}
        for x in ssdic.keys():
            if x[-4:] == '_bar':
                ssdic2[x] = ssdic[x]
        for x in ssdic2.keys():
            ssdic2[x] = np.log(ssdic2[x])
        vlist2 = []
        for x in self.vlist:
            ma = self.re_var2.search(x)
            tvar = ma.group('var')
            tvar = tvar.upper()
            tvar = tvar+'_bar'
            vlist2.append(tvar)
        invec = []
        for x in vlist2:
            invec.append(ssdic2[x])
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,'cd '+path_opts['MlabPath']['path'])
        mlabraw.eval(sess1,'cd Klein')
        mlabraw.put(sess1,'x0',invec)
        mlabraw.eval(sess1,"grdm=centgrad('mfunc',x0)")
        gradm = np.matrix(mlabraw.get(sess1,'grdm'))
        self.gradm = gradm

    def mkAB(self):
        nexost = len(self.varvecs['z'])
        nendost = len(self.varvecs['k'])
        nendocon = len(self.varvecs['y'])
        A = self.gradm[:,:nexost+nendost+nendocon]
        B = -self.gradm[:,nexost+nendost+nendocon:]
        self.AA = A
        self.BB = B

    def solab(self):
        A = self.AA
        B = self.BB
        tstates = len(self.varvecs['k'])+len(self.varvecs['z'])
        (F,P,retcon) = isolab(A,B,tstates,MAT.shape(A)[0])
        if MAT.sum(P.reshape(-1,1)) == 0.0:
            return
        else:
            self.P = np.matrix(P)
            self.F = np.matrix(F)

    def mkhess(self):
        self.mkfunc()
        sess1 = self.sess1
        sstate = self.sstate
        ssdic = sstate[0]
        ssdic.update(sstate[1])
        ssdic2 = {}
        for x in ssdic.keys():
            if x[-4:] == '_bar':
                ssdic2[x] = ssdic[x]
        for x in ssdic2.keys():
            ssdic2[x] = np.log(ssdic2[x])
        vlist2 = []
        for x in self.vlist:
            ma = self.re_var2.search(x)
            tvar = ma.group('var')
            tvar = tvar.upper()
            tvar = tvar+'_bar'
            vlist2.append(tvar)
        invec = []
        for x in vlist2:
            invec.append(ssdic2[x])
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,'cd '+path_opts['MlabPath']['path'])
        mlabraw.eval(sess1,'cd Klein')
        mlabraw.put(sess1,'x0',invec)
        mlabraw.eval(sess1,"hessm=centhess('mfunc',x0)")
        hessm = np.matrix(mlabraw.get(sess1,'hessm'))
        self.hessm = hessm

    def solab2(self):
        sess1 = self.sess1
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,'cd '+path_opts['MlabPath']['path'])
        mlabraw.eval(sess1,'cd Klein')
        mlabraw.put(sess1,'gradm',self.gradm)
        mlabraw.put(sess1,'hessm',self.hessm)
#----------------------------------------------------------------------------------------------------------------------
class MatWood:

    def __init__(self,intup):
        self.jAA = intup[0]
        self.jBB = intup[1]
        self.nexo = intup[2]
        self.ncon = intup[3]
        self.nendo = intup[4]
        self.sess1 = intup[5]
        self.NY = self.ncon+self.nendo
        self.NK = self.nendo
        self.NX = self.nexo
        self.mkwoodm()

    def solve(self):
        self.reds()
        self.solds()

    def mkwoodm(self):
        AA = self.jAA
        BB = self.jBB
        nendo = self.nendo
        nexo = self.nexo
        ncon = self.ncon
        nstates = nexo+nendo

        wAA = MAT.hstack((AA[:-nexo,nstates:],AA[:-nexo,nexo:nstates]))
        #wAA = MAT.hstack((AA[:,nstates:],AA[:,:nstates]))
        wBB = MAT.hstack((BB[:-nexo,nstates:],BB[:-nexo,nexo:nstates]))
        #wBB = MAT.hstack((BB[:,nstates:],BB[:,:nstates]))
        wCC = BB[:-nexo,:nexo]
        #wCC = MAT.zeros((ncon+nstates,1))


        self.wAA = wAA
        self.wBB = wBB
        self.wCC = wCC

    def reds(self):
        A = self.wAA
        B = self.wBB
        C = self.wCC
        NX = self.NX
        NK = self.NK
        NY = self.NY
        sess1 = self.sess1
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,'cd '+mlabpath)
        mlabraw.eval(sess1,'cd Woodford')
        mlabraw.put(sess1,'wAA',A)
        mlabraw.put(sess1,'wBB',B)
        mlabraw.put(sess1,'wCC',C)
        mlabraw.put(sess1,'NX',NX)
        mlabraw.put(sess1,'NY',NY)
        mlabraw.put(sess1,'NK',NK)
        mlabraw.eval(sess1,'[Br,Cr,Lr,NF] = redsf(wAA,wBB,wCC,NY,NX,NK)')
        self.Br = np.matrix(mlabraw.get(sess1,'Br'))
        self.Cr = np.matrix(mlabraw.get(sess1,'Cr'))
        self.Lr = np.matrix(mlabraw.get(sess1,'Lr'))
        self.NF = np.matrix(mlabraw.get(sess1,'NF'))

    def solds(self):
        Br = self.Br
        Cr = self.Cr
        Lr = self.Lr
        NF = self.NF
        NX = self.NX
        NK = self.NK
        NY = self.NY
        sess1 = self.sess1
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,'cd '+mlabpath)
        mlabraw.eval(sess1,'cd Woodford')
        mlabraw.put(sess1,'Br',Br)
        mlabraw.put(sess1,'Cr',Cr)
        mlabraw.put(sess1,'Lr',Lr)
        mlabraw.put(sess1,'NF',NF)
        mlabraw.put(sess1,'NX',NX)
        mlabraw.put(sess1,'NY',NY)
        mlabraw.put(sess1,'NK',NK)
        mlabraw.eval(sess1,'[D,F,G,H] = soldsf(Br,Cr,Lr,NY,NX,NK,NF)')
        self.D = np.matrix(mlabraw.get(sess1,'D'))
        self.F = np.matrix(mlabraw.get(sess1,'F'))
        self.G = np.matrix(mlabraw.get(sess1,'G'))
        self.H = np.matrix(mlabraw.get(sess1,'H'))
#----------------------------------------------------------------------------------------------------------------------
class ForKlein:

    def __init__(self,intup):
        self.uhlig = PyUhlig(intup)
        self.uhlig.mkmats()
        self.AA = self.uhlig.AA
        self.BB = self.uhlig.BB
        self.CC = self.uhlig.CC
        self.DD = self.uhlig.DD
        self.FF = self.uhlig.FF
        self.GG = self.uhlig.GG
        self.HH = self.uhlig.HH
        self.JJ = self.uhlig.JJ
        self.KK = self.uhlig.KK
        self.LL = self.uhlig.LL
        self.MM = self.uhlig.MM
        self.NN = self.uhlig.NN
        self.mkKleinMats()
        self.outdic = {}

    def mkKleinMats(self):
        # Make uhlig matrices locally available for computations
        AA = self.AA
        BB = self.BB
        CC = self.CC
        DD = self.DD
        FF = self.FF
        GG = self.GG
        HH = self.HH
        JJ = self.JJ
        KK = self.KK
        LL = self.LL
        MM = self.MM
        NN = self.NN

        # Determine size of states, endogenous
        exo_st = MAT.shape(NN)[1]
        endo_st = MAT.shape(BB)[1]
        endo_cn = MAT.shape(CC)[1]
        n_deteq = MAT.shape(AA)[0]
        n_expeq = MAT.shape(JJ)[0]
        tot_st = exo_st+endo_st
        self.tstates = tot_st
        tot_var = tot_st+endo_cn

        klein_A_rtwo = MAT.hstack((LL,GG))
        klein_A_rtwo = MAT.hstack((klein_A_rtwo,JJ))
        klein_A_rtwo_rows = MAT.shape(klein_A_rtwo)[0]
        klein_A_rtwo_cols = MAT.shape(klein_A_rtwo)[1]

        klein_A_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))
        klein_A_rone = MAT.hstack((MAT.identity(exo_st),klein_A_rone[:,exo_st:]))

        klein_A_rthree = MAT.hstack((MAT.zeros((n_deteq,exo_st)),AA))
        klein_A_rthree = MAT.hstack((klein_A_rthree,MAT.zeros((n_deteq,endo_cn))))

        klein_A = MAT.vstack((klein_A_rone,klein_A_rtwo))
        klein_A = MAT.vstack((klein_A,klein_A_rthree))

        klein_B_rone = MAT.zeros((exo_st,klein_A_rtwo_cols))

        klein_B_rone = MAT.hstack((NN,klein_B_rone[:,exo_st:]))
        klein_B_rtwo = MAT.hstack((-MM,-HH))
        klein_B_rtwo = MAT.hstack((klein_B_rtwo,-KK))
        klein_B_rthree = MAT.hstack((-DD,-BB))
        klein_B_rthree = MAT.hstack((klein_B_rthree,-CC))

        klein_B = MAT.vstack((klein_B_rone,klein_B_rtwo))
        klein_B = MAT.vstack((klein_B,klein_B_rthree))  

        self.A = klein_A
        self.B = klein_B

    def solve(self):
        A = self.A
        B = self.B
        tstates = MAT.shape(self.AA)[1] + MAT.shape(self.DD)[1]
        (F,P,retcon) = isolab(A,B,tstates,MAT.shape(A)[0])
        if MAT.sum(P.reshape(-1,1)) == 0.0:
            return
        else:
            self.P = np.matrix(P)
            self.F = np.matrix(F)
#----------------------------------------------------------------------------------------------------------------------
class PyKlein2D:

    def __init__(self,intup):
        self.gra = intup[0]
        self.hes = intup[1]
        self.nendo = intup[2]
        self.nexo = intup[3]
        self.ncon = intup[4]
        self.sigma = intup[5]
        self.A = intup[6]
        self.B = intup[7]
        self.vardic = intup[8]
        self.vdic = intup[9]
        self.modname = intup[10]
        self.audic = intup[11]
        self.nacon = len(self.audic['con']['var'])
        self.naendo = len(self.audic['endo']['var'])
        if self.vardic['other']['var']:
            self.numjl = intup[12]
            self.numhl = intup[13]
            self.nother = intup[14]
            self.oswitch = 1
            kintup = (self.gra,self.nendo,
                  self.nexo,self.ncon,
                  self.sigma,self.A,self.B,
                  self.vardic,self.vdic,
                  self.modname,self.audic,
                  self.numjl,self.nother)
        else:
            self.oswitch = 0
            kintup = (self.gra,self.nendo,
                  self.nexo,self.ncon,
                  self.sigma,self.A,self.B,
                  self.vardic,self.vdic,
                  self.modname,self.audic)
        self.tstates = self.nendo+self.nexo
        self.forkleind = ForKleinD(kintup)
        self.ssigma = self.mkssigma(self.tstates,self.nexo,self.sigma)

    def mkssigma(self,tstates,nexo,sigma):
        ssigma = MAT.zeros((tstates,tstates))
        for i in range(nexo):
            ssigma[i,i] = sigma[i,i]
        return ssigma

    def solab2(self):

        # Tracem function
        def tracem(xmat):
            n = MAT.shape(xmat)[1]
            m = MAT.shape(xmat)[0]/n
            ymat = MAT.zeros((m,1))
            for i1 in xrange(int(m)):
                ymat[i1,0]=MAT.trace(xmat[n*i1:i1*n+1,:n])
            return ymat

        ssigma = self.ssigma
        pp = self.P
        ff = self.F
        gra = self.gra
        hes = self.hes
        m = self.ncon+self.nendo+self.nexo
        ny = self.ncon
        nx = self.tstates

        f1 = gra[:,:nx]
        f2 = gra[:,nx:m]
        f4 = gra[:,m+nx:2*m]

        mm = MAT.vstack((pp,ff*pp,MAT.eye(nx),ff))
        aa1 = MAT.kron(MAT.eye(m),mm.T)*hes*mm
        bb1 = MAT.kron(f1,MAT.eye(nx))
        bb2 = MAT.kron(f2,MAT.eye(nx))
        bb4 = MAT.kron(f4,MAT.eye(nx))
        cc1 = MAT.kron(MAT.eye(ny),pp.T)
        cc2 = MAT.kron(ff,MAT.eye(nx))

        aa_1 = MAT.kron(MAT.eye(nx),bb4)+MAT.kron(pp.T,bb2*cc1)
        aa_2 = MAT.kron(MAT.eye(nx),bb1+bb2*cc2)
        aa = MAT.hstack((aa_1,aa_2))
        sol = np.linalg.solve(-aa,HLP.ctr(aa1.flatten(1)))
        ee = HLP.ctr(sol[:nx**2*ny].reshape(-1,nx*ny))
        gg = HLP.ctr(sol[nx**2*ny:].reshape(-1,nx**2))
        ma = 2.0*MAT.hstack((f1+f2*ff,f2+f4))
        eyeff = MAT.vstack((MAT.eye(nx),ff,MAT.zeros((m,nx))))
        ve_1 = f2*tracem(MAT.kron(MAT.eye(ny),ssigma)*ee)
        ve_2 = tracem(MAT.kron(MAT.eye(m),HLP.ctr(eyeff))*hes*eyeff*ssigma)
        ve = ve_1 + ve_2
        kxy = -(ma.I)*ve
        kx = kxy[:nx]
        ky = kxy[nx:]
        self.E = ee
        self.G = gg
        self.KX = kx
        self.KY = ky

    def solve(self):
        self.forkleind.solve()
        self.P = self.forkleind.P
        self.F = self.forkleind.F
        self.solab2()

    def sim(self,tlen,sntup=None):
        # Add 1000 observations, as they will be thrown away
        # Add one more observation to start first-order vector
        exoli = [x[1] for x in self.vardic['exo']['var']]
        indx=[]
        if sntup == None:
            indx = range(len(exoli))
        else:
            for name in sntup:
                if name not in exoli:
                    return 'Error: '+name+' is not a valid exo variable for this model!'
                else:
                    indx.append(exoli.index(name))
        indx.sort()
        ncon = self.ncon
        nexo = self.nexo
        nendo = self.nendo
        tstates = self.tstates
        tlena = 1000+tlen
        sigma = self.sigma
        ssigma = self.ssigma
        if self.oswitch:
            numjl = self.numjl
            numhl = self.numhl
        kx = self.KX
        ky = self.KY
        pp = self.P
        ff = self.F
        ee = self.E
        gg = self.G
        count = 0
        for varia in MAT.diag(sigma):
            if locals().has_key('ranvec'):
                if count in indx:
                    ranvec = MAT.vstack((ranvec,N.sqrt(varia)*MAT.matrix(N.random.standard_normal(tlena))))
                else:
                    ranvec = MAT.vstack((ranvec,MAT.zeros((1,tlena))))
            else:
                if count in indx:
                    ranvec = np.sqrt(varia)*MAT.matrix(N.random.standard_normal(tlena))
                else:
                    ranvec = MAT.zeros((1,tlena))
            count = count + 1

        ranvec = MAT.vstack((ranvec,MAT.zeros((nendo,tlena))))

        x_one_m1 = ranvec[:,0]
        x_one_0 = pp*x_one_m1
        x_two_m1 = kx + ranvec[:,0]
        y_one_0 = ff*x_one_m1
        y_one_m1 = MAT.zeros(y_one_0.shape)
        y_two_0 = ky+0.5*MAT.kron(MAT.eye(ncon),x_one_m1.T)*ee*x_one_m1
        if self.oswitch:
            nother = self.nother
            o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))
            o_two_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))+\
                0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1)).T)*\
                numhl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

        x_one_c = COP.deepcopy(x_one_m1)
        y_one_c = COP.deepcopy(y_one_0)
        x_two_c = COP.deepcopy(x_two_m1)
        y_two_c = COP.deepcopy(y_two_0)
        if self.oswitch:
            o_one_c = COP.deepcopy(o_one_0)
            o_two_c = COP.deepcopy(o_two_0)

        x_one = COP.deepcopy(x_one_c)
        y_one = COP.deepcopy(y_one_c)
        x_two = COP.deepcopy(x_two_c)
        y_two = COP.deepcopy(y_two_c)
        if self.oswitch:
            o_one = COP.deepcopy(o_one_c)
            o_two = COP.deepcopy(o_two_c)

        for i1 in range(1,ranvec.shape[1],1):

            x_one_n = pp*x_one_c+ranvec[:,i1]
            y_one_n = ff*x_one_c

            x_two_n = kx+pp*x_two_c+\
                0.5*MAT.kron(MAT.eye(tstates),x_one_c.T)*\
                gg*x_one_c+ranvec[:,i1]

            y_two_n = ky+ff*x_two_c+\
                0.5*MAT.kron(MAT.eye(ncon),x_one_c.T)*\
                ee*x_one_c

            if self.oswitch:
                o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))
                o_two_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))+\
                    0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c)).T)*\
                    numhl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

            x_one = MAT.hstack((x_one,x_one_c))
            y_one = MAT.hstack((y_one,y_one_n))
            x_two = MAT.hstack((x_two,x_two_c))
            y_two = MAT.hstack((y_two,y_two_n))
            if self.oswitch:
                o_one = MAT.hstack((o_one,o_one_n))
                o_two = MAT.hstack((o_two,o_two_n))

            x_one_c = x_one_n
            y_one_c = y_one_n
            x_two_c = x_one_n
            y_two_c = y_two_n
            if self.oswitch:
                o_one_c = o_one_n
                o_two_c = o_two_n

        # Throw away first 1000 obs
        x_one = x_one[:,1000:]
        y_one = y_one[:,1000:]
        y_two = y_two[:,1000:]
        x_two = x_two[:,1000:]
        if self.oswitch:
            o_one = o_one[:,1000:]
            o_two = o_two[:,1000:]

        self.sim_y_two = y_two
        self.sim_y_one = y_one
        self.sim_x_two = x_two
        self.sim_x_one = x_one
        self.insim = [self.sim_x_two,self.sim_y_two]
        if self.oswitch:
            self.sim_o_one = o_one
            self.sim_o_two = o_two
            self.insim = self.insim + [self.sim_o_two,]

    def mkcvar(self,vname=False,insim='insim'):
        # Now also produce unconditional variances, and normal. variances with v
        # Check if simul data is available
        if insim not in dir(self):
            return "Error: Can't produce uncond. variance table without simulation data!"

        exoli = [x[1] for x in self.vardic['exo']['var']]
        endoli = [x[1] for x in self.vardic['endo']['var']]
        stateli = exoli + endoli
        conli = [x[1] for x in self.vardic['con']['var']]
        alli = stateli+conli
        if self.oswitch:
            otherli = [x[1] for x in self.vardic['other']['var']]
            alli = alli + otherli

        # Check if vname is valid
        if vname not in alli:
            return 'Error: '+vname+' is not a valid variable for this model!'
        insim = eval('self.'+insim)
        osim_x = COP.deepcopy(insim[0])
        osim_y = COP.deepcopy(insim[1])
        nendo = self.nendo
        nexo = self.nexo
        ncon = self.ncon

        # Now hp filter the simulations before graphing according to filtup
        for i1 in xrange(osim_y.shape[0]):
            yy = osim_y[i1,:].__array__().T
            woy = np.zeros((osim_y.shape[1],3))
            lam = 1600
            yyf = MAT.matrix(hpfilter(data=yy,lam=1600))
            osim_y[i1,:] = yyf[0]
        # Now filter the state variables!
        for i1 in xrange(osim_x.shape[0]):
            xx = osim_x[i1,:].__array__().T
            wox = np.zeros((osim_x.shape[1],3))
            lam = 1600
            xxf = MAT.matrix(hpfilter(data=xx,lam=1600))
            osim_x[i1,:] = xxf[0]

        if self.oswitch:
            osim_o = COP.deepcopy(insim[2])
            # Now hp filter the other variables!
            for i1 in xrange(osim_o.shape[0]):
                oo = osim_o[i1,:].__array__().T
                woo = np.zeros((osim_o.shape[1],3))
                lam = 1600
                oof = MAT.matrix(hpfilter(data=oo,lam=1600))
                osim_o[i1,:] = oof[0]

        # Stack all matrices into one
        allmat = MAT.vstack((osim_x,osim_y))
        if self.oswitch:
            allmat = MAT.vstack((allmat,osim_o))
        varmat = MAT.zeros((allmat.shape[0],1))
        for i1 in range(allmat.shape[0]):
            varmat[i1,0] = np.var(allmat[i1,:].A.flatten(1))

        if vname:
            relmat = MAT.zeros((allmat.shape[0],1))
            indx = alli.index(vname)
            varvari = np.var(allmat[indx,:].A.flatten(1))
            for i1 in range(allmat.shape[0]):
                relmat[i1,0] = varmat[i1,0]/varvari
            self.cvarm = MAT.hstack((varmat,relmat))
        else:
            self.cvarm = varmat



    def mkact(self,vname,intup=None,insim='insim'):
        # Now also produce autocorrelation table
        # Check if simul data is available
        if insim not in dir(self):
            return "Error: Can't produce autocorrelation table without simulation data!"

        exoli = [x[1] for x in self.vardic['exo']['var']]
        endoli = [x[1] for x in self.vardic['endo']['var']]
        stateli = exoli + endoli
        conli = [x[1] for x in self.vardic['con']['var']]
        alli = stateli+conli
        if self.oswitch:
            otherli = [x[1] for x in self.vardic['other']['var']]
            alli = alli + otherli
        # Check if vname is valid
        if vname not in alli:
            return 'Error: '+vname+' is not a valid variable for this model!'

        if intup == None:
            lags = 3
            leads = 3
        else:
            lags = intup[0]
            leads = intup[1]

        insim = eval('self.'+insim)
        osim_x = COP.deepcopy(insim[0])
        osim_y = COP.deepcopy(insim[1])

        # Now hp filter the simulations before graphing according to filtup
        for i1 in xrange(osim_y.shape[0]):
            yy = osim_y[i1,:].__array__().T
            woy = np.zeros((osim_y.shape[1],3))
            lam = 1600
            yyf = MAT.matrix(hpfilter(data=yy,lam=1600))
            osim_y[i1,:] = yyf[0]
        # Now filter the state variables!
        for i1 in xrange(osim_x.shape[0]):
            xx = osim_x[i1,:].__array__().T
            wox = np.zeros((osim_x.shape[1],3))
            lam = 1600
            xxf = MAT.matrix(hpfilter(data=xx,lam=1600))
            osim_x[i1,:] = xxf[0]

        sim_x = osim_x[:,lags:-leads]
        sim_xf = osim_x[:,leads+1:]
        sim_yf = osim_y[:,leads+1:]
        sim_y = osim_y[:,lags:-leads]
        sim_xl = osim_x[:,:-lags-1]
        sim_yl = osim_y[:,:-lags-1]

        if self.oswitch:
            osim_o = COP.deepcopy(insim[2])
            # Now hp filter the other variables!
            for i1 in xrange(osim_o.shape[0]):
                oo = osim_o[i1,:].__array__().T
                woo = np.zeros((osim_o.shape[1],3))
                lam = 1600
                oof = MAT.matrix(hpfilter(data=oo,lam=1600))
                osim_o[i1,:] = oof[0]
            sim_o = osim_o[:,lags:-leads]
            sim_of = osim_o[:,leads+1:]
            sim_ol = osim_o[:,:-lags-1]
            alli = alli + otherli

        posa = alli.index(vname)
        actm = MAT.zeros((len(alli),lags+1+leads))
        sima = MAT.vstack((sim_x,sim_y))
        sima_f = MAT.vstack((sim_xf,sim_yf))
        sima_l = MAT.vstack((sim_xl,sim_yl))
        if self.oswitch:
            sima = MAT.vstack((sima,sim_o))
            sima_f = MAT.vstack((sima_f,sim_of))
            sima_l = MAT.vstack((sima_l,sim_ol))
        for i1 in range(sima.shape[0]):
            # Do lags
            for i2 in range(lags):
                actm[i1,i2] = np.round(np.corrcoef(sima[posa,:].A.flatten(1),sima_l[i1,i2:sima_l.shape[1]-lags+i2+1].A.flatten(1))[1,0],3)
            # Do current
            actm[i1,lags] = np.round(np.corrcoef(sima[posa,:].A.flatten(1),sima[i1,:].A.flatten(1))[1,0],3)
            # Do leads
            for i2 in range(leads):
                actm[i1,lags+1+i2] = np.round(np.corrcoef(sima[posa,:].A.flatten(1),sima_f[i1,i2:sima_f.shape[1]-leads+i2+1].A.flatten(1))[1,0],3)

        self.actm = actm
        self.actname = vname
        if intup:
            self.actin = intup

    def show_act(self):
        if 'actname' not in dir(self):
            return 'Error: You have not produced any autocorrelation table yet!'
        nexo = self.nexo
        nendo = self.nendo
        ncon = self.ncon
        if 'naendo' in dir(self):
            naendo = self.naendo
        else:
            naendo = 0
        if 'nacon' in dir(self):
            nacon = self.nacon
        else:
            nacon = 0
        nother = len(self.vardic['other']['var'])
        respva = COP.deepcopy(self.actname)
        if 'actin' in dir(self):
            lags = COP.deepcopy(int(self.actin[0]))
            leads = COP.deepcopy(int(self.actin[1]))
        else:
            lags = 3
            leads = 3
        actm = COP.deepcopy(self.actm)
        vdic = self.vdic
        allvari = [x[1] for x in vdic['exo']]
        allvari = allvari + [x[1] for x in vdic['endo']]
        allvari = allvari + [x[1] for x in vdic['con']]
        if self.oswitch:
            allvari = allvari + [x[1] for x in vdic['other']]
        allvari = allvari[:nexo+nendo-naendo]+allvari[nexo+nendo:nexo+nendo+ncon-nacon]
        actm = MAT.vstack((actm[:nexo+nendo-naendo,:],actm[nexo+nendo:nexo+nendo+ncon-nacon,:]))

        # determine longest varname
        vnmax = 0
        for vari in allvari:
            if len(vari) > vnmax:
                vnmax = len(vari)
        print ''
        print 'Autocorrelation table, current '+respva
        print '='*65
        for i1,vari in enumerate(allvari):
            print vari+(vnmax-len(vari))*' '+'  |'+str(actm[i1,:])[2:-2]

    def show_sim(self,intup,filtup=None,insim='insim'):
        # Check if simulations have been carried out
        if insim not in dir(self):
            return 'Error: You have not produced any simulations yet! Nothing to graph!'
        mname = self.modname
        vardic = self.vardic
        insim = eval('self.'+insim)
        sim_x = COP.deepcopy(insim[0])
        sim_y = COP.deepcopy(insim[1])
        tlen = sim_x.shape[1]
        if self.oswitch:
            sim_o = COP.deepcopy(insim[2])
            nother = self.nother
            otherli = [x[1] for x in vardic['other']['var']]
        nexo = self.nexo
        nendo = self.nendo
        ncon = self.ncon
        conli = [x[1] for x in vardic['con']['var']]
        endoli = [x[1] for x in vardic['endo']['var']]
        exoli = [x[1] for x in vardic['exo']['var']]

        # If filtup is None, create filtup from vardic
        if filtup == None:
            filli = [0]*len(intup)
            for i1,x in enumerate(intup):
                for y in vardic.keys():
                    if x in [z[1] for z in vardic[y]['var']]:
                        indx = [z[1] for z in vardic[y]['var']].index(x)
                        if 'hp' in [z for z in vardic[y]['mod']][indx]:
                            filli[i1] = 1
                        elif 'bk' in [z for z in vardic[y]['mod']][indx]:
                            filli[i1] = 2
                        elif 'cf' in [z for z in vardic[y]['mod']][indx]:
                            filli[i1] = 3
            filtup = tuple(filli)


        stateli = exoli+endoli
        alli = stateli+conli
        if self.oswitch:
            alli = alli+otherli
        # Check if all name in intup are available to graph
        for name in intup:
            if name not in alli:
                return 'Error: '+name+' is not a valid variable name for this model!'
        # Create x, y and z indeces
        indx = []
        indy = []
        indo = []
        for name in intup:
            if name in stateli:
                indx.append(stateli.index(name))
            elif name in conli:
                indy.append(conli.index(name))
            elif self.oswitch and name in otherli:
                indo.append(otherli.index(name))
        leg = []
        indx.sort()
        indy.sort()
        indo.sort()     

        # Now hp filter the simulations before graphing according to filtup
        for i1 in xrange(sim_y.shape[0]):
            if indy and (i1 in indy) and filtup[list(intup).index(conli[i1])]:
                if filtup[list(intup).index(conli[i1])] == 1:
                    conli[i1] = conli[i1]+'(hpf)'
                    yy = sim_y[i1,:].__array__().T
                    woy = np.zeros((tlen,3))
                    lam = 1600
                    yyf = MAT.matrix(hpfilter(data=yy,lam=1600))
                    sim_y[i1,:] = yyf[0]*100.0
                elif filtup[list(intup).index(conli[i1])] == 2:
                    conli[i1] = conli[i1]+'(bkf)'
                    yy = sim_y[i1,:].__array__().T
                    woy = np.zeros((tlen,1))
                    up = 6
                    dn = 32
                    kkl = 12
                    yyf = MAT.matrix(bkfilter(data=yy,up=up,dn=dn,kkl=kkl))
                    sim_y[i1,:] = yyf[0]*100.0
                elif filtup[list(intup).index(conli[i1])] == 3:
                    conli[i1] = conli[i1]+'(cff)'
                    yy = sim_y[i1,:].__array__().T
                    woy = np.zeros((tlen,1))
                    low = 6
                    high = 32
                    drift = True
                    yyf = MAT.matrix(cffilter(data=yy,low=low,high=high,drift=drift))
                    sim_y[i1,:] = yyf[0]*100.0
        # Now filter the state variables!
        for i1 in xrange(sim_x.shape[0]):
            if indx and (i1 in indx) and filtup[list(intup).index(stateli[i1])]:
                if filtup[list(intup).index(stateli[i1])] == 1:
                    stateli[i1] = stateli[i1]+'(hpf)'
                    xx = sim_x[i1,:].__array__().T
                    wox = np.zeros((tlen,3))
                    lam = 1600
                    xxf = MAT.matrix(hpfilter(data=xx,lam=1600))
                    sim_x[i1,:] = xxf[0]*100.0
                elif filtup[list(intup).index(stateli[i1])] == 2:
                    stateli[i1] = stateli[i1]+'(bkf)'
                    xx = sim_x[i1,:].__array__().T
                    wox = np.zeros((tlen,1))
                    up = 6
                    dn = 32
                    kkl = 12
                    xxf = MAT.matrix(bkfilter(data=xx,up=up,dn=dn,kkl=kkl))
                    sim_x[i1,:] = xxf[0]*100.0
                elif filtup[list(intup).index(stateli[i1])] == 3:
                    stateli[i1] = stateli[i1]+'(cff)'
                    xx = sim_x[i1,:].__array__().T
                    wox = np.zeros((tlen,1))
                    low = 6
                    high = 32
                    drift = True
                    xxf = MAT.matrix(cffilter(data=xx,low=low,high=high,drift=drift))
                    sim_x[i1,:] = xxf[0]*100.0 
        # Now hp filter the other variables!
        if self.oswitch:
            for i1 in xrange(sim_o.shape[0]):
                if indo and (i1 in indo) and filtup[list(intup).index(otherli[i1])]:
                    if filtup[list(intup).index(otherli[i1])] == 1:                   
                        otherli[i1] = otherli[i1]+'(hpf)'
                        oo = sim_o[i1,:].__array__().T
                        woo = np.zeros((tlen,3))
                        lam = 1600
                        oof = MAT.matrix(hpfilter(data=oo,lam=1600))
                    elif filtup[list(intup).index(otherli[i1])] == 2:
                        otherli[i1] = otherli[i1]+'(bkf)'
                        oo = sim_o[i1,:].__array__().T
                        woo = np.zeros((tlen,1))
                        up = 6
                        dn = 32
                        kkl = 12
                        oof = MAT.matrix(bkfilter(data=oo,up=up,dn=dn,kkl=kkl))
                    elif filtup[list(intup).index(otherli[i1])] == 3:
                        otherli[i1] = otherli[i1]+'(cff)'
                        oo = sim_o[i1,:].__array__().T
                        woo = np.zeros((tlen,1))
                        low = 6
                        high = 32
                        drift = True
                        oof = MAT.matrix(cffilter(data=oo,low=low,high=high,drift=drift)) 

        if indx and indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_x.T[:,indx],sim_y.T[:,indy],sim_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and indo:
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_y.T[:,indy],sim_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_x.T[:,indx],sim_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_x.T[:,indx],sim_y.T[:,indy])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(sim_x.T[:,indx].A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and not indo:
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(sim_y.T[:,indy].A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and not indy and indo:
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(sim_o.T[:,indo].A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        return figo

    def save_sim(self,intup,fpath=None,filtup=None,insim='insim'):
        # Check if simulations have been carried out
        if insim not in dir(self):
            return 'Error: You have not produced any simulations yet! Nothing to graph!'
        if fpath == None:
            return 'Error: You have not supplied the function with an fpath= argument where to save the plot!'
        mname = self.modname
        vardic = self.vardic
        insim = eval('self.'+insim)
        sim_x = COP.deepcopy(insim[0])
        sim_y = COP.deepcopy(insim[1])
        tlen = sim_x.shape[1]
        if self.oswitch:
            sim_o = COP.deepcopy(insim[2])
            nother = self.nother
            otherli = [x[1] for x in vardic['other']['var']]
        nexo = self.nexo
        nendo = self.nendo
        ncon = self.ncon
        conli = [x[1] for x in vardic['con']['var']]
        endoli = [x[1] for x in vardic['endo']['var']]
        exoli = [x[1] for x in vardic['exo']['var']]

        # If filtup is None, create filtup from vardic
        if filtup == None:
            filli = [0]*len(intup)
            for i1,x in enumerate(intup):
                for y in vardic.keys():
                    if x in [z[1] for z in vardic[y]['var']]:
                        indx = [z[1] for z in vardic[y]['var']].index(x)
                        if 'hp' in [z for z in vardic[y]['mod']][indx]:
                            filli[i1] = 1
                        elif 'bk' in [z for z in vardic[y]['mod']][indx]:
                            filli[i1] = 2
            filtup = tuple(filli)


        stateli = exoli+endoli
        alli = stateli+conli
        if self.oswitch:
            alli = alli+otherli
        # Check if all name in intup are available to graph
        for name in intup:
            if name not in alli:
                return 'Error: '+name+' is not a valid variable name for this model!'
        # Create x, y and z indeces
        indx = []
        indy = []
        indo = []
        for name in intup:
            if name in stateli:
                indx.append(stateli.index(name))
            elif name in conli:
                indy.append(conli.index(name))
            elif self.oswitch and name in otherli:
                indo.append(otherli.index(name))
        leg = []
        indx.sort()
        indy.sort()
        indo.sort()     

        # Now hp filter the simulations before graphing according to filtup
        for i1 in xrange(sim_y.shape[0]):
            if indy and (i1 in indy) and filtup[list(intup).index(conli[i1])]:
                if filtup[list(intup).index(conli[i1])] == 1:
                    conli[i1] = conli[i1]+'(hpf)'
                    yy = sim_y[i1,:].__array__().T
                    woy = np.zeros((tlen,3))
                    lam = 1600
                    yyf = MAT.matrix(hpfilter(data=yy,lam=1600))
                    sim_y[i1,:] = yyf[0]*100
                elif filtup[list(intup).index(conli[i1])] == 2:
                    conli[i1] = conli[i1]+'(bkf)'
                    yy = sim_y[i1,:].__array__().T
                    woy = np.zeros((tlen,1))
                    up = 6
                    dn = 32
                    kkl = 12
                    yyf = MAT.matrix(bkfilter(yy,up,dn,kkl,tlen))
                    sim_y[i1,:] = yyf[0]*100
        # Now filter the state variables!
        for i1 in xrange(sim_x.shape[0]):
            if indx and (i1 in indx) and filtup[list(intup).index(stateli[i1])]:
                if filtup[list(intup).index(stateli[i1])] == 1:
                    stateli[i1] = stateli[i1]+'(hpf)'
                    xx = sim_x[i1,:].__array__().T
                    wox = np.zeros((tlen,3))
                    lam = 1600
                    xxf = MAT.matrix(hpfilter(data=xx,lam=1600))
                    sim_x[i1,:] = xxf[0]*100
                elif filtup[list(intup).index(stateli[i1])] == 2:
                    stateli[i1] = stateli[i1]+'(bkf)'
                    xx = sim_x[i1,:].__array__().T
                    wox = np.zeros((tlen,1))
                    up = 6
                    dn = 32
                    kkl = 12
                    xxf = MAT.matrix(bkfilter(xx,up,dn,kkl,tlen))
                    sim_x[i1,:] = xxf[0]*100                  
            # Now hp filter the other variables!
        if self.oswitch:
            for i1 in xrange(sim_o.shape[0]):
                if indo and (i1 in indo) and filtup[list(intup).index(otherli[i1])]:
                    if filtup[list(intup).index(otherli[i1])] == 1:                   
                        otherli[i1] = otherli[i1]+'(hpf)'
                        oo = sim_o[i1,:].__array__().T
                        woo = np.zeros((tlen,3))
                        lam = 1600
                        oof = MAT.matrix(hpfilter(data=oo,lam=1600))
                    elif filtup[list(intup).index(otherli[i1])] == 2:
                        otherli[i1] = otherli[i1]+'(bkf)'
                        oo = sim_o[i1,:].__array__().T
                        woo = np.zeros((tlen,1))
                        up = 6
                        dn = 32
                        kkl = 12
                        oof = MAT.matrix(bkfilter(oo,up,dn,kkl,tlen)) 

        if indx and indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_x.T[:,indx],sim_y.T[:,indy],sim_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and indo:
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_y.T[:,indy],sim_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_x.T[:,indx],sim_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(MAT.hstack((sim_x.T[:,indx],sim_y.T[:,indy])).A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(sim_x.T[:,indx].A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and not indo:
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(sim_y.T[:,indy].A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and not indy and indo:
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            P.plot(sim_o.T[:,indo].A)
            P.title(str(tlen)+' simulated periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage deviation from SS')
            P.grid()
            P.legend(leg)
        PLT.savefig(fpath,dpi=100)
        return "Your plot has been saved in: "+fpath


    def irf(self,tlen,sntup):
        tlen = tlen + 1
        ncon = self.ncon
        nexo = self.nexo
        nendo = self.nendo
        tstates = self.tstates
        sigma = self.sigma
        ssigma = self.ssigma
        if self.oswitch:
            numjl = self.numjl
            numhl = self.numhl
        kx = self.KX
        ky = self.KY
        pp = self.P
        ff = self.F
        ee = self.E
        gg = self.G

        sposli=[]
        exoli = [x[1] for x in self.vardic['exo']['var']]
        # Check if names are valid
        for name in sntup:
            if name not in exoli:
                return 'Error: '+name+' is not a valid exoshock name for this model!'
        for name in sntup:
            sposli.append(exoli.index(name))
        # Expose spos choice to self for show_irf
        self.spos = (sntup,sposli)

        shock = MAT.zeros((nexo,1))
        sendo = MAT.zeros((nendo,1))

        for spos in sposli:
            shock[spos,0] = 1.0
        shock = MAT.vstack((shock,sendo))

        x_one_m1 = shock
        x_one_0 = pp*x_one_m1
        x_two_m1 = shock
        y_one_0 = ff*x_one_m1
        y_one_m1 = MAT.zeros(y_one_0.shape)
        y_two_0 = 0.5*MAT.kron(MAT.eye(ncon),x_one_m1.T)*ee*x_one_m1
        if self.oswitch:
            nother = self.nother
            o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))
            o_two_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))+\
                0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1)).T)*\
                numhl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

        x_one_c = COP.deepcopy(x_one_m1)
        y_one_c = COP.deepcopy(y_one_0)
        x_two_c = COP.deepcopy(x_two_m1)
        y_two_c = COP.deepcopy(y_two_0)
        if self.oswitch:
            o_one_c = COP.deepcopy(o_one_0)
            o_two_c = COP.deepcopy(o_two_0)

        x_one = COP.deepcopy(x_one_c)
        y_one = COP.deepcopy(y_one_c)
        x_two = COP.deepcopy(x_two_c)
        y_two = COP.deepcopy(y_two_c)
        if self.oswitch:
            o_one = COP.deepcopy(o_one_c)
            o_two = COP.deepcopy(o_two_c)

        for i1 in range(tlen):
            x_one_n = pp*x_one_c
            y_one_n = ff*x_one_c
            x_two_n = pp*x_two_c+0.5*MAT.kron(MAT.eye(tstates),x_one_c.T)*gg*x_one_c
            y_two_n = ff*x_two_c+0.5*MAT.kron(MAT.eye(ncon),x_one_c.T)*ee*x_one_c
            if self.oswitch:
                o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))
                o_two_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))+\
                    0.5*MAT.kron(MAT.eye(nother),MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c)).T)*\
                    numhl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

            x_one = MAT.hstack((x_one,x_one_c))
            y_one = MAT.hstack((y_one,y_one_n))
            x_two = MAT.hstack((x_two,x_two_c))
            y_two = MAT.hstack((y_two,y_two_n))
            if self.oswitch:
                o_one = MAT.hstack((o_one,o_one_n))
                o_two = MAT.hstack((o_two,o_two_n))

            x_one_c = COP.deepcopy(x_one_n)
            y_one_c = COP.deepcopy(y_one_n)
            x_two_c = COP.deepcopy(x_one_n)
            y_two_c = COP.deepcopy(y_two_n)

        self.irf_x_two = x_two[:,1:-1]
        self.irf_y_two = y_two[:,1:-1]
        self.inirf = [self.irf_x_two,self.irf_y_two]
        if self.oswitch:
            self.irf_o_one = o_one[:,2:]
            self.irf_o_two = o_two[:,2:]
            self.inirf = self.inirf + [self.irf_o_two,]

    def show_irf(self,intup,inirf='inirf'):
        # Check if simulations have been carried out
        if inirf not in dir(self):
            return 'Error: You have not produced any IRFs yet! Nothing to graph!'
        inirf = eval('self.'+inirf)
        irf_x = COP.deepcopy(inirf[0])
        tlen = irf_x.shape[1]
        irf_x = MAT.hstack((MAT.zeros((irf_x.shape[0],20)),irf_x))
        irf_y = COP.deepcopy(inirf[1])
        irf_y = MAT.hstack((MAT.zeros((irf_y.shape[0],20)),irf_y))
        irf_x = irf_x*100.0
        irf_y = irf_y*100.0
        mname = self.modname
        vardic = self.vardic
        time_axis = np.arange(-20,tlen,1)
        if self.oswitch:
            irf_o = COP.deepcopy(inirf[2])
            irf_o = MAT.hstack((MAT.zeros((irf_o.shape[0],20)),irf_o))
            nother = self.nother
            varother = self.vardic['other']['var']
            otherli = [x[1] for x in varother]
        nexo = self.nexo
        nendo = self.nendo
        ncon = self.ncon
        conli = [x[1] for x in vardic['con']['var']]
        endoli = [x[1] for x in vardic['endo']['var']]
        exoli = [x[1] for x in vardic['exo']['var']]

        stateli = exoli+endoli
        alli = stateli+conli
        if self.oswitch:
            alli = alli+otherli
        # Check if all name in intup are available to graph
        for name in intup:
            if name not in alli:
                return 'Error: '+name+' is not a valid variable name for this model!'
        # Create x, y and z indeces
        indx = []
        indy = []
        indo = []
        for name in intup:
            if name in stateli:
                indx.append(stateli.index(name))
            elif name in conli:
                indy.append(conli.index(name))
            elif self.oswitch and name in otherli:
                indo.append(otherli.index(name))
        leg = []
        indx.sort()
        indy.sort()
        indo.sort()

        if indx and indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_x.T[:,indx],irf_y.T[:,indy],irf_o.T[:,indo])).A)
            P.title(str(tlen)+' IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and indo:
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_y.T[:,indy],irf_o.T[:,indo])).A)
            P.title(str(tlen)+' IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_x.T[:,indx],irf_o.T[:,indo])).A)
            P.title(str(tlen)+' IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_x.T[:,indx],irf_y.T[:,indy])).A)
            P.title(str(tlen)+' IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,irf_x.T[:,indx].A)
            P.title(str(tlen)+' IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and not indo:
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,irf_y.T[:,indy].A)
            P.title(str(tlen)+' IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and not indy and indo:
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,irf_o.T[:,indo].A)
            P.title(str(tlen)+' IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        return figo

    def save_irf(self,intup,fpath=None,inirf='inirf'):
        # Check if simulations have been carried out
        if inirf not in dir(self):
            return 'Error: You have not produced any IRFs yet! Nothing to graph!'
        if fpath == None:
            return 'Error: You have not supplied the function with an fpath= argument where to save the plot!'
        inirf = eval('self.'+inirf)
        irf_x = COP.deepcopy(inirf[0])
        tlen = irf_x.shape[1]
        irf_x = MAT.hstack((MAT.zeros((irf_x.shape[0],20)),irf_x))
        irf_y = COP.deepcopy(inirf[1])
        irf_y = MAT.hstack((MAT.zeros((irf_y.shape[0],20)),irf_y))
        irf_x = irf_x*100
        irf_y = irf_y*100
        mname = self.modname
        vardic = self.vardic
        time_axis = np.arange(-20,tlen,1)
        if self.oswitch:
            irf_o = COP.deepcopy(inirf[2])
            nother = self.nother
            varother = self.vardic['other']['var']
            otherli = [x[1] for x in varother]
        nexo = self.nexo
        nendo = self.nendo
        ncon = self.ncon
        conli = [x[1] for x in vardic['con']['var']]
        endoli = [x[1] for x in vardic['endo']['var']]
        exoli = [x[1] for x in vardic['exo']['var']]

        stateli = exoli+endoli
        alli = stateli+conli
        if self.oswitch:
            alli = alli+otherli
        # Check if all name in intup are available to graph
        for name in intup:
            if name not in alli:
                return 'Error: '+name+' is not a valid variable name for this model!'
        # Create x, y and z indeces
        indx = []
        indy = []
        indo = []
        for name in intup:
            if name in stateli:
                indx.append(stateli.index(name))
            elif name in conli:
                indy.append(conli.index(name))
            elif self.oswitch and name in otherli:
                indo.append(otherli.index(name))
        leg = []
        indx.sort()
        indy.sort()
        indo.sort()

        if indx and indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_x.T[:,indx],irf_y.T[:,indy],irf_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and indo:
            for y in indy:
                leg.append(conli[y])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_y.T[:,indy],irf_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and indo:
            for x in indx:
                leg.append(stateli[x])
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_x.T[:,indx],irf_o.T[:,indo])).A)
            P.title(str(tlen)+' simulated IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,MAT.hstack((irf_x.T[:,indx],irf_y.T[:,indy])).A)
            P.title(str(tlen)+' simulated IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif indx and not indy and not indo:
            for x in indx:
                leg.append(stateli[x])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,irf_x.T[:,indx].A)
            P.title(str(tlen)+' simulated IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS')
            P.grid()
            P.legend(leg)
        elif not indx and indy and not indo:
            for y in indy:
                leg.append(conli[y])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,irf_y.T[:,indy].A)
            P.title(str(tlen)+' simulated IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS, hp-filtered')
            P.grid()
            P.legend(leg)
        elif not indx and not indy and indo:
            for o in indo:
                leg.append(otherli[o])
            leg = tuple(leg)
            figo = P.figure()
            ax = figo.add_subplot(111)
            ax.set_xlim([-20,tlen])
            ax.plot(time_axis,irf_o.T[:,indo].A)
            P.title(str(tlen)+' simulated IRF periods, '+mname)
            P.xlabel('Time')
            P.ylabel('Percentage Deviation from SS, hp-filtered')
            P.grid()
            P.legend(leg)
        PLT.savefig(fpath,dpi=100)
        return "Your plot has been save in: "+fpath

#----------------------------------------------------------------------------------------------------------------------
class MatKlein2D(PyKlein2D):

    def __init__(self,intup):
        self.sess1 = intup[-1]
        intup2 = intup[:-1]
        # Get all PyKlein2D attributes
        PyKlein2D.__init__(self, intup2)

    def solve(self):
        tstates = self.tstates
        sess1 = self.sess1
        ssigma = self.ssigma
        mlabraw.eval(sess1,'clear all;')
        mlabraw.eval(sess1,'cd '+mlabpath)
        mlabraw.eval(sess1,'cd Klein')
        mlabraw.put(sess1,'jac',self.gra)
        mlabraw.put(sess1,'hess',self.hes)
        mlabraw.put(sess1,'nstates',self.tstates)
        mlabraw.put(sess1,'ssigma',self.ssigma)
        mlabraw.eval(sess1,'[ff,pp,ee,gg,kx,ky] = solab2(jac,hess,ssigma,nstates)')
        self.F = np.matrix(mlabraw.get(sess1,'ff'))
        self.P = np.matrix(mlabraw.get(sess1,'pp'))
        self.E = np.matrix(mlabraw.get(sess1,'ee'))
        self.G = np.matrix(mlabraw.get(sess1,'gg'))
        self.KX = np.matrix(mlabraw.get(sess1,'kx'))
        self.KY = np.matrix(mlabraw.get(sess1,'ky'))
#----------------------------------------------------------------------------------------------------------------------
class ForKleinD(PyKlein2D):

    def __init__(self,intup):
        self.gra = intup[0]
        self.nendo = intup[1]
        self.nexo = intup[2]
        self.ncon = intup[3]
        self.sigma = intup[4]
        self.A = intup[5]
        self.B = intup[6]
        self.vardic = intup[7]
        self.vdic = intup[8]
        self.modname = intup[9]
        self.audic = intup[10]
        if self.vardic['other']['var']:
            self.numjl = intup[11]
            self.nother = intup[12]
            self.oswitch = 1
        else:
            self.oswitch = 0
        self.tstates = self.nendo+self.nexo
        self.ssigma = self.mkssigma(self.tstates,self.nexo,self.sigma)

    def solve(self):
        A = self.A
        B = self.B
        tstates = self.tstates
        (F,P,retcon) = isolab(A,B,tstates,MAT.shape(A)[0])
        if MAT.sum(P.reshape(-1,1)) == 0.0:
            return
        else:
            self.P = np.matrix(P)
            self.F = np.matrix(F)

    def sim(self,tlen,sntup=None,shockvec=None):
        # Add 1000 observations, as they will be thrown away
        # Add one more observation to start first-order vector
        exoli = [x[1] for x in self.vardic['exo']['var']]
        indx=[]
        if sntup == None:
            indx = range(len(exoli))
        else:
            for name in sntup:
                if name not in exoli:
                    return 'Error: '+name+' is not a valid exo variable for this model!'
                else:
                    indx.append(exoli.index(name))
        indx.sort()
        ncon = self.ncon
        nexo = self.nexo
        nendo = self.nendo
        tstates = self.tstates
        tlena = 1000+tlen
        sigma = self.sigma
        ssigma = self.ssigma
        if self.oswitch:
            numjl = self.numjl

        pp = self.P
        ff = self.F

        if shockvec == None:
            count = 0
            for varia in MAT.diag(sigma):
                if locals().has_key('ranvec'):
                    if count in indx:
                        ranvec = MAT.vstack((ranvec,np.sqrt(varia)*MAT.matrix(np.random.standard_normal(tlena))))
                    else:
                        ranvec = MAT.vstack((ranvec,MAT.zeros((1,tlena))))      
                else:
                    if count in indx:
                        ranvec = np.sqrt(varia)*MAT.matrix(np.random.standard_normal(tlena))
                    else:
                        ranvec = MAT.zeros((1,tlena))
                count = count + 1
            ranvec = MAT.vstack((ranvec,MAT.zeros((nendo,tlena))))
            # Save the random shocks vector
            self.shockvec = COP.deepcopy(ranvec)
        else:
            ranvec = COP.deepcopy(shockvec)

        x_one_m1 = ranvec[:,0]
        y_one_0 = ff*x_one_m1
        y_one_m1 = MAT.zeros(y_one_0.shape)
        x_one_0 = pp*x_one_m1
        if self.oswitch:
            nother = self.nother
            o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

        x_one_c = COP.deepcopy(x_one_m1)
        y_one_c = COP.deepcopy(y_one_0)
        if self.oswitch:
            o_one_c = COP.deepcopy(o_one_0)

        x_one = COP.deepcopy(x_one_c)
        y_one = COP.deepcopy(y_one_c)
        if self.oswitch:
            o_one = COP.deepcopy(o_one_c)

        for i1 in range(1,ranvec.shape[1],1):

            x_one_n = pp*x_one_c+ranvec[:,i1]
            y_one_n = ff*x_one_c
            if self.oswitch:
                o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

            x_one = MAT.hstack((x_one,x_one_c))
            y_one = MAT.hstack((y_one,y_one_n))
            if self.oswitch:
                o_one = MAT.hstack((o_one,o_one_n))

            x_one_c = x_one_n
            y_one_c = y_one_n

        # Throw away first 1000 obs
        x_one = x_one[:,1000:]
        y_one = y_one[:,1000:]
        if self.oswitch:
            o_one = o_one[:,1000:]

        self.sim_x_one = x_one
        self.sim_y_one = y_one
        self.insim = [self.sim_x_one,self.sim_y_one]
        if self.oswitch:
            self.sim_o_one = o_one
            self.insim = self.insim + [self.sim_o_one,]

    def irf(self,tlen,sntup):
        tlen = tlen + 1
        ncon = self.ncon
        nexo = self.nexo
        nendo = self.nendo
        tstates = self.tstates
        sigma = self.sigma
        ssigma = self.ssigma
        tlen = tlen
        if self.oswitch:
            numjl = self.numjl

        pp = self.P
        ff = self.F


        sposli=[]
        exoli = [x[1] for x in self.vardic['exo']['var']]
        # Check if names are valid
        for name in sntup:
            if name not in exoli:
                return 'Error: '+name+' is not a valid exoshock name for this model!'
        for name in sntup:
            sposli.append(exoli.index(name))
        # Expose spos choice to self for show_irf
        self.spos = (sntup,sposli)

        shock = MAT.zeros((nexo,1))
        sendo = MAT.zeros((nendo,1))

        for spos in sposli:
            shock[spos,0] = 1.0
        shock = MAT.vstack((shock,sendo))

        x_one_m1 = shock
        y_one_0 = ff*x_one_m1
        y_one_m1 = MAT.zeros(y_one_0.shape)
        x_one_0 = pp*x_one_m1
        if self.oswitch:
            nother = self.nother
            o_one_0 = numjl*MAT.vstack((x_one_0,y_one_0,x_one_m1,y_one_m1))

        x_one_c = COP.deepcopy(x_one_m1)
        y_one_c = COP.deepcopy(y_one_0)
        if self.oswitch:
            o_one_c = COP.deepcopy(o_one_0)

        x_one = COP.deepcopy(x_one_c)
        y_one = COP.deepcopy(y_one_c)
        if self.oswitch:
            o_one = COP.deepcopy(o_one_c)

        for i1 in range(tlen):
            x_one_n = pp*x_one_c
            y_one_n = ff*x_one_c
            if self.oswitch:
                o_one_n = numjl*MAT.vstack((x_one_n,y_one_n,x_one_c,y_one_c))

            x_one = MAT.hstack((x_one,x_one_c))
            y_one = MAT.hstack((y_one,y_one_n))
            if self.oswitch:
                o_one = MAT.hstack((o_one,o_one_n))

            x_one_c = x_one_n
            y_one_c = y_one_n

        # Throw away first observation
        self.irf_x_one = x_one[:,1:-1]
        self.irf_y_one = y_one[:,1:-1]
        self.inirf = [self.irf_x_one,self.irf_y_one]
        if self.oswitch:
            self.irf_o_one = o_one[:,2:]
            self.inirf = self.inirf + [self.irf_o_one,]

#------------------------currently not called and really not working, removed from init----------------
class FairTaylor:

    def __init__(self,algtype='type1',intup=None):
        self.algtype = algtype
        self.data = intup[0]
        self.param = intup[1]
        self.sstates_list = intup[2]
        self.vardic = intup[3]
        self.vardic2 = intup[4]
        self.modeq = intup[5]
        self.maxlags = intup[6]
        self.maxleads = intup[7]
        self.maxshocks = intup[8]
        self.maxinfosets = intup[9]
        self.shock_vars = intup[10]
        self.iid_vars = intup[11]
        self.vars = intup[12]
        self.re_var = intup[13]
        self.re_var2 = intup[14]
        self.sstate = intup[15]

        self.prepdat(self.data)
        self.prepmodeq(self.modeq)

    def prepdat(self,data=None):
        self.cur_pres = {}
        self.cur_past = {}
        self.cur_fut = {}
        for x1 in data.items():
            self.cur_pres[x1[0]] = x1[1][self.maxlags:len(x1[1])-self.maxleads]
            self.cur_past[x1[0]] = x1[1][0:-2]
            self.cur_fut[x1[0]] = x1[1][2:]

    def prepmodeq(self,modeq=None):
        re_var2 = self.re_var2
        self.modeq1 = COP.deepcopy(modeq)
        tmp_list = []
        i1 = 0
        for x1 in modeq:
            itera = re_var2.finditer(x1)
            for x2 in itera:
                if x2.group('svar'):
                    tmp_list = tmp_list + [[i1,x2.group('svar'),x2.start(),x2.end()],]
            i1 = i1 + 1
        tmp_list.reverse()
        if tmp_list:
            for x2 in tmp_list:
                i1 = x2[0]
                expr = x2[1]
                estart = x2[2]
                eend = x2[3]
                x1 = self.modeq1[i1][0:estart]+expr.capitalize()+\
                   '_bar'+self.modeq1[i1][eend:]
                self.modeq1[i1] = x1
        self.modeq2 = self.modeq1[0:len(modeq)-self.maxshocks]


        self.modeq3 = COP.deepcopy(self.modeq2)
        # Create varlist from vardic to establish fixed order
        self.varlist = self.vardic.items()
        # Create ordered vars and varnames
        varnames = []
        variables = []
        sub_list1 = []
        sub_list2 = []
        tmp_index = {}
        for x1 in self.varlist:
            varnames = varnames + [x1[1],]
            variables = variables + [x1[0],]
        self.varnames = varnames
        self.variables = variables
        re_var2 = self.re_var2
        i1 = 0
        i2 = 0
        for x1 in self.modeq2:
            itera = re_var2.finditer(x1)
            for x2 in itera:
                if x2.group('vvar') and not x2.group('exp') and not x2.group('lagt'):
                    tmp_match = variables.index(x2.group('var'))
                    if not tmp_index.has_key(x2.group('var')):
                        tmp_index[x2.group('var')] = (i2,x2.group(0))
                        i2 = i2 + 1
                    sub_list1 = sub_list1 + [['invar['+str(tmp_index[x2.group('var')][0])+']',
                                  x2.start(),x2.end(),x2.group(0),tmp_match,i1],]

            sub_list1.reverse()
            if sub_list1:
                for x4 in sub_list1:
                    str_tmp = x4[0]
                    pstart = x4[1]
                    pend = x4[2]
                    self.modeq3[i1] = self.modeq3[i1][0:pstart]+str_tmp+self.modeq3[i1][pend:]
                sub_list2 = sub_list2 + [sub_list1,]
                sub_list1 = []
            i1 = i1 + 1
        self.index = tmp_index
        # Create inverted index2
        index2 = {}
        for x1 in self.index.items():
            index2[x1[1][0]] = (x1[0],x1[1][1])
        self.index2 = index2
        return sub_list2

    def makecur(self,tpos,curlag,curpres,curfut):
        lag_dic={}
        pres_dic={}
        fut_dic={}
        all_dic={}
        for x1 in curlag.items():
            if float(x1[1][tpos].data) > 0:
                lag_dic[self.vardic2[x1[0]]+'_L'] = float(x1[1][tpos].data)
            else:
                lag_dic[self.vardic2[x1[0]]+'_L'] = 0.01
        for x1 in curpres.items():
            if float(x1[1][tpos].data) > 0:
                pres_dic[self.vardic2[x1[0]]+'_N'] = float(x1[1][tpos].data)
            else:
                pres_dic[self.vardic2[x1[0]]+'_N'] = 0.01
        for x1 in curfut.items():
            if float(x1[1][tpos].data) > 0:
                fut_dic['INFO_N_'+self.vardic2[x1[0]]+'_F'] = float(x1[1][tpos].data)
            else:
                fut_dic['INFO_N_'+self.vardic2[x1[0]]+'_F'] = 0.01
        all_dic.update(lag_dic)
        all_dic.update(pres_dic)
        all_dic.update(fut_dic)
        self.alldic =all_dic
        return all_dic

    def solveone(self,tpos,curlag,curpres,curfut):

        #Define the function to be handed over
        #to fsolve

        def func(invar):
            locals().update(self.param)
            locals().update(self.sstate[0])
            locals().update(self.makecur(tpos,self.cur_past,self.cur_pres,self.cur_fut))
            fdot1 = np.zeros((len(self.modeq3)),'d')
            i1=0
            for x in self.modeq3:
                fdot1[i1] = eval(x)
                i1 = i1 + 1
            return fdot1


        #def func(invar):
            #locals().update(self.param)
            #locals().update(self.sstate[0])
            #locals().update(self.makecur(tpos,self.cur_past,self.cur_pres,self.cur_fut))
            #fdot1 = np.zeros((len(self.modeq3)),'d')
            #i1=0
            #for x in self.modeq3:
                #fdot1[i1] = abs(eval(x))
                #i1 = i1 + 1
            #return sum(fdot1)

        #def fprime(invar):
            #return np.array(([1,-1]))


        # Define the initial values and
        # start the non-linear solver
        init_val = len(self.index)*[0,]
        for x1 in self.index.items():
            varia = self.vardic[x1[0]]['var']
            init_val[x1[1][0]] = float(curpres[varia][tpos].data)

        # Determine bounds
        #in_bounds = len(self.index)*[(0.01,None),]

        (output,infodict,ier,mesg) = O.fsolve(func,init_val,full_output=1)

        #(output,f_out,d_out) = O.fmin_l_bfgs_b(func,init_val,
                                                            #fprime,approx_grad=True,bounds=in_bounds)


        # Attach the outputs of the solver as attributes
        fsout={}
        self.output = output
        i1 = 0
        for x1 in self.index2.items():
            fsout[x1[1][1]] = (output[x1[0]],x1[1][0])
            i1 = i1 + 1
        self.fsout = fsout
        self.output = output
        self.infodict = infodict
        self.ier = ier
        self.mesg = mesg

        return fsout

    def solveall(self):
        # Calculate initial convergence criterion
        criterion = 10
        loop_c = 0
        while criterion > 0.01:
            cur_past_t = {}
            cur_pres_t = {}
            cur_fut_t = {}
            for x1 in self.vardic2.items():
                cur_past_t[x1[0]] = np.zeros(len(self.cur_pres.items()[0][1]))
                cur_past_t[x1[0]][0] = float(self.cur_past[x1[0]][0].data)
                cur_pres_t[x1[0]] = np.zeros(len(self.cur_pres.items()[0][1]))
                cur_fut_t[x1[0]] = np.zeros(len(self.cur_pres.items()[0][1]))
                cur_fut_t[x1[0]][-1] = float(self.cur_fut[x1[0]][-1].data)
            for i1 in range(0,len(self.cur_pres.items()[0][1])):
                out_tmp = self.solveone(i1,self.cur_past,self.cur_pres,self.cur_fut)
                for x2 in out_tmp.items():
                    # Standard
                    if self.algtype == 'type1':
                        cur_pres_t[self.vardic[x2[1][1]['var']]][i1] = x2[1][0]
                    # Fast-Gauss Seidel
                    elif self.algtype == 'type2':
                        cur_pres_t[self.vardic[x2[1][1]['var']]][i1] = x2[1][0]
                        try:
                            self.cur_past[self.vardic[x2[1][1]['var']]][i1+1] = x2[1][0]
                        except:
                            pass

            # Recalculate convergence criterion
            cr_list = [0,]*len(self.vardic)
            va_list = [0,]*len(self.vardic)
            for x1,x2 in zip(self.vardic.items(),range(0,len(self.vardic),1)):
                cr_list[x2] = cr_list[x2] + sum(abs(cur_fut_t[x1[1]][0:-1]-cur_pres_t[x1[1]][1:]))
                va_list[x2] = x1[1]
            criterion = sum(cr_list)
            loop_c = loop_c + 1
            print ('Loop: ',loop_c)
            print ('Var:  ',va_list)
            print ('Crit: ',cr_list)


            for x1 in self.vardic2.items():
                cur_past_t[x1[0]][1:] = cur_pres_t[x1[0]][:-1]
                cur_fut_t[x1[0]][0:-1] = cur_pres_t[x1[0]][1:]

            for x1,x2,x3 in zip(cur_past_t.items(),cur_pres_t.items(),cur_fut_t.items()):
                self.cur_past[x1[0]] = TSS.time_series\
                    (x1[1],start_date=TSS.Date(freq=self.cur_past[x1[0]].start_date.freqstr,
                                   year=self.cur_past[x1[0]].start_date.year,
                                   quarter=self.cur_past[x1[0]].start_date.quarter),
                                   freq=self.cur_past[x1[0]].start_date.freqstr)

                self.cur_pres[x1[0]] = TSS.time_series\
                    (x2[1],start_date=TSS.Date(freq=self.cur_pres[x1[0]].start_date.freqstr,
                                   year=self.cur_pres[x1[0]].start_date.year,
                                   quarter=self.cur_pres[x1[0]].start_date.quarter),
                                   freq=self.cur_pres[x1[0]].start_date.freqstr)

                self.cur_fut[x1[0]] = TSS.time_series\
                    (x3[1],start_date=TSS.Date(freq=self.cur_fut[x1[0]].start_date.freqstr,
                                   year=self.cur_fut[x1[0]].start_date.year,
                                   quarter=self.cur_fut[x1[0]].start_date.quarter),
                                   freq=self.cur_fut[x1[0]].start_date.freqstr)
