# Import the pymaclab module into its namespace, also import os module
import pymaclab as pm
import pymaclab.modfiles.models as models
from copy import deepcopy
import numpy.matlib as MAT

# This test works through a number of linear models to check if they solve
def test9():
    
    #####################################
    ## STANDARD COOLEY HANSE CIA MODEL ##
    #####################################
    # Load and solve the manually linearized model
    rbc1 = pm.newMOD(models.abcs_rbcs.cooley_hansen_cia_linear,mesg=True)
    rbc1.modsolvers.pyuhlig.solve()
    rbc1.modsolvers.forklein.solve()
    
    # Load and solve the automatically linearized model
    rbc2 = pm.newMOD(models.abcs_rbcs.cooley_hansen_cia_cf,mesg=True)
    rbc2.modsolvers.forkleind.solve()
    
    # Check equivalence of steady states
    for keyo in rbc1.sstate.keys():
        if keyo in rbc2.sstate.keys():
            assert round(rbc1.sstate[keyo],5) == round(rbc2.sstate[keyo],5)
    
    # Check equivalence of results
    nexo = len(rbc2.vardic['exo']['var'])
    nendo = len(rbc2.vardic['endo']['var'])
    
    modlin1 = MAT.hstack((rbc1.modsolvers.pyuhlig.Q,rbc1.modsolvers.pyuhlig.P))
    modlin1 = [round(modlin1[0,i1],5) for i1 in range(modlin1.shape[1])]
    modnlin1 = rbc2.modsolvers.forkleind.P[-nendo:,:]
    modnlin1 = [round(modnlin1[0,i1],5) for i1 in range(modnlin1.shape[1])]
    print 'Comparison: Standard CIA model'
    print "Linear is: ",modlin1
    print '----------------------'
    print "Nonlinear is: ",modnlin1
    assert modlin1 == modnlin1
    
    modlin2 = MAT.hstack((rbc1.modsolvers.pyuhlig.S,rbc1.modsolvers.pyuhlig.R))
    modlin2 = [[round(modlin2[i2,i1],5) for i1 in range(modlin2.shape[1])] for i2 in range(modlin2.shape[0])]
    modnlin2 = rbc2.modsolvers.forkleind.F
    modnlin2 = [[round(modnlin2[i2,i1],5) for i1 in range(modnlin2.shape[1])] for i2 in range(modnlin2.shape[0])]
    print modlin2
    print '----------------------'
    print modnlin2    
    assert modlin2 == modnlin2
    

    #####################################################
    ## STANDARD COOLEY HANSE CIA MODEL WITH SEIGNORAGE ##
    #####################################################
    # Load and solve the manually linearized model
    rbc1 = pm.newMOD(models.abcs_rbcs.cooley_hansen_cia_seignorage_linear,mesg=True)
    rbc1.modsolvers.pyuhlig.solve()
    rbc1.modsolvers.forklein.solve()
    
    # Load and solve the automatically linearized model
    rbc2 = pm.newMOD(models.abcs_rbcs.cooley_hansen_cia_seignorage_cf,mesg=True)
    rbc2.modsolvers.forkleind.solve()
    
    # Check equivalence of steady states
    for keyo in rbc1.sstate.keys():
        if keyo in rbc2.sstate.keys():
            assert round(rbc1.sstate[keyo],5) == round(rbc2.sstate[keyo],5)
    
    # Check equivalence of results
    nexo = len(rbc2.vardic['exo']['var'])
    nendo = len(rbc2.vardic['endo']['var'])
    
    modlin1 = MAT.hstack((rbc1.modsolvers.pyuhlig.Q,rbc1.modsolvers.pyuhlig.P))
    modlin1 = [round(modlin1[0,i1],5) for i1 in range(modlin1.shape[1])]
    modnlin1 = rbc2.modsolvers.forkleind.P[-nendo:,:]
    modnlin1 = [round(modnlin1[0,i1],5) for i1 in range(modnlin1.shape[1])]
    print 'Comparison: Standard CIA model with seignorage'
    print "Linear is: ",modlin1
    print '----------------------'
    print "Nonlinear is: ",modnlin1
    assert modlin1 == modnlin1
    
    modlin2 = MAT.hstack((rbc1.modsolvers.pyuhlig.S,rbc1.modsolvers.pyuhlig.R))
    modlin2 = [[round(modlin2[i2,i1],5) for i1 in range(modlin2.shape[1])] for i2 in range(modlin2.shape[0])]
    modnlin2 = rbc2.modsolvers.forkleind.F
    modnlin2 = [[round(modnlin2[i2,i1],5) for i1 in range(modnlin2.shape[1])] for i2 in range(modnlin2.shape[0])]
    print modlin2
    print '----------------------'
    print modnlin2    
    assert modlin2 == modnlin2


    #####################################################################
    ## STANDARD COOLEY HANSE CIA MODEL WITH SEIGNORAGE AND CES UTILITY ##
    #####################################################################
    # Load and solve the manually linearized model
    rbc1 = pm.newMOD(models.abcs_rbcs.cooley_hansen_cia_seignorage_ces_linear,mesg=True)
    rbc1.modsolvers.pyuhlig.solve()
    rbc1.modsolvers.forklein.solve()
    
    # Load and solve the automatically linearized model
    rbc2 = pm.newMOD(models.abcs_rbcs.cooley_hansen_cia_seignorage_ces_cf,mesg=True)
    rbc2.modsolvers.forkleind.solve()
    
    # Check equivalence of steady states
    for keyo in rbc1.sstate.keys():
        if keyo in rbc2.sstate.keys():
            assert round(rbc1.sstate[keyo],5) == round(rbc2.sstate[keyo],5)
    
    # Check equivalence of results
    nexo = len(rbc2.vardic['exo']['var'])
    nendo = len(rbc2.vardic['endo']['var'])
    
    modlin1 = MAT.hstack((rbc1.modsolvers.pyuhlig.Q,rbc1.modsolvers.pyuhlig.P))
    modlin1 = [[round(modlin1[i2,i1],5) for i1 in range(modlin1.shape[1])] for i2 in range(modlin1.shape[0])]
    modnlin1 = rbc2.modsolvers.forkleind.P[-nendo:,:]
    modnlin1 = [[round(modnlin1[i2,i1],5) for i1 in range(modnlin1.shape[1])] for i2 in range(modnlin1.shape[0])]
    print 'Comparison: Standard CIA model with seignorage and CES utility'
    print "Linear is: ",modlin1
    print '----------------------'
    print "Nonlinear is: ",modnlin1
    assert modlin1 == modnlin1
    
    modlin2 = MAT.hstack((rbc1.modsolvers.pyuhlig.S,rbc1.modsolvers.pyuhlig.R))
    modlin2 = [[round(modlin2[i2,i1],5) for i1 in range(modlin2.shape[1])] for i2 in range(modlin2.shape[0])]
    modnlin2 = rbc2.modsolvers.forkleind.F
    modnlin2 = [[round(modnlin2[i2,i1],5) for i1 in range(modnlin2.shape[1])] for i2 in range(modnlin2.shape[0])]
    print modlin2
    print '----------------------'
    print modnlin2    
    assert modlin2 == modnlin2