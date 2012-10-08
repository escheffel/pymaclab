# Import the pymaclab module into its namespace, also import os module
import pymaclab as pm
from pymaclab.modfiles import models
from copy import deepcopy

def test8():

    # Test the steady state method in which USE_FOCS was declared inside the model
    rbc_focs = pm.newMOD(models.stable.rbc1_focs,mesg=True)
    # Did it work?
    assert 'sstate' in dir(rbc_focs)
    
    # Now copy the sstate dictionary and test the steady state method in which it is passed from outside the model
    statedic = deepcopy(rbc_focs.sstate)
    rbc_ext = pm.newMOD(models.stable.rbc1_extss,mesg=True,sstate=statedic)
    # Did it work?
    assert 'sstate' in dir(rbc_ext)
    
