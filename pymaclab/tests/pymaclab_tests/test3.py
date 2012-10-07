import pymaclab as pm
from pymaclab.modfiles import models

def test3():
    print "Now Loading and differentiating model rbc1..."
    rbc1 = pm.newMOD(models.stable.rbc1_res,mesg=True,ncpus='auto')
    print
    print "Now loading and differentiating model rbc2..."
    rbc2 = pm.newMOD(models.stable.rbc2,mesg=True,ncpus='auto')
    print
    print "Now loading and differentiating model merz..."
    merz = pm.newMOD(models.stable.merz,mesg=True,ncpus='auto')
    print
    print "Now loading and differentiating model cee..."
    cee  = pm.newMOD(models.testing.cee,mesg=True,ncpus='auto')
