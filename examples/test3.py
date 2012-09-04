import pymaclab as pm
from pymaclab.modfiles import models

rbc1 = pm.newMOD(models.rbc1,mesg=True,ncpus=4)
# rbc2 = pm.newMOD(models.rbc2,mesg=True,ncpus=4)
merz = pm.newMOD(models.merz,mesg=True,ncpus=4)
