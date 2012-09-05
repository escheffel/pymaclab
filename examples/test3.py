import pymaclab as pm
from pymaclab.modfiles import models

rbc1 = pm.newMOD(models.rbc1,mesg=True,ncpus='auto')
rbc2 = pm.newMOD(models.rbc2,mesg=True,ncpus='auto')
merz = pm.newMOD(models.merz,mesg=True,ncpus='auto')
