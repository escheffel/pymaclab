import pymaclab as pm
from pymaclab.modfiles import models

rbc1 = pm.newMOD(models.stable.rbc1_res,mesg=True,ncpus='auto')
rbc2 = pm.newMOD(models.stable.rbc2,mesg=True,ncpus='auto')
merz = pm.newMOD(models.stable.merz,mesg=True,ncpus='auto')
