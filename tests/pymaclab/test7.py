# Import the pymaclab module into its namespace, also import os module
import pymaclab as pm
from pymaclab.modfiles import models

# Define the ssidic of initial guesses or starting values
alpha = 0.36
tau   = 5.0
g     = 0.005
delta = 0.025
betta = 1.014
bettaa = betta*(1+g)**(-tau)
ssidic = {}

ssidic['k_bar']     = 30.0                                                  #1
ssidic['inv_bar']   = 1.0                                                   #2
ssidic['c_bar']     = 2.0                                                   #3
ssidic['mu_bar']    = 8.0                                                   #4

# Pick first 10 FOC equations to leave out exogenous law of motion
focs = range(4)

# Instantiate a new DSGE model instance like so
asset = pm.newMOD(models.stable.jermann98,mesg=True,use_focs=focs,ssidic=ssidic)
