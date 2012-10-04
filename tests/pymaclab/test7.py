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

# We cheat here and use values from the steady state which was already solved for
ssidic['k_bar']     = 36.5                                                  #1
ssidic['y_bar']     = 3.644                                                 #2
ssidic['inv_bar']   = 1.089                                                 #3
ssidic['c_bar']     = 2.555                                                 #4
ssidic['mu_bar']    = 8.0                                                   #5
ssidic['w_bar']     = 2.332                                                 #6
ssidic['div_bar']   = 0.222                                                 #7
ssidic['rf_bar']    = 1.011                                                 #8
ssidic['r_bar']     = 1.011                                                 #9
ssidic['eqp_bar']   = 0.0                                                   #10

# Pick first 10 FOC equations to leave out exogenous law of motion
focs = range(10)

# Instantiate a new DSGE model instance like so
asset = pm.newMOD(models.stable.jermann98,mesg=True,use_focs=focs,ssidic=ssidic)
