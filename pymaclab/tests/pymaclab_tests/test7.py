# Import the pymaclab module into its namespace, also import os module
import pymaclab as pm
from pymaclab.modfiles import models

def test7():
    # Define the ssidic of initial guesses or starting values
    alpha      = 0.36
    chi        = 0.819
    xi         = 1.0/4.3
    g          = 0.005
    tau        = 5.0
    betta      = 1.014
    bettaa     = betta*(1+g)**(-tau)
    delta      = 0.025
    rho        = 0.95
    a1         = (g+delta)**(1.0/xi)
    a2         = (g+delta)/(1.0-xi)
    z_bar      = 1.0
    sigma_eps  = 0.01
    
    ssidic = {}
    ssidic['k_bar']     = 36.0                                                  #1
    ssidic['inv_bar']   = 1.0                                                   #2
    ssidic['c_bar']     = 2.5                                                   #3
    ssidic['mu_bar']    = 8.0                                                   #4

    # Pick first FOC equations to leave out exogenous law of motion
    focs = range(4)

    # Instantiate a new DSGE model instance like so
    asset = pm.newMOD(models.stable.jermann98_ext,mesg=True,use_focs=focs,ssidic=ssidic)
    
    return asset
