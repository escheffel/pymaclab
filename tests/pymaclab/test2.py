# Import the pymaclab module into its namespace, also import os module
import pymaclab as pm
from pymaclab.modfiles import models
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

# Instantiate a new DSGE model instance like so
rbc1 = pm.newMOD(models.stable.rbc1_res)

# Create an array representing a finely-spaced range of possible impatience values
# Then convert to corresponding steady state gross real interest rate value
betarr = np.arange(0.8,0.99,0.001)
betarr = 1.0/betarr


ss_capital = []
for betar in betarr:
    rbc1.updaters.paramdic['R_bar'] = betar
    rbc1.sssolvers.fsolve.solve()
    ss_capital.append(rbc1.sssolvers.fsolve.fsout['k_bar'])

# Create a nice figure
fig1 = plt.figure()
plt.grid()
plt.title('Plot of steady state physical capital against R\_bar')
plt.xlabel(r'Steady state gross real interest rate')
plt.ylabel(r'Steady State of physical capital')
plt.plot(betarr,ss_capital,'k-')
plt.show()
