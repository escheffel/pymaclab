# Import the pymaclab module into its namespace, also import os module
import pymaclab as pm
from pymaclab.modfiles import models

# Also import matplotlib.pyplot for showing the graph
from matplotlib import pyplot as plt
from copy import deepcopy

# Instantiate a new DSGE model instance like so
rbc1 = pm.newMOD(models.stable.rbc1_res,mk_hessian=False)

# Now solve and simulate the model
rbc1.modsolvers.forkleind.solve()
rbc1.modsolvers.forkleind.sim(200)

# Plot the simulation and show it on screen
rbc1.modsolvers.forkleind.show_sim(('output','consumption'))
plt.show()

# Now save the shocks, by saving a clone or copy, instead of a reference
shockv = deepcopy(rbc1.modsolvers.forkleind.shockvec)

# Change the filterin assumption of output and consumption using the queued updater branch
rbc1.updaters_queued.vardic['con']['mod'][0][1] = 'hp'
rbc1.updaters_queued.vardic['con']['mod'][1][1] = 'hp'
rbc1.updaters_queued.process_queue()

# Now we could run the simulation again, this time passing the randomly drawn shocks
rbc1.modsolvers.forkleind.solve()
rbc1.modsolvers.forkleind.sim(200,shockvec=shockv)

# Plot the simulation and show it on screen
rbc1.modsolvers.forkleind.show_sim(('output','consumption'))
plt.show()

# Change the filterin assumption of output and consumption using the queued updater branch
rbc1.updaters_queued.vardic['con']['mod'][0][1] = 'cf'
rbc1.updaters_queued.vardic['con']['mod'][1][1] = 'cf'
rbc1.updaters_queued.process_queue()

# Now we could run the simulation again, this time passing the randomly drawn shocks
rbc1.modsolvers.forkleind.solve()
rbc1.modsolvers.forkleind.sim(200,shockvec=shockv)

# Plot the simulation and show it on screen
rbc1.modsolvers.forkleind.show_sim(('output','consumption'))
plt.show()
