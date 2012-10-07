# Import the pymaclab module into its namespace, also import os module
import pymaclab as pm
from pymaclab.modfiles import models

# Also import matplotlib.pyplot for showing the graph
from matplotlib import pyplot as plt

def test4():
    # Instantiate a new DSGE model instance like so
    rbc1 = pm.newMOD(models.stable.rbc1_res)

    # Now simulate the model
    rbc1.modsolvers.forkleind.solve()
    rbc1.modsolvers.forkleind.sim(200)

    # Plot the simulation and show it on screen
    rbc1.modsolvers.forkleind.show_sim(('output','consumption'))
    #plt.show()
