.. index:: cloud; sphinx theme, sphinx theme; cloud

=======================
PyMacLab Tutorial
=======================

Simulating DSGE models
======================

*Introduction*

  In the previous tutorial we discovered the general structure of the PyMacLab DSGE model instance and saw how this programming approach lent
  itself well to the idea of inspecting and exploring the instantiated models' current state, summarized by its data fields and supplied
  instance methods equipping them with functionality. This section finally discusses how the library equips DSGE model instances with methods
  which make use of the models' computed Jacobian and Hessian, which are evaluated at the models' numerical steady state. As a short reminder,
  we may recall here that it is often this step of obtaining a steady state can prove difficult in the case of a number of well-known models.
  Notwithstanding, for the remainder of this section we will assume that a steady state has successfully been attained and that the model's
  Jacobian and Hessian have been computed. Let's first start with our usual setup of lines of code entered into an IPython shell:

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: from pymaclab.modfiles import models

    # Instantiate a new DSGE model instance like so
    In [4]: rbc1 = pm.newMOD(models.rbc1)

  Recall that instantiating a DSGE model without any additional parameters means that an automatic attempt of finding the steady state is being
  made and if successfully found the model is already solved dynamically using a preferred (1st-order approximate) method. This means that you 
  now have all the information you need to simulate the model as well as to generate impulse-response functions and plot them. Let's focus on
  the simulations first, they are being generate using the following command:

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: from pymaclab.modfiles import models

    # Instantiate a new DSGE model instance like so
    In [3]: rbc1 = pm.newMOD(models.rbc1)

    # Now simulate the model
    In [4]: rbc1.modsolvers.forkleind.sim(200,('productivity'))

  This would simulate the model for 200 time periods using only iid shocks hitting the law of motion of the total factor productivity term in the
  model. Notice that here RBC1 is a very simple model only containing one structural shock, but more complicated models may possess more than one
  exogenous state variable. In that case, if you called instead:

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: from pymaclab.modfiles import models

    # Instantiate a new DSGE model instance like so
    In [3]: rbc1 = pm.newMOD(models.rbc1)

    # Now simulate the model
    In [4]: rbc1.modsolvers.forkleind.sim(200)

  The model would be simulated using all shocks (exogenous state variables) specified in the model. However, since RBC1 only contains one shock,
  the two variants shown here of simulating the model would yield the same results as "productivity" is is the only exogenous state available here
  anyway. We can then also graph the simulation to get a better understanding of the model by running the command:

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: from pymaclab.modfiles import models

    # Also import matplotlib.pyplot for showing the graph
    In [3]: from matplotlib import pyplot as plt

    # Instantiate a new DSGE model instance like so
    In [4]: rbc1 = pm.newMOD(models.rbc1)

    # Now solve and simulate the model
    In [5]: rbc1.modsolvers.forkleind.solve()
    In [6]: rbc1.modsolvers.forkleind.sim(200)

    # Plot the simulation and show it on screen
    In [7]: rbc1.modsolvers.forkleind.show_sim(('output','consumption'))
    In [8]: plt.show()

  This produces the following nice graph. Notice that you must specify the variables to be graphed and all simulated data is filtered according
  to the argument passed to each variable in the model file. So hp gave hp-filtered data while bk gave Baxter-King-filtered data.

  .. plot:: ../../pymaclab/examples/test4.py


  That was nice and simple, was it not? Notice that filtered simulations are always stored in data fields which means that statistics such as
  correlations at leads and lags can easily be computed as well.

Generating impulse-response functions
=====================================

*Introduction*

  Dynamic solutions obtained to first-order approximated DSGE models using the method of perturbations have a lot in common with standard
  Vector Autoregression models commonly used in applied Macroeconometrics. This in turn implies that solved DSGE models can be described using
  so-called impulse-responses functions (also abbreviated as IRFs) or impulse-response graphs which show how the solved model responds to a
  one-off shock to a particular exogenous state variable. In PyMacLab this can easily be achieved as follows:

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: from pymaclab.modfiles import models

    # Also import matplotlib.pyplot for showing the graph
    In [3]: from matplotlib import pyplot as plt

    # Instantiate a new DSGE model instance like so
    In [4]: rbc1 = pm.newMOD(models.rbc1)

    # Now solve and simulate the model
    In [5]: rbc1.modsolvers.forkleind.solve()
    In [6]: rbc1.modsolvers.forkleind.irf(100,('productivity',))

    # Plot the simulation and show it on screen
    In [7]: rbc1.modsolvers.forkleind.show_irf(('output','consumption'))
    In [8]: plt.show()

  This produces the following nice graph. Notice that here the shock to total productivity has been normalized to 100%.

  .. plot:: ../../pymaclab/examples/test5.py