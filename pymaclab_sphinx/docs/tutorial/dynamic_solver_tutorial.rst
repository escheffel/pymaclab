.. index:: cloud; sphinx theme, sphinx theme; cloud

=======================
PyMacLab Tutorial
=======================

Exploring dynamic solution methods
==================================

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

  Since we have not supplied the instantiation call with the `initlev` argument, the model has been solved out completely, which includes the
  computation of a preferred dynamic solution, which in the current library's version is Paul Klein's 1st-order accurate method based on the
  Schur Decomposition. This means that one particular solution - if existing and found - has already been computed by default. Readers familar
  with this solution method will know that this method requires a partitioning of the evaluated model's Jacobian into two separate matrices, we
  call for simplicity `A` and `B`. They are computed by default and are for instance inspectable at ``rbc1.jAA`` and ``rbc1.jBB``.

  One of the great advantages of PyMacLab is that the Jacobian and Hessian are NOT computed numerically using the method of finite differences,
  but are calculated analytically using the special-purposes Python library Sympy. You can also inspect the analytical counterparts which have not yet
  been evaluated numerically at ``rbc1.jdic`` and ``rbc1.hdic``, where the latter reference actually refers to the 3-dimensional analytical
  Hessian. Notice that at this point in time, the analytical expression contain strange variable names, which are replacements for variables such
  as `c(t)` or `E(t)|z(t+1)`. The reason why PyMacLab internally operates with replacement definitions for variables has to do with the simple
  fact that computer algebra systems such as Sympy obviously interpret round brackets in a special mathematical way and thus are not permissible
  as part of variable names.

  to be finished...
