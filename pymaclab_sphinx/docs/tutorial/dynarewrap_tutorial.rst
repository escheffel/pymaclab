.. index:: tutorial; DSGE instance; simulation; impulse-response; plotting; filtering

.. raw:: latex

   \newpage

Tutorial 8 - PyMacLab-to-Dynare++ Converter
===========================================

Introduction
------------

  PyMacLab comes shipped with a very useful dynare++ wrapper/translator which essentially possesses two types of functionality. First of all,
  it can translate a PyMacLab-formatted DSGE model files into a valid dynare++ model file. Secondly, it internally can pass this translated
  model file to dynare++, instruct it to solve the model, and then pass all of the solution objects back to PyMacLab where they get stored in
  the solution branch ``model.modsolvers.dynare.``.

  .. note::

    Notice that in order for PyMacLab's dynare++ wrapper/translator to function correctly, in a standard \*nix environment you *need* to have a
    functioning (compiled or otherwise obtained) dynare++ binary (executable) program file in your local PATH. This means that dynare++ on your
    system can be called from any location inside you shell. Also important is that the Python library ``Mako`` is installed which is a templating
    library used in generating the dynare++ conformable model file. This library should however be automatically installed when you installed
    PyMacLab on your system.
  
  The most basic functionality of the wrapper/translator works as follows calling the method ``model.mk_dynare()`` without any additional arguments:

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: from pymaclab.modfiles import models

    # Instantiate a new DSGE model instance like so
    In [4]: rbc1 = pm.newMOD(models.stable.rbc1_focs)
    
    # Use the dynare++ translator/wrapper as follows, by calling
    In [5]: rbc1.mk_dynare()
    
    # If the model solved, then check the contents of the dynare branch
    In [6]: dir(rbc1.modsolvers.dynare) 
    ['FF',
     'PP',
     '__class__',
     '__delattr__',
     '__dict__',
     '__doc__',
     '__format__',
     '__getattribute__',
     '__hash__',
     '__init__',
     '__module__',
     '__new__',
     '__reduce__',
     '__reduce_ex__',
     '__repr__',
     '__setattr__',
     '__sizeof__',
     '__str__',
     '__subclasshook__',
     '__weakref__',
     'attdic',
     'dyn_g_0',
     'dyn_g_1',
     'dyn_i_R',
     'dyn_i_c',
     'dyn_i_eps',
     'dyn_i_k',
     'dyn_i_y',
     'dyn_i_z',
     'dyn_irfm_eps_mean',
     'dyn_irfm_eps_var',
     'dyn_irfp_eps_mean',
     'dyn_irfp_eps_var',
     'dyn_mean',
     'dyn_nboth',
     'dyn_nforw',
     'dyn_npred',
     'dyn_nstat',
     'dyn_shocks',
     'dyn_ss',
     'dyn_state_vars',
     'dyn_steady_states',
     'dyn_vars',
     'dyn_vcov',
     'dyn_vcov_exo',
     'dynorder',
     'dynsorder',
     'modfile',
     'sstate']
     
  Notice that *all* of dynare++'s computed solution objects (which are mostly matrices containing various types of information characterizing the
  solution to the model) have been attached to this solution branch as standard numpy matrices/scalars, which can then be inspected inside a normal
  scientific Python stack environment. The *only* two objects which are being computed in a *derived* fashion using dynare++'s results are the matrices
  ``model.modsolvers.dynare.PP`` and ``model.modsolvers.dynare.FF`` which are the implied solutions from dynare++'s output in the solution format which
  get internally computed in PyMacLab using the ``model.modsolvers.forkleind.solve()`` method instead.
  
  This is useful as they can be directly compared and so checked for consistency between the different solution methods. Also important is the attribute
  ``model.modsolvers.dynare.modfile`` which is a string containing the PyMacLab-to-dynare++ translated model file which was passed to dynare++ in order
  to find the solution. Notice that in test runs it has not been uncommon to encounter situations in which PyMacLab's internal algorithms successfully
  computed a model's steady state, while dynare++ failed in doing so.
  
  Finally, bear in mind that dynare++'s standard solution method uses an approximation around the *level* of the model's steady state, and not
  the *logarithm*! This means that when you wish to compare solutions from dynare++ with those computed internally using PyMacLab's
  ``model.modsolvers.forkleind.solve()`` method, to take one specific example, then first make sure to *omit* the [log] qualifier in the variable section
  for each variable defined in the PyMacLab DSGE model file. That way you are instructing PyMacLab to also employ approximations around the level and not
  the logarithm of the model's SS.

The mk_dynare() solution method's options
-----------------------------------------

  Sometimes it may be useful to instruct PyMacLab *not* to solve the model using dynare++ but instead just to generate the dynare++ conformable file for
  inspection by the model-builder. This can be done as follows:
  
  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: from pymaclab.modfiles import models

    # Instantiate a new DSGE model instance like so
    In [4]: rbc1 = pm.newMOD(models.stable.rbc1_focs)
    
    # Use the dynare++ translator/wrapper as follows, by calling
    In [5]: rbc1.mk_dynare(fpath='./test.mod')
    
  In this case, no finding of a model solution would be attempted, but instead the PyMacLab-to-dynare++ translated model file would be copied into the
  location passed to the `fpath=` option. In the above example we are saving the translated dynare++ file into the current directory using the
  filename `test.mod`. For the above model file the generated dynare++ model file would look as follows:
  
  .. sourcecode:: ipython
  
    var k, c, y, R, z;
    varexo eps;

    parameters z_bar, psi, sigma_eps, betta, eta, rho, delta, c_bar, R_bar, k_bar, y_bar;

    z_bar = 1.0;
    psi = 0.95;
    sigma_eps = 0.052;
    betta = 0.990099009901;
    eta = 2.0;
    rho = 0.36;
    delta = 0.025;
    c_bar = 2.0;
    R_bar = 1.01;
    k_bar = 30.0;
    y_bar = 3.40222985854;


    model;
    0=y - (k-(1-delta)*k(-1)) - c;
    0=betta * ((((c(+1)^(-eta))))/((c^(-eta)))) * R(+1) - 1;
    0=R - (1+(z(-1)*k(-1)^(-1 + rho)*rho)-delta);
    0=y - (z(-1)*k(-1)^rho);
    0=log(z) - psi*log(z(-1)) - eps;
    end;


    initval;
    R = 1.01; 
    c = 2.0; 
    k = 30.0; 
    y = 3.40222985854; 
    z = 1.0; 
    z = 1.0;
    end;


    order = 1;

    vcov = [ 0.002704 ];

  Notice here that the starting values passed to dynare++'s non-linear steady state solver are already the solutions found by the non-linear solution
  algorithm used internally within the PyMacLab library itself. This particular function is useful in the sense that it may help model builders
  identify bugs or translation mistakes which may have occurred when converting the PyMacLab DSGE model to a dynare++ conformable format. Finally the full
  set of arguments which can be passed to the ``mk_dynare()`` method are given by ``model.mk_dynare(order=1,centralize=False,fpath=None,focli=None)``.
  
  The first argument `order=` is clearly used to specify the order of
  approximation used, `centralize=` can be used to specify whether the approximation should be computed around the deterministic or the stochastic
  (or risky) steady state [#f1]_, while `focli=` accepts a list or tuple of index numbers in order to pick out a subset of specific equations from PyMacLab's
  DSGE model file's declaration of the nonlinear system of equations of the model's first-order conditions of optimality. `fpath=` has already been
  discussed in the above and is an option which can be used in order to output the dynare++ conformable model file to a specific location in your
  computer's file system.


.. rubric:: Footnotes

.. [#f1] Dynare++ usually as default behaviour computes approximations around the stochastic or risk steady state. In order to be able to make comparisons
    between dynare++'s solutions and those obtained using most of the solution algorithms which are available internally to PyMacLab, it is important to
    pass the `centralize=False` option, which is the standard choice for ``mk_dynare``'s method call. That way dynare++ computes the solutions around the
    deterministic steady state of the model.