===============================================
PyMacLab - Python Macroeconomics Laboratory
===============================================

Introduction
============

Description
-----------
  PyMacLab is the Python Macroeconomics Laboratory which currently primarily serves the purpose
  of providing a convenience framework written in form of a `Python <http://www.python.org/>`_ library with the ability to solve
  non-linear DSGE models. At the time of writing these words, the library supports solving DSGE models
  using 1st and 2nd order perturbation methods which are computed around the steady state. If you want to learn about PyMacLab
  as quickly as possible, skip reading this and instead start reading through the tutorial series available in the Documentation
  section to this site.

  The library provides wrapper functions for `Paul Klein's <http://paulklein.ca/newsite/start/start.php>`_ 1st-order
  accurate method based on the Schur Decomposition :cite:`Kle:2000` as well a more recently published method by the same author
  (co-authored with Paul Gomme) which computes 2nd-order accurate solutions without using Tensor Algebra :cite:`GomKle:2012`
  (using the Magnus & Neudecker 1999 definition of the Hessian matrix). PyMacLab possesses the added advantage of being equipped with
  an advanced model file parser module, similar to the one available in `Dynare <http://www.dynare.org>`_, which automates cumbersome
  and error-prone (log-)linearization by hand. PyMacLab is also written entirely in Python, is free and incredibly flexible to use and extend.

Features at a Glance
--------------------
  * No log-linearization by hand required, done automatically py parsing a DSGE model file.
  * All solutions based on analytical/symbolical computation of Jacobian and Hessian of models using `Sympy <http://www.sympy.org/>`_.
  * DSGE models are Python DSGE class instances, treat them as if they were data structures, pass them around, copy them, stack them into arrays,
    and work with many of them simultaneously!
  * Loop over a DSGE model instance thousands of times to alter the parameter space, each time re-computing the solution.
  * Choose from closed form or non-linear steady state solvers or a combination of both.
  * Choose from a number of tried and tested perturbation methods, such as Klein's 1st order accurate and Klein & Gomme's 2nd order accurate methods.
  * Solving models is as fast as using optimized compiled C or Fortran code, expensive computation of Jacobian and Hessian employs parallelized multi-core CPU approach.
  * DSGE example models are provided, including very complex ones such as the one based on Christiano, Eichenbaum and Evans (2001) :cite:`ChrEicEva:2005`.
  * Benefit from a large and growing set of convience methods to simulate models and graph filtered simulated series as well as impulse-response functions.
  * Use PyMacLab as a free Python library within a rich and rapidly evolving Python software ecosystem for scientists.
  * Enjoy the power, flexibility and extensibility of the Python programming language and the open-source transparency of PyMacLab.
  * PyMacLab is free as in freedom and distributed under a `Apache 2.0 license <http://www.apache.org/licenses/LICENSE-2.0.html>`_.

.. note::

    PyMacLab is currently known to work well but continues to mature. This documentation site is well under way but still work-in-progress.
    If you have used PyMacLab already and spotted some bugs or felt that some other important features are missing, you can head over to the
    library's `Github <http://github.com/escheffel/pymaclab/>`_ repository to submit an Issue item. We are currently in the process of adding
    more example DSGE model files (and eliminating mistakes in already existing ones). If you have used PyMacLab yourself and want to contribute
    your own DSGE model files we are happy to include them!

Documentation
=============

Introduction
------------

:doc:`What is PyMacLab? <pymaclab_intro>`
    This is a succinct introduction to PyMacLab including an explanation of its current features.
:doc:`Philosophy behind PyMacLab <pymaclab_philo>`
    Here I discuss the basic Philosophy behind PyMacLand and what it sets out to do now and in the near future.
:doc:`Why Macroeconomics in Python? <pymaclab_python>`
    In this section I touch upon the the pros and cons of doing Macroeconomics or scientific computing using Python in general.


Tutorial Series
---------------

1) :doc:`Basic DSGE tutorial <tutorial/started_tutorial>`
    Brief tutorial on how to use PyMacLab to work with DSGE models.
2) :doc:`PyMacLab DSGE instance tutorial <tutorial/dsge_instance_tutorial>`
    Succinct tutorial facilitating the understaning of the DSGE OOP data structure in PyMacLab.
3) :doc:`PyMacLab DSGE instance updater tutorial <tutorial/dsge_instance_updater_tutorial>`
    Succinct tutorial facilitating the understaning of the DSGE OOP data structure in PyMacLab.
4) :doc:`PyMacLab DSGE steady state solver tutorial <tutorial/steady_solver_tutorial>`
    This section finally shows how dynamic solution to the PyMacLab DSGE models are obtained.
5) :doc:`PyMacLab DSGE dynamic solver tutorial <tutorial/dynamic_solver_tutorial>`
    This section finally shows how dynamic solution to the PyMacLab DSGE models are obtained.
6) :doc:`PyMacLab DSGE simulation and plotting tutorial <tutorial/simirf_plotting_tutorial>`
    Short tutorial on using convenience functions for simulations, IRFs and plotting.
7) :doc:`Description of all template DSGE models <tutorial/started_allmodels>`
    Detailed description of all of the template DSGE models which come supplied with PyMacLab.


API Documentation
------------------

:doc:`api_doc`
    The auto-generated documentation of pymaclab's main modules and classes

Reference
---------

:doc:`bibliography`
    Reference list of academic articles related to the solution of DSGE models.

:doc:`history`
    History of current and past releases

Download & Installation
=======================

  PyMacLab is known to work with any of Python version greater than or equal to 2.4 and smaller than 3.0.
  In the future we will consider implementing a compatibility branch for versions of Python greater
  than or equal to 3.0, once all core dependencies are known to have been migrated as well.

Dependencies
-------------

  Proper functioning of PyMacLab depends on a number of additional Python libraries already being installed on
  your system, such as `Numpy <http://numpy.scipy.org/>`_, `Scipy <http://www.scipy.org/>`_,
  `Sympy <http://www.sympy.org>`_, `Matplotlib <http://matplotlib.sourceforge.net/>`_ and
  `scikits.timeseries <http://pytseries.sourceforge.net/>`_. All of these are great libraries by themselves and
  should be checked out by any serious scientist interested in doing work in Python. However, the installation
  of `Parallel Python <http://www.parallelpython.com/>`_ is also highly recommended, as this allows exploiting
  multi-core CPUs in the computation of DSGE models' Jacobian and Hessian significantly speeding up execution speed.

  Also, if you want to enjoy a Matlab-style interactive environment in which to execute and inspect DSGE and other
  data structures, you'd be hard-pressed to pass over the brilliant and now extra features-ladden
  `IPython <http://ipython.org/>`_. When downloading and installing pymaclab using ``pip`` all of these dependencies
  should be installed automatically for you, if they are not already present on your system.

  To use some convience plotting methods included in PyMacLab, you need to have Python's most advanced plotting
  library installed, which is called `Matplotlib <http://matplotlib.sourceforge.net/>`_. Besides being indispensable
  for any scientist working with Python for graphical analysis, it is also used to quickly produce plots of simulated
  solved DSGE models as well as impulse response functions. Following right below is a list of options users have to
  install PyMacLab on their Python-ready computers.

Option 1
----------

  You can download the source code of PyMacLab right here. Alternatively, PyMacLab is also hosted at
  `PyPI <http://pypi.python.org/pypi/pymaclab/>`_ and can be installed in the usual way by executing the
  command inside a Linux shell using ``pip``::

    sudo pip install pymaclab

Option 2
---------

  Otherwise get the latest source code compressed as a tarball here:

    `pymaclab.tar.gz <http://pypi.python.org/packages/source/p/pymaclab/pymaclab-0.88.1.tar.gz>`_

  And install it in the usual way by running in a Linux shell the command::

    sudo python setup.py install

Option 3
---------

  Alternatively, for the brave-hearted and bleeding-edge aficionados, they can also navigate over to our open
  Github repository where PyMacLab is currently being maintained, and clone the most up-to-date version and/or
  nightly build, by having git installed on your system and calling::

    git clone git://github.com/escheffel/pymaclab.git

  This will create a new folder called pymaclab containing the latest version of the source code as well as the
  installation script ``setup.py`` which you can then use in the usual way to install the module on your system.

Credit & Thanks
================

  Thanks must go to all members of the Python scientific community without whose efforts projects like PyMacLab
  would be much harder to implement. We are all standing on the shoulders of giants! Special thanks go to
  Eric Jones, Travis Oliphant and Pearu Peterson, the founding coders of the `Numpy/Scipy <http://www.scipy.org>`_ Suite
  which PyMacLab heavily makes use of.

  I would also like to give a special mention to `Skipper Seabold <http://github.com/jseabold>`_, lead coder of another
  unique and outstanding Python library, `Statsmodels <http://statsmodels.sourceforge.net/>`_, who has kindly helped me
  clean up some of the rough edges of my code. I would also like to thank colleagues at Nottingham University Business
  School China, especially `Gus Hooke <http://www.nottingham.edu.cn/en/business/staff/staffprofile/angushooke.aspx>`_
  and `Carl Fey <http://www.nottingham.edu.cn/en/business/people/staffprofile/carlfey.aspx>`_ for their kind support.

Online Resources
================

    .. rst-class:: html-plain-table

    ====================== ===================================================
    Author Homepage:       `<http://www.ericscheffel.com>`_
    Github Homepage:       `<http://github.com/escheffel/pymaclab>`_
    Scipy Homepage:        `<http://www.scipy.org>`_
    Download & PyPI:       `<http://pypi.python.org/pypi/pymaclab>`_
    ====================== ===================================================
