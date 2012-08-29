===============================================
PyMacLab - Python Macroeconomics Laboratory
===============================================

Introduction
============

Description
-----------
  :mod:`!PyMacLab` is the Python Macroeconomics Laboratory which currently primarily serves the purpose
  of providing a convenience framework written in form of a `Python <http://www.python.org/>`_ library with the ability to solve
  non-linear DSGE models. At the time of writing this, the library supports solving DSGE models
  using 1st and 2nd order perturbation methods which are computed around the steady state. In particular, the library provides wrapper
  function for `Paul Klein's <http://paulklein.ca/newsite/start/start.php>`_ 1st-order accurate method based on the
  Schur Decomposition :cite:`Kle:2000` as well a more recently published method by the same author (co-authored with Paul Gomme)
  which provides 2nd-order accurate solutions without using Tensor Algebra :cite:`GomKle:2012` (using the Magnus and Neudecker
  1999 definition of the Hessian matrix). :mod:`!PyMacLab` however possesses the added advantage of being equipped with an advanced
  model file parser module, similar to the one available in `Dynare <http://www.dynare.org>`_, which automates cumbersome and error-prone
  (log-)linearization by hand. :mod:`!PyMacLab` is also written entirely in Python, which makes it free and incredibly flexible to use and extend.

Features at a Glance
--------------------
  * No log-linearisation by hand required, done automatically py parsing a DSGE model file.
  * All solutions based on analytical/symbolical computation of Jacobian and Hessian of models using `Sympycore <http://code.google.com/p/sympycore/>`_.
  * DSGE models are Python DSGE class instances, treat them as if they were data structures, pass them around, copy them, stack them into arrays,
    and work with many of them simultaneously!
  * Loop over a DSGE model instance thousands of times to alter the parameter space, each time re-computing the solution.
  * Choose from closed form or non-linear steady state solvers.
  * Choose from a number of tried and tested perturbation methods, such as Klein's 1st order accurate and Klein & Gomme's 2nd order accurate methods.
  * Solving models is as fast as using optimized compiled C or Fortran code, most CPU time spent parsing the DSGE model file.
  * DSGE model examples are provided, including very complex ones such as the one based on Christiano, Eichenbaum and Evans (2001) :cite:`ChrEicEva:2005`.
  * Benefit from a large and growing set of convience methods to simulate models and graph filtered simulated series as well as impulse-response functions.
  * Enjoy the power, flexibility and extensibility of the Python programming language.
  * PyMacLab is free as in freedom and distributed under a `Apache 2.0 license <http://www.apache.org/licenses/LICENSE-2.0.html>`_.

  PyMacLab is currently known to work well but continues to mature. If you have used PyMacLab already and spotted some bugs or felt that some other important
  features are missing, you can head over to the library's `Github <https://github.com/escheffel/pymaclab/>`_ repository to submit an Issue item. We are
  eager to receive your feedback!

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


Getting Started
---------------
:doc:`A basic DSGE tutorial <started_tutorial>`
    This section provides a succinct tutorial on how to use PyMacLab to work with DSGE models.
:doc:`A description of all template DSGE models <started_allmodels>`
    This section gives detailed descriptions of all of the template DSGE models which come supplied with PyMacLab, including the derivation of all necessary algebraic results.


API Documentation
------------------
:doc:`api_doc`
    The auto-generated documentation of pymaclab's main modules and classes

Reference
---------
:doc:`bibliography`
    requirements and installations instructions

:doc:`history`
    history of current and past releases

Download & Installation
=======================

PyMacLab is known to work with any of Python version greater than or equal to 2.4 and smaller than 3.0.
In the future we will consider implementing a compatibility branch for versions of Python greater
than or equal to 3.0, once all core dependencies are known to have been migrated as well.

Option 1
----------
You can download the source code of PyMacLab right here. Alternatively, PyMacLab is also hosted at PyPI and
can be installed in the usual way by executing the command inside a Linux shell using ``pip``::

    sudo pip install pymaclab

Option 2
---------
Otherwise get the latest source code compressed as a tarball here:

`pymaclab.tar.gz <https://github.com/escheffel/pymaclab/tarball/v0.8>`_

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

Dependencies
-------------
Proper functioning of PyMacLab depends on a number of additional Python libraries already being installed on
your system, such as `Numpy <http://numpy.scipy.org/>`_, `Scipy <http://www.scipy.org/>`_,
`Sympycore <http://code.google.com/p/sympycore/>`_, `Sympy <http://sympy.org/en/index.html>`_ and
`scikits.timeseries <http://pytseries.sourceforge.net/>`_. All of these are great libraries by themselves and
should be checked out by any serious scientist interested in doing work in Python.

Also, if you want to enjoy a Matlab-style interactive environment in which to execute and inspect DSGE and other
data structures, you'd be hard-pressed to pass over the brilliant and now extra features-ladden
`IPython <http://ipython.org/>`_. When downloading and installing pymaclab using ``pip`` all of these dependencies
should be installed automatically for you, if they are not already present on your system.

To use some convience plotting methods included in PyMacLab, you need to have Python's most advanced plotting
library installed, which is called `Matplotlib <http://matplotlib.sourceforge.net/>`_. Besides being indispensable
for any scientist working with Python for graphical analysis, it is also used to quickly produce plots of simulated
solved DSGE models as well as impulse response functions.

Credits & Thanks
================
Thanks and kudos must go to all members of the Python scientic community without whose efforts projects like PyMacLab
would be much harder to implement. We are all standing on the shoulders of giants! Special thanks go to
Eric Jones, Travis Oliphant and Pearu Peterson, the leading coders of the `Numpy/Scipy <http://www.scipy.org>`_ Suite
which PyMacLab heavily makes use of, as well as `Skipper Seabold <https://github.com/jseabold>`_, lead coder of another
unique and outstanding Python library, `Statsmodels <http://statsmodels.sourceforge.net/>`_, who has kindly helped me
clean up some of the rough edges of my code.

Online Resources
================

    .. rst-class:: html-plain-table

    ====================== ===================================================
    Author Homepage:       `<http://www.ericscheffel.com>`_
    Github Homepage:       `<https://github.com/escheffel/pymaclab>`_
    Scipy Homepage:        `<http://www.scipy.org>`_
    Download & PyPI:       `<http://pypi.python.org/pypi/pymaclab>`_
    ====================== ===================================================
