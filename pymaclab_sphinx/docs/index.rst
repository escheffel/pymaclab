===============================================
PyMacLab - Python Macroeconomics Laboratory
===============================================

Introduction
============

Description
-----------
  PyMacLab is the Python Macroeconomics Laboratory which currently primarily serves the purpose
  of providing a convenience framework written in form of a `Python <http://www.python.org/>`_ library with the ability to solve
  non-linear DSGE models using a DSGE model class from which to instantiate instances. At the time of writing these words, the library supports solving DSGE models
  using 1st and 2nd order perturbation methods which are computed around the steady state. Apart from that, the library also contains two
  advanced macroeconometric classes, the VAR class and the FAVAR class which can be employed for empirical work or in combination with DSGE models
  in order to estimate instead of calibrate deep parameters. If you want to learn about PyMacLab as quickly as possible, skip reading this and
  instead start reading through the tutorial series available in the Documentation section to this site.

  The DSGE model class provides wrapper functions for `Paul Klein's <http://paulklein.ca/newsite/start/start.php>`_ 1st-order
  accurate method based on the Schur Decomposition :cite:`Kle:2000` as well a more recently published method by the same author
  (co-authored with Paul Gomme) which computes 2nd-order accurate solutions without using Tensor Algebra :cite:`GomKle:2012`
  (using the Magnus & Neudecker 1999 definition of the Hessian matrix). PyMacLab possesses the added advantage of being equipped with
  an advanced model file parser module, similar to the one available in `Dynare <http://www.dynare.org>`_, which automates cumbersome
  and error-prone (log-)linearization by hand. PyMacLab is also written entirely in Python, is free and incredibly flexible to use and extend.
  
  At this moment, what PyMacLab does not *yet* provide are any methods suitable for pulling in data and estimating deep-structure parameters based on
  some specific estimation framework, such as Bayesian estimation, Maximum Likelihood, Method of Moments or some Limited Information estimation method.
  Having said that, especially LI methods should be easy to "bolt onto" the library as it stands right now by experienced Python programmers using up
  a comparatively little amount of their time. In the near future PyMacLab will provide a few estimation methods for users to work with.

First things first
------------------

  * Documentation at `http://www.pymaclab.com <http://www.pymaclab.com>`_ or `http://packages.python.org/pymaclab/ <http://packages.python.org/pymaclab/>`_
  * Latest development documentation at `http://www.development.pymaclab.com <http://www.development.pymaclab.com>`_
  * Latest source tar ball at `http://pypi.python.org/pypi/pymaclab/ <http://pypi.python.org/pypi/pymaclab/>`_
  * Latest bleeding-edge source via git at `http://github.com/escheffel/pymaclab <http://github.com/escheffel/pymaclab>`_
  * Source code issues tracker at `http://github.com/escheffel/pymaclab/issues/ <http://github.com/escheffel/pymaclab/issues/>`_
  * Download this online documentation as a PDF document :download:`here <PyMacLab.pdf>`.
  
.. note::

    If you want to give PyMacLab a try without installing it onto your own computer you can access an IPython web server frontend in which
    you can create notebooks (as in Mathematica and Maple) from which calls to Python, Numpy, Scipy and PyMacLab can be made. Plots are
    directly rendered to screen. Access this experimental web portal at:
    
                                   `http://www.notebook.pymaclab.com <http://notebook.pymaclab.com>`_
    
    There you will also find provided example scripts which can be run inside your webbrowser.

Features at a Glance
--------------------
  * No "paper-and-pencil" linearization required, done automatically by parsing a DSGE model file.
  * Solutions based on analytical computation of Jacobian and Hessian of models using `Sympycore <http://www.sympy.org/>`_.
  * DSGE models are Python DSGE class instances, treat them as if they were data structures, pass them around, copy them, stack them into arrays,
    and work with many of them simultaneously!
  * Loop over a DSGE model instance thousands of times to alter the parameter space, each time re-computing the solution.
  * Choose from closed form or non-linear steady state solvers or a combination of both.
  * Choose from a number of tried and tested perturbation methods, such as Klein's 1st order accurate and Klein & Gomme's 2nd order accurate methods.
  * Solving models is as fast as using optimized compiled C or Fortran code, expensive computation of analytical Jacobian and Hessian employs parallelized multi-core CPU approach.
  * DSGE example models are provided, including very complex ones such as the one based on Christiano, Eichenbaum and Evans (2001) :cite:`ChrEicEva:2005`.
  * Benefit from a large and growing set of convenience methods to simulate models and plot filtered simulated series as well as impulse-response functions.
  * Carry out advanced empirical macroeconometric analyses using the VAR and FAVAR classes which come provided.
  * Use PyMacLab as a free Python library within a rich and rapidly evolving Python software ecosystem for scientists.
  * Enjoy the power, flexibility and extensibility of the Python programming language and the open-source transparency of PyMacLab.
  * PyMacLab is free as in freedom and distributed under a `Apache 2.0 license <http://www.apache.org/licenses/LICENSE-2.0.html>`_.

.. note::

    PyMacLab is currently known to work well but continues to mature. This documentation site is well under way but still work-in-progress.
    If you have used PyMacLab already and spotted some bugs or felt that some other important features are missing, you can head over to the
    library's `Github <http://github.com/escheffel/pymaclab/>`_ repository to submit an Issue item. We are currently in the process of adding
    more example DSGE model files (and eliminating mistakes in already existing ones). If you have used PyMacLab yourself and want to contribute
    your own DSGE model files we are happy to include them! Finally, to better understand PyMacLab's inner workings, take a look at the API
    documentation.

.. raw:: latex

   \newpage

Documentation
=============

Introduction
------------

:doc:`What is PyMacLab? <pymaclab_intro>`
    This is a succinct introduction to PyMacLab including an explanation of its current features.
:doc:`Philosophy behind PyMacLab <pymaclab_philo>`
    Here I discuss the basic Philosophy behind PyMacLab and what it sets out to do now and in the near future.
:doc:`Why Macroeconomics in Python? <pymaclab_python>`
    In this section I touch upon the the pros and cons of doing Macroeconomics or scientific computing using Python in general.


Series of Brief Tutorials
-------------------------

1) :doc:`Basic DSGE tutorial <tutorial/started_tutorial>`
    Brief tutorial on how to use PyMacLab to work with DSGE models.
2) :doc:`PyMacLab DSGE instance tutorial <tutorial/dsge_instance_tutorial>`
    Succinct tutorial facilitating the understanding of the DSGE OOP data structure in PyMacLab.
3) :doc:`PyMacLab DSGE instance updater tutorial <tutorial/dsge_instance_updater_tutorial>`
    Tutorial on how to use DSGE model instance's intelligent runtime update features.
4) :doc:`PyMacLab DSGE steady state solver tutorial <tutorial/steady_solver_tutorial>`
    This section illustrates various options available to solve DSGE models' steady state.
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

:doc:`linsci_scratch`
    Building a Linux scientific environment from scratch.

:doc:`bibliography`
    Reference list of academic articles and books related to the solution of DSGE models or Python programming.

:doc:`history`
    History of current and past releases

.. raw:: latex

   \newpage

Download & Installation
=======================

Introduction
------------

  PyMacLab is known to work with any of Python version greater than or equal to 2.4 and smaller than 3.0.
  In the future we will consider implementing a compatibility branch for versions of Python greater
  than or equal to 3.0, once all core dependencies are known to have been migrated as well. PyMacLab is always
  extensively tested on Linux and is therefore well supported on this platform. In particular, the author of
  PyMacLab is running his hardware on `Slackware 14.0 <http://www.slackware.com>`_, but other distributions such
  as `Ubuntu <http://www.ubuntu.com>`_ should also work.
  
  PyMacLab will also work on Windows and MacOS so long as users are capable and willing to navigate the
  `murky waters <http://www.scipy.org/Installing_SciPy>`_ of getting a Numpy/Scipy environment set up on their operating
  systems, which because of BLAS and LAPACK dependencies can on occasion be tricky. The internet is littered with explanations
  of how to do this so I will refrain from repeating it here. I should point out however that any Python/Numpy/Scipy system
  definitely requires system-wide available BLAS and LAPACK installations as well as available C++ and Fortran compilers.
  At least one reason for this is that PyMacLab compiles and links in Klein & Gomme's solution routines during installation,
  which they provide as Fortran source files and which come packaged with PyMacLab. Obviously without an installed Fortran
  compiler and a correctly configured system this part of PyMacLab's installation routine is prone to failure.
  
  In Linux these features may come installed by default, in other "user-oriented" operating systems this may not be the case.
  In particular, using Windows, users are best advised to employ the `MinGW32 <http://mingw.org/>`_ Linux system clone and
  to set up a scientific Python environment there. Again, the Numpy/Scipy website contains `help pages <http://scipy.github.com/building/windows.html>`_
  which describe how to do this. Macintosh users are encouraged to take a look at `Scipy Superpack <http://fonnesbeck.github.com/ScipySuperpack/>`_
  which appears to be the better choice over the alternative `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_,
  which is also available for Windows (EPD).
  
  No matter which route users choose to install PyMacLab, the rule of thumb is that so long as they manage to *compile* both
  Numpy and Scipy from their *source files* without problems, installing PyMacLab should also pose no further difficulties. The
  key to success is to have detectable BLAS and LAPACK libraries as well as required compilers installed on the system, where
  *in particular* a good (free) `Fortran compiler <http://gcc.gnu.org/fortran/>`_ will be *absolutely* necessary. In the long run,
  I may consider making pre-built binaries for various platforms available so that users can bypass the error-prone setup using
  compilation from source.

Dependencies
-------------

  Proper functioning of PyMacLab depends on a number of additional Python libraries already being installed on
  your system, such as:

  * `Numpy <http://numpy.scipy.org/>`_
  * `Scipy <http://www.scipy.org/>`_,
  * `Sympycore <http://www.sympy.org>`_,
  * `Parallel Python <http://www.parallelpython.com/>`_
  * `Matplotlib <http://matplotlib.sourceforge.net/>`_
  * `Pandas <http://pandas.pydata.org/>`_

  Sympycore and Parallel Python come distributed with PyMacLab and will be installed along with the main library; the other
  required Python libraries need to be installed separately before an installation of PyMacLab is attempted. All of the
  mentioned scientific packages are great libraries by themselves and should be checked out by any serious scientist interested
  in doing work in Python.
  
  The Pandas data library is *not* needed by the DSGE-modelling features of PyMacLab itself, but is instead required in the experimentally
  made available modules used to estimated and work with VAR and FAVAR models. These modules are in the ``pymaclab.stats.`` branch and
  some test files are included in the test/stats directory.

  If you want to enjoy a Matlab-style interactive environment in which to execute and inspect DSGE and other data structures,
  you'd be hard-pressed to pass over the brilliant and now extra features-ladden `IPython <http://ipython.org/>`_. When downloading
  and installing pymaclab using ``pip`` all of these dependencies should be installed automatically for you, if they are not already
  present on your system. Following right below is a list of options users have to install PyMacLab on their Python-ready computers.

  If you already have a working Python programming environment with some of the above libraries installed, you may want to consider
  installing PyMacLab in its own isolated execution environment using `virtualenv <http://pypi.python.org/pypi/virtualenv>`_ which would
  ensure that your existing system Python installation would remain untouched by PyMacLab's setup routine and its dependency resolution.

Option 1
----------

  You can download the source code of PyMacLab right here. Alternatively, PyMacLab is also hosted at
  `PyPI <http://pypi.python.org/pypi/pymaclab/>`_ and can be installed in the usual way by executing the
  command inside a Linux shell using ``pip``::

    sudo pip install pymaclab
    
  Using this option will also automatically take care of the dependencies by downloading and installing them on-the-fly whenever they
  are not already encountered on the system.

Option 2
---------

  Otherwise get the latest source code compressed as a tarball here:

    `pymaclab.tar.gz <http://pypi.python.org/packages/source/p/pymaclab/pymaclab-0.90.1.tar.gz>`_

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
  Alternatively you can also download a zip file containing the latest "bleeding-edge" version of PyMacLab by
  clicking `here <https://github.com/escheffel/pymaclab/zipball/master>`_.

.. raw:: latex

   \newpage

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

  Last but most certainly not least, my expression of thanks go to my former PhD supervisor `Max Gillman <http://www.maxgillman.com>`_
  who has introduced me to the world of general equilibrium macroeconomics and to monetary macroeconomics more deeply.
  Also, many of the lectures once delivered by `Martin Ellison <http://www.economics.ox.ac.uk/members/martin.ellison/>`_
  formerly at the Economics Department at Warwick now at Oxford made a lasting impression on me.

Online Resources
================

    .. rst-class:: html-plain-table

    ====================== ===================================================
    Author Homepage:       `<http://www.ericscheffel.com>`_
    Github Homepage:       `<http://github.com/escheffel/pymaclab>`_
    Scipy Homepage:        `<http://www.scipy.org>`_
    Download & PyPI:       `<http://pypi.python.org/pypi/pymaclab>`_
    Python Tutorial:       `<http://docs.python.org/tutorial/>`_
    ====================== ===================================================
