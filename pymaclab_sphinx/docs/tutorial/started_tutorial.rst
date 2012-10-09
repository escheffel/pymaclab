.. raw:: latex

   \newpage

Tutorial 1 - Getting started
============================

Introduction
------------

  PyMacLab's strength or orginal design goal has been that of providing users with a rich and flexible DSGE data structure (also called
  class in object-oriented programming speak) which allows them to do lots of interesting things with DSGE models. Don't forget that another
  stated goal of PyMacLab is to permit researchers to work with DSGE models who know how to and enjoy to use the programming language Python.
  This tutorial assumes that you already know how to use Python, if you are completely new to it I suggest you consult one of the
  `many tutorials <http://docs.python.org/tutorial/>`_ available on the internet to get started with the language. Don't forget that Python is
  free as in free as freedom, no proprietary software such as Matlab is tainting the freedom of the environment you will be working in when using PyMacLab.
  
  The easiest and best way to get started and explore all of the available features of PyMacLab, or for that matter any other Python library, is to
  work with `IPython <http://ipython.org/>`_, which is an abbreviation for Interactive Python. With this program, Python users who are scientifically
  inclined, can more or less emulate the behaviour of programs such as Matlab, Octave or Scilab. This means that sessions can be launched in which commands
  entered into a command prompt are executed immediately and resulting output is printed to screen. It also means that generated objects resulting from such
  computations carried out at the command prompt (such as matrices or DSGE models) begin to exist within the session's variable scope and become ready for
  further manipulation, i.e. they persist. Before being able to work with PyMacLab we have to launch an IPython session and import the library into it's
  scope:

    .. sourcecode:: ipython

      # Import the pymaclab module into its namespace
      In [1]: import pymaclab as pm

      # Get the version and author's name
      In [2]: pm.__version__
      '0.95.1'

      # Get the library's author's name
      In [3]: pm.__author__
      'Eric M. Scheffel'

  Here we simply have imported the PyMacLab module and inspected some of its attributes, such as the current version numbering as well as them
  module's author's name. Let's look deeper into the recesses of the module though to better understand who it is organized

    .. sourcecode:: ipython

      # Import the pymaclab module into its namespace
      In [1]: import pymaclab as pm

      # Use the dir() command to view all available attributes and method calls,
      # this command returns a list
      In [2]: dir(pm)
      ['OPS',
      '__builtins__',
      '__doc__',
      '__file__',
      '__name__',
      '__package__',
      '__path__',
      'dattrans',
      'db_graph',
      'dsge',
      'explain',
      'favar',
      'filters',
      'ldbs',
      'linalg',
      'lmods',
      'lvars',
      'macrolab',
      'make_modfile',
      'modedit',
      'modfiles',
      'modinfo',
      'modsolve',
      'newDB',
      'newFAVAR',
      'newMOD',
      'newVAR',
      'pyed',
      'stats',
      'sys',
      'texedit',
      'var']


  As you can see the module contains a quite a few attributes, many of which are still experimental and perhaps best not called at this stage. The most mature
  and arguable most interesting method call is that called ``pm.newMOD``, which allows users to instantiate a DSGE model instance, which would be done like so:

    .. sourcecode:: ipython

      # Import the pymaclab module into its namespace, also import os module
      In [1]: import pymaclab as pm
      In [2]: import os

      # Define the relative path to your modfiles
      In [3]: modpath = "../pymaclab/modfiles/models/stable"

      # Instantiate a new DSGE model instance like so
      In [4]: rbc1 = pm.newMOD(os.path.join(modpath,"rbc1_res.txt"))

      # As an example, check the models computed steady stated
      In [5]: print rbc1.sstate
      {'betta': 0.99009900990099009,
      'c_bar': 2.7560505909330626,
      'k_bar': 38.160700489842398,
      'y_bar': 3.7100681031791227}

  Alternatively, you can also test some of the DSGE model files which come supplied with PyMacLab's standard installation. For this to work all you have to do is
  to import a provided handler module, ``pymaclab.modfiles.models``, which contains all of the DSGE models' names and their correspoding full file paths.
  Notice however that the models themselves are further classified into three categories, ``models.stable``, ``models.testing`` and ``models.development``
  which helps to distinguish between models which are in the process of being added and such which are known to work correctly:
    
    .. sourcecode:: ipython

      # Import the pymaclab module into its namespace, also import os module
      In [1]: import pymaclab as pm
      # Import the DSGE models' filepath handle
      In [2]: from pymaclab.modfiles import models
      
      #Check all of the available models in the stable branch
      In [3]: dir(models.stable)
      ['__builtins__',
      '__doc__',
      '__file__',
      '__name__',
      '__package__',
      '__path__',
      'jermann98',
      'jermann98_ext',
      'merz',
      'prog',
      'rbc1_cf',
      'rbc1_ext',
      'rbc1_extss',
      'rbc1_focs',
      'rbc1_num',
      'rbc1_res',
      'rbc1_sug',
      'rbc2',
      'sims']
      
      #Check all of the available models in the development branch
      In [4]: dir(models.development)
     ['RBC_Romer',
      '__builtins__',
      '__doc__',
      '__file__',
      '__name__',
      '__package__',
      '__path__',
      'max1',
      'max2',
      'mbc1',
      'model2',
      'model3',
      'nk_nocapital',
      'nkm',
      'nkm_nocapital']
      
      #Check all of the available models in the testing branch
      In [5]: dir(models.testing)
      ['__builtins__',
      '__doc__',
      '__file__',
      '__name__',
      '__package__',
      '__path__',
      'cee']


      # The DSGE models objects in pymaclab.modfiles.models
      # are just references to full file paths, i.e.

      In [6]: models.stable.rbc1_res
      '/usr/lib/python2.7/site-packages/pymaclab/modfiles/rbc1_res.txt'

      #Instantiate a new DSGE model instance like so
      In [7]: rbc1 = pm.newMOD(models.stable.rbc1_res)

      #As an example, check the models computed steady stated
      In [8]: print rbc1.sstate
      {'betta': 0.99009900990099009,
      'c_bar': 2.7560505909330626,
      'k_bar': 38.160700489842398,
      'y_bar': 3.7100681031791227}


  Now we have already seen some of the power and simplicity we can leverage by using PyMacLab. Before learning some of its additional power, we do however
  still need to take a quick detour to study the model file ``rbc1.txt`` which we had to pass as an argument to the ``pm.newMOD`` method call, as its
  structure is determined by a set of conventions which are important to adhere to in order to enable PyMacLab to parse and employ the information contained
  therein correctly and efficiently.
