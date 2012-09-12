

=======================
PyMacLab Tutorial
=======================

Getting started - Basics
========================

*Introduction*

PyMacLab's strength or orginal design goal has been that of providing users with a rich and flexible DSGE data structure (also called
class in object-oriented programming speak) which allows them to do lots of interesting things with DSGE models. Don't forget that another
stated goal of PyMacLab is to permit researchers to work with DSGE models who know how to and enjoy to use the programming language Python.
This tutorial assumes that you already know how to use Python, if you are completely new to it I suggest you consult one of the many tutorials
available on the internet to get started with the language. Don't forget that Python is free as in free as freedom, no proprietary software such
as Matlab is tainting the freedom of the environment you will be working in when using PyMacLab. The easiest way to get started and explore the
features of PyMacLab is to launch a IPython session and to import PyMacLab into it

  ::

    # import the pymaclab module into its namespace
    import pymaclab as pm

    # get the version and author's name
    pm.__version__
    '0.8'
    pm.__author__
    'Eric M. Scheffel'

Here we simply have imported the PyMacLab module and inspected some of its attributes, such as the current version numbering as well as them
module's author's name. Let's look deeper into the recesses of the module though to better understand who it is organized

  ::

    # import the pymaclab module into its namespace
    import pymaclab as pm

    # use the dir() command to view all available attributes and method calls,
    # this command returns a list
    dir(pm)
    ['OPS',
    '__author__',
    '__builtins__',
    '__date__',
    '__doc__',
    '__file__',
    '__loader__',
    '__name__',
    '__package__',
    '__path__',
    '__revision__',
    '__version__',
    'db_graph',
    'dsge',
    'explain',
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
    'newMOD',
    'newVAR',
    'pyed',
    'stats',
    'texedit']

As you can see the module contains a quite a few attributes, many of which are still expermimental and perhaps best not called at this stage. The most mature
and arguable most interesting method call is that called ``pm.newMOD``, which allows users to instantiate a DSGE model instance, which would be done like so:

  ::

    #import the pymaclab module into its namespace, also import os module
    import pymaclab as pm
    import os

    #Define the relative path to your modfiles
    modpath = "../pymaclab/modfiles/"

    #Instantiate a new DSGE model instance like so
    rbc1 = pm.newMOD(os.path.join(modpath,"rbc1.txt"))

    #As an example, check the models computed steady stated
    print rbc1.sstate

    {'betta': 0.99009900990099009,
    'c_bar': 2.7560505909330626,
    'k_bar': 38.160700489842398,
    'y_bar': 3.7100681031791227}


Now we have already seen some of the power and simplicity we can leverage by using PyMacLab. Before learning some of its additional power, we do however
still need to take a quick detour to study the model file ``rbc1.txt`` which we had to pass as an argument to the ``pm.newMOD`` method call, as its
structure is determined by a set of conventions which are important to adhere to in order to enable PyMacLab to parse and employ the information contained
therein correctly and efficiently.


The PyMacLab DSGE model file
============================
  In order to be able to load or instantiate your first DSGE model and work with it, you have to make sure to first fill in a so-called PyMacLab
  DSGE model file. The idea behing this is the same as the Dynare model file which typically ends in .mod. PyMacLab already comes provided with a
  number of such files pre-compiled for you to experiment with. For instance the most basic real business cycle model is described in the model file
  ``rbc1.txt`` which looks as follows::

    %Model Description+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    This is just a standard RBC model, as you can see.
    
    
    %Model Information+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Name = Standard RBC Model;
    
    
    %Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    rho       = 0.36;
    delta     = 0.025;
    R_bar     = 1.01; 
    eta       = 2.0; 
    psi       = 0.95;
    z_bar     = 1.0;
    sigma_eps = 0.052; 
    
    
    %Variable Vectors+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    [1]  k(t):capital{endo}[log,bk]
    [2]  c(t):consumption{con}[log,bk]
    [4]  y(t):output{con}[log,bk]      
    [5]  z(t):eps(t):productivity{exo}[log,bk]
    [6]  @inv(t):investment[log,bk]
    [7]  @R(t):rrate
    
    %Boundary Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    None
    
    
    %Variable Substitution Non-Linear System++++++++++++++++++++++++++++++++++++++++++++++++
    [1]   @inv(t) = k(t)-(1-delta)*k(t-1);
    [2]   @R(t) = rho*z(t)**psi*k(t)**(rho-1)+(1-delta);
    [3]   @y(t) = z(t)*k(t-1)**(rho);



    %Non-Linear First-Order Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Insert here the non-linear FOCs in format g(x)=0

    [1]   @y(t)-@inv(t)-c(t) = 0;
    [2]   betta*E(t)|c(t+1)**(-eta)*c(t)**(eta)*@R(t)-1 = 0;
    [3]   z(t)*k(t-1)**(rho)-y(t) = 0;
    [4]   psi*LOG(z(t))-LOG(E(t)|z(t+1)) = 0;


    %Steady States [Closed Form]++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    betta   = 1.0/R_bar;
    k_bar   = ((rho*z_bar)/(R_bar - 1 + delta))**(1.0/(1 - rho));
    y_bar   = (z_bar*k_bar)**rho;
    c_bar   = y_bar - delta*k_bar;

    %Steady State Non-Linear System [Manual]+++++++++++++++++++++++++++++++++++++++++++++++++
    [1]   z_bar*k_bar**(rho)-delta*k_bar-c_bar = 0;
    [2]   rho*z_bar*k_bar**(rho-1)+(1-delta)-R_bar = 0;
    [3]   (betta*R_bar)-1 = 0;
    [4]   z_bar*k_bar**(rho)-y_bar = 0;

    c_bar = 1.0;
    k_bar = 1.0;
    y_bar = 1.0;
    betta = 1.0;


    %Log-Linearized Model Equations++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    None


    %Variance-Covariance Matrix++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Sigma = [sigma_eps**2];


    %End Of Model File+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  So what does this file mean, and in particular, what it the meaning and purpose of the individual sections?

.. testcode::
    import pymaclab as pm
    print pm.__version__
   

.. testoutput::
    :hide:
    :options: -ELLIPSIS, +NORMALIZE_WHITESPACE


*Model Description and Information Section*
  This is this

*Parameters Section*
  This is this

*Variable Vectors Section*
  This is this

*Boundary Conditions Section*
  This is this

*Variable Substitution Non-Linear System*
  This is this

*Non-Linear First-Order Conditions Section*
  This is this

*Steady States [Closed Form] Section*
  This is this

*Steady State Non-Linear System [Manual] Section*
  This is this

*Log-Linearized Model Equations Section*
  This is this

*Variance-Covariance Matrix Section*
  This is this
