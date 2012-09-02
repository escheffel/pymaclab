.. index:: cloud; sphinx theme, sphinx theme; cloud

========================
PyMacLab Tutorial Series
========================

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

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace
    In [1]: import pymaclab as pm

    # Get the version and author's name
    In [2]: pm.__version__
    '0.8'

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

  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    In [2]: import os

    # Define the relative path to your modfiles
    In [3]: modpath = "../pymaclab/modfiles/"

    # Instantiate a new DSGE model instance like so
    In [4]: rbc1 = pm.newMOD(os.path.join(modpath,"rbc1.txt"))

    # As an example, check the models computed steady stated
    In [5]: print rbc1.sstate
    {'betta': 0.99009900990099009,
    'c_bar': 2.7560505909330626,
    'k_bar': 38.160700489842398,
    'y_bar': 3.7100681031791227}

Alternatively, you can also test some of the DSGE model files which come supplied with PyMacLab's standard installation. For this to work all you have to do is
to import a provided handler module, ``pymaclab.modfiles.models``, which contains all of the DSGE models' names and their correspoding full file paths:
    
  .. sourcecode:: ipython

    # Import the pymaclab module into its namespace, also import os module
    In [1]: import pymaclab as pm
    # Import the DSGE models' filepath handle
    In [2]: from pymaclab.modfiles import models
    
    #Check all of the available models
    In [3]: dir(models)
    ['__builtins__',
     '__doc__',
     '__file__',
     '__name__',
     '__package__',
     '__path__',
     'cee',
     'max1',
     'max2',
     'mbc1',
     'merz',
     'model2',
     'model3',
     'nk_nocapital',
     'rbc1',
     'rbc2',
     'sims']

    # The DSGE models objects in pymaclab.modfiles.models
    # are just references to full file paths, i.e.

    In [4]: pm.modfiles.models.rbc1
    '/usr/lib/python2.7/site-packages/pymaclab/modfiles/rbc1.txt'

    #Instantiate a new DSGE model instance like so
    In [5]: rbc1 = pm.newMOD(models.rbc1)

    #As an example, check the models computed steady stated
    In [6]: print rbc1.sstate
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
``rbc1.txt`` which looks as follows

  ::

      %Model Description+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      This is just a standard RBC model, as you can see.


      %Model Information+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Name = Standard RBC Model;


      %Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      rho       = 0.36;
      delta     = 0.025;
      R_bar     = 1.01; 
      eta	= 2.0; 
      psi	= 0.95;
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
      [1]   @inv(t)   = k(t)-(1-delta)*k(t-1);
      [2]   @F(t)     = z(t)*k(t-1)**rho;
      [3]   @R(t)     = 1+DIFF{@F(t),k(t-1)}-delta;
      [3]   @R(t+1)   = FF_1{@R(t)};
      [4]   @U(t)     = c(t)**(1-eta)/(1-eta);
      [5]   @MU(t)    = DIFF{@U(t),c(t)};
      [6]   @MU(t+1)  = FF_1{@MU(t)};



      %Non-Linear First-Order Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Insert here the non-linear FOCs in format g(x)=0

      [1]   @F(t)-@inv(t)-c(t) = 0;
      [2]   betta*(@MU(t+1)/@MU(t))*@R(t+1)-1 = 0;
      [3]   @F(t)-y(t) = 0;
      [4]   LOG(E(t)|z(t+1))-psi*LOG(z(t)) = 0;


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

      [1]  c_bar = 1.0;
      [2]  k_bar = 1.0;
      [3]  y_bar = 1.0;
      [4]  betta = 1.0;


      %Log-Linearized Model Equations++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      None


      %Variance-Covariance Matrix++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Sigma = [sigma_eps**2];


      %End Of Model File+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  So what does this file mean, and in particular, what is the meaning and purpose of the individual sections?
  These and related questions are addressed in the sections to follow below. They mostly discuss the syntax
  conventions model builders have to adhere to in order to use PyMacLab correctly.

A Description of the model file's individual sections
=====================================================

*Model Description Section*

  In the model description section of the model file you can use plain text in order to described more verbosely
  the type of the model summarized in the file, perhaps added with references to important academic journal articles
  in which the model first appeared.

*Information Section*

  This section allows you to add more succinct model properties, including a shorter denominator given by `Name=`
  qualifier. These shorter attributes will then be attached to the model instance where they help to uniquely identify
  the model. In contrast to the information contained in the previous section these qualifiers should be short.

*Parameters Section*

  As the name suggests, this section provides space for writing down the model's deep and presumably invariable parameters
  which are important as they appear in functionals such as the household's utility or the firm's production function. Don't
  forget to close each declaration with a semi-colon, as this is one of the text parser's conventions.

*Variable Vectors Section*

  This section is very important as it contains a summary of all of the (time-subscripted) variables of the model. The general format
  of this section for each variable is:

   ::

    [1] x(t):var_name{endo|con|exo}[log,hp|bk]

  The first element is a descriptor of how the time-subscripted variable will appear in the system of nonlinear equations. The second
  descriptor is a more revealing but still short name, such as `capital` or `consumption`. It is preferable to write longer variable names
  with an underscore, such as for example `physical_capital` or `human_capital`. Thirdly, the descriptor in curly brackets allows you to
  specifically mark of each variable as either, control variable, endogenous state or exogenous state variable, using optimal control theory
  language. These are inserted in abbreviated style using either `con`, `endo` or `exo`. Finally, the last option given enclosed in squared
  brackets allows for two additional options to be specified. Supplying the keyword `log` means that the approximation of the model showed be
  formed about the log of the variable, while the last option allows to supply a filtering option which is applied to the computation of results
  based on simulations of the solved model. Currently available choices are either `hp` for the HP-Filter or `bk` for the Baxter-King-Filter.

*Boundary Conditions Section*

  This section is currently not in use but has been included for future compatibility with solution methods which are not based on the perturbation
  paradigm.

*Variable Substitution Non-Linear System*
  This is perhaps one of the most useful and convenient sections of the model file. In the section right after this one users are asked to insert
  the DSGE model's firs-order conditions of optimality which can often be quite tedious and long algebraically. One way of giving users a more
  convenient and intuitive way of writing down the model's FOCs is to work with a subsitution system which can be declared in this section. So for
  example if one wanted to write down the expression for output or the Euler equation for physical capital, one could resort to the following
  useful replacement definitions:

   ::

      [1]   @inv(t)   = k(t)-(1-delta)*k(t-1);
      [2]   @F(t)     = z(t)*k(t-1)**rho;
      [3]   @R(t)     = 1+DIFF{@F(t),k(t-1)}-delta;
      [3]   @R(t+1)   = FF_1{@R(t)};
      [4]   @U(t)     = c(t)**(1-eta)/(1-eta);
      [5]   @MU(t)    = DIFF{@U(t),c(t)};
      [6]   @MU(t+1)  = FF_1{@MU(t)};

  These can then be used in the following section instead of having to work with the full expressions instead. Additionally, convience operators
  are accessible, given by:

   ::

      DIFF{EXPRESSION,x(t)}  # is replaced by first derivate if expression w.r.t. x(t)

      FF_X{EXPRESSION}       # is replaced with expression forwarded in time by X periods.
                             # The timing of the information set for expectations is unchanged!

      BB_X{EXPRESSION}       # is replaced with expression lagged in time by X periods.
                             # The timing if the information set for expectations is unchanged!

  When declaring replacement items in this section make sure to adhere to the syntax of always naming them beginning with a @. Also, within this
  section substitutions within substitutions are permitted. Replacement items for steady-state calculations in the subsequent sections can also
  be supplied here, but have to be of the form such as:

   ::

      [1]   @F_bar   = z_bar*k_bar**rho;

  In PyMacLab steady state expressions of variables strictly have to adhere to the `x_bar` naming convention, i.e. be expressed by the stem variable
  name abbreviation followed by and underscore and the word `bar`.

*Non-Linear First-Order Conditions Section*

  In this section users can supply the model's first order conditions of optimality which are passed to PyMacLab for differentiation and evaluation.
  So to use the example from the RBC1 example file given above, filling in this section would look as follows:

   ::

      [1]   @F(t)-@inv(t)-c(t) = 0;
      [2]   betta*(@MU(t+1)/@MU(t))*@R(t+1)-1 = 0;
      [3]   @F(t)-y(t) = 0;
      [4]   LOG(E(t)|z(t+1))-psi*LOG(z(t)) = 0;

  where we have made ample use of the convenient substitution definitions declared in the previous section. Expressions, such as the law of
  motion for the productivity shock, can be supplied in logs for the sake of readability, but otherwise could also alternatively be written as:

   ::

      [4]   E(t)|z(t+1)/(z(t)**psi) = 0;

  Also, for the exogenous state variable such as, again, the productivity shock, it does not matter whether we write down the law of motion as
  in the previous example or alternatively as:

   ::

      [4]   z(t)/(z(t-1)**psi) = 0;

*Steady States [Closed Form] Section*
  For relatively simple models, closed form solutions for the steady state may exist and can be entered here as follows:

   ::

      betta   = 1.0/R_bar;
      k_bar   = ((rho*z_bar)/(R_bar - 1 + delta))**(1.0/(1 - rho));
      y_bar   = (z_bar*k_bar)**rho;
      c_bar   = y_bar - delta*k_bar;

  Note that not only steady-state variables like `x_bar` can be supplied here, but indeed any variable who's steady-state value has to be
  determined endogenously withing the model. Sometimes, depending on the model builder's assumptions taken, this could also involve the'
  determination of a parameter such as `betta`. Sometimes the model's full steady-state can best be determined as a combination of closed form
  expressions AND the additional numerical solution of a system on nonlinear equations. 
   

*Steady State Non-Linear System [Manual] Section*

  In this section a partial list of or the entire model's variables' steady states can be determined numerically here using good starting values
  and a Newton-like root-finder algorithm. So this section would something like this:

   ::

      [1]   z_bar*k_bar**(rho)-delta*k_bar-c_bar = 0;
      [2]   rho*z_bar*k_bar**(rho-1)+(1-delta)-R_bar = 0;
      [3]   (betta*R_bar)-1 = 0;
      [4]   z_bar*k_bar**(rho)-y_bar = 0;

      [1]  c_bar = 1.0;
      [2]  k_bar = 1.0;
      [3]  y_bar = 1.0;
      [4]  betta = 1.0;

  Very often, this section is simply a restatement of the first order conditions of optimality but with time subscripts removed and instead
  replaced with the steady state `x_bar` notation. This section and the previous can often be the most difficult ones to specify well, as many
  more complex DSGE models' steady states are not easy to determine and often require some good judegement, experience and good starting values
  for the root-finding algorithm.

*Log-Linearized Model Equations Section*

  In this section you could theoretically also supply the first-order log-linearized equations manually, such as was necessary in Harald Uhlig's
  toolbox. But this feature is perhaps best relegated to compatibility tests and proof-of-concept experiments to show that PyMacLab's computed
  solutions based on automatic differentiation are identical with the ones computed from this section. An example would be:

   ::

      # foc consumption
      [1]   (1/C_bar)**Theta*X_bar**(Psi*(1-Theta))*x(t)...
           -(1/C_bar)**Theta*X_bar**(Psi*(1-Theta))*c(t)=...
             LAM_bar*lam(t)+A_bar*MU_bar*mu(t);
      # foc leisure
      [2]   (1-Theta)*c(t)+(Psi*(1-Theta)-1)*x(t)=lam(t)+...
             z(t)+(1-alpha)*k(t-1)-(1-alpha)*l(t);

  In this case all variables already have to be interpreted as percentage deviations from steady state. Both in this and in the nonlinear FOCs
  section, model equations DO NOT necessarily have to be expressed as `g(x)=0`, but can also be written as `f(x)=g(x)`. In this case the PyMacLab
  parser simply internally generates `f(x)-g(x) = 0` and works with this expression instead.

*Variance-Covariance Matrix Section*

  The standard way of supplying information on the variance-convariance structure of the iid shocks hitting the laws of motions of the exogenous
  state variables. So this section would look something like this:

   ::

           Sigma = [sigma_eps**2];

  or for more elaborate models like this:

   ::

     Sigma = [sigma_eps**2   0;
              0    sigma_xi**2];

*All sections*

  If in any of the lines of one of the sections the keyword `None` is inserted, even in a section which has otherwise been declared in the correct
  way as described above, then the entire section will be ignored and treated as empty, such as for instance:

   ::

       %Log-Linearized Model Equations++++++++++++++++++++++++++++
       None

  If alebraic expression become to long, one can also employ a line-breaking syntax using the elipsis, such as:

   ::


      [1]   (1-Theta)*c(t)+(Psi*(1-Theta)-1)*x(t)=lam(t)+...
             z(t)+(1-alpha)*k(t-1)-(1-alpha)*l(t);

  Finally, as is customary from other programming languages, comments can also be inserted into DSGE model files. However, in contrast to other
  languages conventions, such as Python itself, at the moment the library will only parse model files correctly if the comments are on a line of
  their own, and not intermingled with model description items. As usual comments are identified by beginning with the hash symbol #.
