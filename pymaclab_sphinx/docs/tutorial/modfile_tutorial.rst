.. index:: tutorial; model file; OOP; class

.. raw:: latex

   \newpage


Tutorial 2 - The PyMacLab model file
====================================

The PyMacLab DSGE model file
----------------------------

  In order to be able to load or instantiate your first DSGE model and work with it, you have to make sure to first fill in a so-called PyMacLab
  DSGE model file. The idea behind this is the same as with the Dynare model file which typically ends with a .mod file suffix. PyMacLab already
  comes provided with a number of such files pre-compiled for you to experiment with. For instance the most basic real business cycle model is
  described in the model file ``rbc1_num.txt`` (and in many other similar files which demonstrate the different ways of computing the steady state).
  This file looks as follows:

    ::

      %Model Description+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      This is just a standard RBC model, as you can see.


      %Model Information+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Name = Standard RBC Model, NUM-SS;


      %Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      rho       = 0.36;
      delta     = 0.025;
      betta     = 1.0/1.01;
      eta	= 2.0; 
      psi	= 0.95;
      z_bar     = 1.0;
      sigma_eps = 0.052; 


      %Variable Vectors+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      [1]  k(t):capital{endo}[log,bk]
      [2]  c(t):consumption{con}[log,bk]
      [4]  y(t):output{con}[log,bk]      
      [5]  z(t):eps(t):productivity{exo}[log,bk]
      [6]  @inv(t):investment[log,bk]
      [7]  @R(t):rrate

      %Boundary Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      None


      %Variable Substitution Non-Linear System++++++++++++++++++++++++++++++++++++++
      [1]   @inv(t)   = k(t)-(1-delta)*k(t-1);
      [2]   @inv_bar  = SS{@inv(t)};
      [2]   @F(t)     = z(t)*k(t-1)**rho;
      [2]   @Fk(t)    = DIFF{@F(t),k(t-1)};
      [2]   @Fk_bar   = SS{@Fk(t)};
      [2]   @F_bar    = SS{@F(t)};
      [3]   @R(t)     = 1+DIFF{@F(t),k(t-1)}-delta;
      [4]   @R_bar    = SS{@R(t)};
      [3]   @R(t+1)   = FF_1{@R(t)};
      [4]   @U(t)     = @I{eta!=1.0}{c(t)**(1-eta)/(1-eta)}+@I{eta==1.0}{LOG(c(t))};
      [5]   @MU(t)    = DIFF{@U(t),c(t)};
      [5]   @MU_bar   = SS{@U(t)};
      [6]   @MU(t+1)  = FF_1{@MU(t)};



      %Non-Linear First-Order Conditions++++++++++++++++++++++++++++++++++++++++++++
      # Insert here the non-linear FOCs in format g(x)=0

      [1]   @F(t)-@inv(t)-c(t) = 0;
      [2]   betta*(@MU(t+1)/@MU(t))*@R(t+1)-1 = 0;
      [3]   @F(t)-y(t) = 0;
      [4]   LOG(E(t)|z(t+1))-psi*LOG(z(t)) = 0;


      %Steady States [Closed Form]+++++++++++++++++++++++++++++++++++++++++++++++++++
      None


      %Steady State Non-Linear System [Manual]+++++++++++++++++++++++++++++++++++++++
      [1]   @F_bar-@inv_bar-c_bar = 0;
      [2]   y_bar-@F_bar = 0;
      [3]   betta*@R_bar-1 = 0;
      [4]   betta*R_bar-1 = 0;

      [1]   c_bar = 1.0;
      [2]   y_bar = 1.0;
      [2]   k_bar = 1.0;
      [3]   R_bar = 1.01;

      %Log-Linearized Model Equations++++++++++++++++++++++++++++++++++++++++++++++++
      None


      %Variance-Covariance Matrix++++++++++++++++++++++++++++++++++++++++++++++++++++
      Sigma = [sigma_eps**2];


      %End Of Model File+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  So what does this file mean, and in particular, what is the meaning and purpose of the individual sections?
  These and related questions are addressed in the sections to follow below. They mostly discuss the syntax
  conventions model builders have to adhere to in order to use PyMacLab correctly.

A Description of the model file's individual sections
-----------------------------------------------------

*Model Description Section*

  In the model description section of the model file you can use plain text in order to described more verbosely
  the type of the model summarized in the file, perhaps added with references to important academic journal articles
  in which the model appeared first.

*Information Section*

  This section allows you to add more succinct model properties, including a shorter denominator given by `Name=`
  qualifier. These shorter attributes will then be attached to the model instance where they help to uniquely identify
  the model. In contrast to the information contained in the previous section these qualifiers should be short. You can
  also add yet another item using the `Desc=` qualifier. So an example of this would be:
  
    ::
  
      %Model Information+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Name = Standard RBC Model, NUM-SS;
      Desc = A fairly canonical RBC model with endogenous labour and physical capital;

*Parameters Section*

  As the name suggests, this section provides space for writing down the model's deep and presumably invariable parameters
  which are important as they appear in functionals such as the household's utility or the firm's production function. Don't
  forget to close each declaration with a semi-colon, as this is one of the text parser's conventions. It is also important,
  at least for the time being, to adhere to the convention of employing only floats in this section and to refrain from using
  integers. So this is discouraged:
  
    ::
  
      %Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      rho        = 1/3;
      ss_labour  = 1/3;
    
  And instead one should use:
  
    ::
  
      %Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      rho        = 1.0/3.0;
      ss_labour  = 1.0/3.0;

*Variable Vectors Section*

  This section is very important as it contains a summary of all of the (time-subscripted) variables of the model. The general format
  of this section for each variable is:

   ::

      [1] x(t):var_name{endo|con|exo}[log,hp|bk|cf]

  The first element is a descriptor of how the time-subscripted variable will appear in the system of nonlinear equations. The second
  descriptor is a more revealing but still short name, such as `capital` or `consumption`. It is preferable to write longer variable names
  with an underscore, such as for example `physical_capital` or `human_capital`. Thirdly, the descriptor in curly brackets allows you to
  specifically mark of each variable as either, control variable, endogenous state or exogenous state variable, using optimal control theory
  language. These are inserted in abbreviated style using either `con`, `endo` or `exo`.
  
  Finally, the last option given enclosed in squared brackets allows for two additional options to be specified. Supplying the keyword `log`
  means that the approximation of the model should be formed about the log of the variable, while the last option allows to supply a filtering
  option which is applied to the computation of results based on simulations of the solved model. Currently available choices are either `hp`
  for the HP-Filter, `bk` for the Baxter-King-Filter or `cf` for the Christiano-Fitzgerald filter. Notice that for exogenous variables you also
  have to specify the name of the iid shock:
  
    ::
    
      [7] x(t):eps(t):var_name{endo|con|exo}[log,hp|bk|cf]

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
     [3]   @F_bar    = SS{@F(t)};
     [4]   @R(t)     = 1+DIFF{@F(t),k(t-1)}-delta;
     [5]   @R(t+1)   = FF_1{@R(t)};
     [6]   @U(t)     = @I{eta!=1.0}{c(t)**(1-eta)/(1-eta)}+@I{eta==1.0}{LOG(c(t))};
     [7]   @MU(t)    = DIFF{@U(t),c(t)};
     [8]   @MU(t+1)  = FF_1{@MU(t)};

  These can then be used in the following section instead of having to work with the full expressions instead. Additionally, convience operators
  are accessible, given by:

   ::

      DIFF{EXPRESSION,x(t)}     # replaced by first derivate if expression w.r.t. x(t)

      SS{EXPRESSION}            # expression is converted to its steady state equivalent

      FF_X{EXPRESSION}          # replaced with expression forwarded in time by X periods.
                                # Timing of the information set for expectations
                                # is unchanged!

      BB_X{EXPRESSION}          # replaced with expression lagged in time by X periods.
                                # Timing if the information set for expectations
                                # is unchanged!
                             
      @ALL{EXPRESSION,[0-1],SS} # short-hand way of declaring many items in one go.
                                # See further below for detailed explanation.
      
      @DISCOUNT                 # Special reserved keyword to define the discount factor.
                                # See further below for detailed explanation.
                                
      @I{CONDITION}{EXPRESSION} # An indicator function which inserts the expression whenever
                                # the condition evaluates to TRUE.

  When declaring replacement items in this section make sure to adhere to the syntax of always naming them beginning with a @. Also, within this
  section substitutions within substitutions are permitted. Replacement items for steady-state calculations in the subsequent sections can also
  be supplied here, but have to be of the form such as:

   ::

      [1]   @F_bar   = z_bar*k_bar**rho;

  In PyMacLab steady state expressions of variables strictly have to adhere to the `x_bar` naming convention, i.e. be expressed by the stem
  variable name abbreviation followed by and underscore and the word `bar`. Finally, the DIFF{EXPRESSION,x(t)} is smart enough to differentiate
  across different time periods. So as an example with habit persistence in consumption our utility function depends on current and past consumption:
  
   ::
    
      [1]   @DISCOUNT = betta;
      [2]   @U(t)     = LOG(c(t)-B*c(t-1));
      [3]   @Uc(t)    = DIFF{@U(t),c(t)};
      
  Here the differentiation operator is smart enough to forward the expression by one period before taking the derivative w.r.t to c(t).
  In fact, internally the above will be replaced with:
   
   ::
     
      [1]   @DISCOUNT = betta;
      [2]   @U(t)     = LOG(c(t)-B*c(t-1));
      [3]   @Uc(t)    = DIFF{LOG(c(t)-B*c(t-1))+betta*LOG(E(t)|c(t+1)-B*c(t)),c(t)};
      
  This feature only works if the special reserved keyword @DISCOUNT is defined at the top of the list. This tells PyMacLab which discount rate to
  apply to future (or past) expressions. Finally, as of version 0.95.1 PyMacLab also supports another keyword which works as a short-cut to declare
  a large number of possible derivatives using only one command. This feature would work as follows:
  
    ::
    
      %Variable Substitution Non-Linear System+++++++++++++
      # The utility function and its derivatives
      [1]   @MU(t)     = LOG(c(t))+em(t-1)**(1-1/ups)/(1-1/ups);
      [2]   @ALL{@MU(t),[0-1],SS};
      
  This command takes all of the partial derivatives (but no cross-partials!) of the supplied function `@MU(t)` both for the current and the future period,
  i.e period running from `[0-1]`. One could also specify this as a list like `[0,1]`. If the additional optional argument `SS` is also supplied then
  the steady state versions of both the original function and the derivatives would be declared. Essentially, the above is just a short-hand for the
  following manually declared version:
  
    ::
    
      %Variable Substitution Non-Linear System+++++++++++++
      # The utility function and its derivatives
      [1]   @MU(t)     = LOG(c(t))+em(t-1)**(1-1/ups)/(1-1/ups);
      [2]   @MU_bar    = SS{@MU(t)};
      [3]   @MUc(t)    = DIFF{@MU(t),c(t)};
      [4]   @MUc_bar   = SS{@MUc(t)};
      [5]   @MUem(t)   = DIFF{@MU(t),em(t-1)};
      [6]   @MUem_bar  = SS{@MUem(t)};
      [7]   @MU(t+1)   = FF_1{@MU(t)};
      [8]   @MUc(t+1)  = DIFF{@MU(t+1),E(t)|c(t+1)};
      [9]   @MUem(t+1) = DIFF{@MU(t+1),em(t)};
      
  Obviously, for reasons of brevity using the `@ALL` command is a much better option, in particular if the derivatives and steady state expressions one works
  with are kind of standard and flow naturally from the functional forms of utlity and production functions, for instance.
  
  .. note::

    The whole point of having the subsitutions section present in the library as a functionality to draw on is to reduce systems to a lower dimensionality
    without having to string together algebraic fragments into enormous mathematical expressions which are hard to read and understand by somebody who has
    not been involved in designing the model. This approach also reduces the likelihood of introducing mistakes. With the substitution systems everything
    looks clean and the intuition is immediately discernable from the simplified first-order conditions containing the substitution declarations.

*Non-Linear First-Order Conditions Section*

  In this section users can supply the model's first order conditions of optimality which are passed to PyMacLab for differentiation and
  evaluation. So to use the example from the RBC1 example file given above, filling in this section would look as follows:

   ::

      [1]   @F(t)-@inv(t)-c(t) = 0;
      [2]   betta*(@MU(t+1)/@MU(t))*@R(t+1)-1 = 0;
      [3]   @F(t)-y(t) = 0;
      [4]   LOG(E(t)|z(t+1))-psi*LOG(z(t)) = 0;

  where we have made ample use of the convenient substitution definitions declared in the previous section. Expressions, such as the law of
  motion for the productivity shock, can be supplied in logs for the sake of readability, but otherwise could also alternatively be written as:

   ::

      [4]   E(t)|z(t+1)/(z(t)**psi) = 0;

   .. deprecated:: 0.85 In previous versions of PyMacLab it was possible to write down the law of motion of exogenous states without expectations, i.e.
      `z(t)/(z(t-1)**psi) = 0;`. This behaviour is now deprecated and no longer supported.

*Steady States [Closed Form] Section*

  For relatively simple models, closed form solutions for the steady state may exist and can be entered here as follows:

   ::

      betta   = 1.0/R_bar;
      k_bar   = ((rho*z_bar)/(R_bar - 1 + delta))**(1.0/(1 - rho));
      y_bar   = (z_bar*k_bar)**rho;
      c_bar   = y_bar - delta*k_bar;

  Note that not only steady-state variables like `x_bar` can be supplied here, but indeed any variable who's steady-state value has to be
  determined endogenously withing the model. Sometimes, depending on the model builder's assumptions taken, this could also involve the'
  determination of a parameter such as `betta`.
  
  Sometimes the model's full steady-state can be best determined using a combination of closed form expressions AND the additional numerical
  solution of a system on nonlinear equations, as is the case in the model file provided as ``rbc1_res.txt``. Notice that here one set of steady state
  variables are calculated in closed from, given the knowledge of a set of other steady state variables, while these in turn are first solved
  for in the section using the nonlinear root-finding algorithm. This make sense as for many DSGE models a core set of steady state variables in
  physical capital and marginal utlity - as an example - can be computed using the non-linear root finder, while all of the other variables' steady
  states follow immediately residually from this.
   

*Steady State Non-Linear System [Manual] Section*

  In this section a partial list of or the entire model's variables' steady states can be determined numerically here using good starting values
  and a Newton-like root-finder algorithm. So this section would something like this:

   ::

      %Steady State Non-Linear System [Manual]+++++++++++++
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
  
  As of version 0.95.1 of the PyMacLab library users can also use symbolic expressions in the starting values subsection following the non-linear
  system of equations, such as for instance:
  
    ::
    
      %Steady State Non-Linear System [Manual]+++++++++++++
      [1]   z_bar*k_bar**(rho)-delta*k_bar-c_bar = 0;
      [2]   rho*z_bar*k_bar**(rho-1)+(1-delta)-R_bar = 0;
      [3]   (betta*R_bar)-1 = 0;
      [4]   z_bar*k_bar**(rho)-y_bar = 0;

      [1]  k_bar = 30.0;
      [2]  y_bar = k_bar**alpha;
      [3]  c_bar = 2.0;
      [4]  betta = 1.0;
      
  Finally, again as of version 0.95.1, users can instead declare in this section the following:

    ::
    
      %Steady State Non-Linear System [Manual]+++++++++++++
      USE_FOCS=[0,1,2,3];

      [1]  k_bar = 30.0;
      [2]  y_bar = k_bar**alpha;
      [3]  c_bar = 2.0;
      [4]  betta = 1.0;
      
  When using this `USE_FOCS` command, users are instructing the DSGE model instance to automatically form steady state versions of the non-linear
  system of equations, but doing this only for the equation numbers provided in the passed vector, i.e. `[0,1,2,3]`, which instructs PyMacLab to pick
  equations 1,2,3,4 out of the system of FOCs declared before this section. Python uses 0-indexed vectors, that is why the list starts with
  0 and not 1. If the FOCs are ordered differently, one can also employ different orderings, such as `[0,2,3,4]`. The point here is to have a way of
  disregarding certain equations we may not want to include in the non-linear root finding algorithm, such as certain exogenous laws of motion for which
  we may have calibrated steady state values and do not have to look for them.
  

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
  their own, and not intermingled with model description items. As usual comments are identified by beginning a new line with the hash symbol #.

  Finally, in all sections where it may be applicable, the operators `LOG(x)` and `EXP(x)` can be employed, where the former takes the natural
  logarithm of expression x while the latter raises e to the power x. An example of this would be:

   ::

      [1]   @U(t)   = LOG(c(t));


More than one way to feed in model properties
---------------------------------------------

  As of PyMacLab version 0.95.1, there now exists more than one way to populate a DSGE model instance with information about the properties/features which
  comprise the model and dictate its ultimate behaviour. These changes have been implemented in order to make PyMacLab's feature set more compatible with a
  programming paradigm often called "Meta-programming" or "Template programming" which encapsulates the idea of allowing programs to change their own
  "source code" or otherwise usually assumed fixed features during runtime.
  
  Or at a more basic level, it simply offers a comfortable way for users of the library to change essential features of DSGE models or alternatively swap
  features between them while a program is running. This makes PyMacLab far more powerful in principle than for instance Dynare. So besides reading in a
  conformable DSGE model file from your computer's file system, which other ways are on offer to populate a DSGE model instance?
  
  Instead of passing the model file's full path as a string to the DSGE model at instantiation time, we could have also alternatively passed the actual
  model file itself as a big triple-quoted string to the DSGE class generating instances. This could be defined inside a Python (batch) script and could for
  instance be done like this:
   
   
    .. sourcecode:: python
    
      modstr='''
      %Model Description+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      This is just a standard RBC model, as you can see.


      %Model Information+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Name = Standard RBC Model, RES-SS;


      %Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      rho       = 0.36;
      delta     = 0.025;
      R_bar     = 1.01;
      betta     = 1.0/R_bar;
      eta	= 2.0; 
      psi	= 0.95;
      z_bar     = 1.0;
      sigma_eps = 0.052; 


      %Variable Vectors+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      [1]  k(t):capital{endo}[log,bk]
      [2]  c(t):consumption{con}[log,bk]
      [4]  y(t):output{con}[log,bk]      
      [5]  z(t):eps(t):productivity{exo}[log,bk]
      [6]  @inv(t):investment[log,bk]
      [7]  @R(t):rrate

      %Boundary Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      None


      %Variable Substitution Non-Linear System++++++++++++++++++++++++++++++++++++++
      # Special discount variable
      [1]   @DISCOUNT = betta;
      [1]   @inv(t)   = k(t)-(1-delta)*k(t-1);
      [2]   @inv_bar  = SS{@inv(t)};
      [2]   @F(t)     = z(t)*k(t-1)**rho;
      [2]   @Fk(t)    = DIFF{@F(t),k(t-1)};
      [2]   @Fk_bar   = SS{@Fk(t)};
      [2]   @F_bar    = SS{@F(t)};
      [3]   @R(t)     = 1+DIFF{@F(t),k(t-1)}-delta;
      [4]   @R_bar    = SS{@R(t)};
      [3]   @R(t+1)   = FF_1{@R(t)};
      [4]   @U(t)     = c(t)**(1-eta)/(1-eta);
      [5]   @MU(t)    = DIFF{@U(t),c(t)};
      [5]   @MU_bar   = SS{@U(t)};
      [6]   @MU(t+1)  = FF_1{@MU(t)};



      %Non-Linear First-Order Conditions+++++++++++++++++++++++++++++++++++++++++++++
      # Insert here the non-linear FOCs in format g(x)=0

      [1]   @F(t)-@inv(t)-c(t) = 0;
      [2]   betta*(@MU(t+1)/@MU(t))*@R(t+1)-1 = 0;
      [3]   @F(t)-y(t) = 0;
      [4]   LOG(E(t)|z(t+1))-psi*LOG(z(t)) = 0;


      %Steady States [Closed Form]+++++++++++++++++++++++++++++++++++++++++++++++++++
      [1]   y_bar = @F_bar;


      %Steady State Non-Linear System [Manual]+++++++++++++++++++++++++++++++++++++++
      [1]   @F_bar-@inv_bar-c_bar = 0;
      [2]   betta*@R_bar-1 = 0;
      [3]   betta*R_bar-1 = 0;

      [1]   c_bar = 1.0;
      [2]   k_bar = 1.0;
      [3]   betta = 0.9;

      %Log-Linearized Model Equations++++++++++++++++++++++++++++++++++++++++++++++++
      None


      %Variance-Covariance Matrix++++++++++++++++++++++++++++++++++++++++++++++++++++
      Sigma = [sigma_eps**2];


      %End Of Model File+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      '''
      
      import pymaclab as pm
      
      rbc = pm.newMOD(modstr,mesg=True,ncpus='auto')

  As you can see, the declared Python object ``modstr`` is just a string which holds a standard PyMacLab model file in its entirety (with line breaks!)
  This is then passed to the DSGE class to instantiate a new model and internally PyMacLab recognizes this not as a full path pointer to a physical model
  file existing in your computer's file system but instead as the contents of the file itself ready for direct processing.
  
  Finally, yet one more way open
  to users to instantiate and populate a new DSGE model with its characteristic features is closely related to the one described immediately above. This
  second way uses a Python templating library called ``wheezy.template`` which allows conformable PyMacLab model files to be generated on the fly from
  within a running Python script using a standard Python dictionary of DSGE model properties. Such a dictionary is always created by default and then
  attached to each DSGE model instance whenever they are created and is held inside the object ``model.template_paramdic``. For a simple RBC model this would
  look like:
  
    .. sourcecode:: python
    
      # Load the library and the models branch
      In [1]: import pymaclab as pm
      In [2]: import pymaclab.modfiles.models as models
      In [3]: import pymaclab.modfiles.templates.wheezy_template as template
      
      # Now instantiate the model
      In [3]: rbc = pm.newMOD(models.stable.rbc1_num, mesg=True)
      
      # Check contents of the template dictionary
      In [4]: rbc.template_paramdic.keys()
      
      # These are only the keys of the dictionary, but check the contents yourself
      # to see that they are all standard Python data structures describing the model
      Out[1]:
      ['use_focs',
       'vardic',
       'sigma',
       'mod_desc',
       'subs_list',
       'focs_list',
       'manss_sys',
       'mod_name',
       'llsys_list',
       'paramdic',
       'ssidic',
       'ssys_list']
       
       # Now use the template to automatically generate
       # a conformable PyMacLab model file string
       In [5]: modstr = template.render(rbc.template_paramdic)
       
       # Now print the modstr and check what it looks like
       In [6]: print modstr
       
       Out[2]:
       '''
       %Model Description++++++++++++++++++++++++++++++++++
       None
 
 
       %Model Information++++++++++++++++++++++++++++++++++
       # Short model name
       Name = Standard RBC Model, NUM-SS;
       # Short model description
 
 
 
       %Parameters+++++++++++++++++++++++++++++++++++++++++
       [1]   z_bar = 1.0;
       [2]   psi = 0.95;
       [3]   sigma_eps = 0.052;
       [4]   betta = 0.990099009901;
       [5]   eta = 2.0;
       [6]   rho = 0.36;
       [7]   delta = 0.025;
 
 
 
       %Variable Vectors++++++++++++++++++++++++++++++++++++
       [1]   k(t):capital{endo}[log,bk]
       [1]   c(t):consumption{con}[log,bk]
       [2]   y(t):output{con}[log,bk]
       [1]   z(t):eps(t):productivity{exo}[log,bk]
       [1]   @inv(t):investment [log,bk]
       [2]   @R(t):rrate


       %Boundary Conditions+++++++++++++++++++++++++++++++++
       None
 
 
       %Variable Substitution Non-Linear System+++++++++++++
       [1]   @inv(t) = k(t)-(1-delta)*k(t-1);
       [2]   @inv_bar = SS{@inv(t)};
       [3]   @F(t) = z(t)*k(t-1)**rho;
       [4]   @Fk(t) = DIFF{@F(t),k(t-1)};
       [5]   @Fk_bar = SS{@Fk(t)};
       [6]   @F_bar = SS{@F(t)};
       [7]   @R(t) = 1+DIFF{@F(t),k(t-1)}-delta;
       [8]   @R_bar = SS{@R(t)};
       [9]   @R(t+1) = FF_1{@R(t)};
       [10]   @U(t) = c(t)**(1-eta)/(1-eta);
       [11]   @MU(t) = DIFF{@U(t),c(t)};
       [12]   @MU_bar = SS{@U(t)};
       [13]   @MU(t+1) = FF_1{@MU(t)};
 
 
       %Non-Linear First-Order Conditions+++++++++++++++++++
       [1]   @F(t)-@inv(t)-c(t) = 0;
       [2]   betta*(@MU(t+1)/@MU(t))*@R(t+1)-1 = 0;
       [3]   @F(t)-y(t) = 0;
       [4]   LOG(E(t)|z(t+1))-psi*LOG(z(t)) = 0;
 
 
       %Steady States [Closed Form]++++++++++++++++++++++++++
       None
 
 
       %Steady State Non-Linear System [Manual]+++++++++++++
       [1]   @F_bar-@inv_bar-c_bar = 0;
       [2]   y_bar-@F_bar = 0;
       [3]   betta*@R_bar-1 = 0;
       [4]   betta*R_bar-1 = 0;
 
       [1]   c_bar = 1.0;
       [2]   k_bar = 1.0;
       [3]   y_bar = 1.0;
       [4]   R_bar = 1.01;
 
 
 
       %Log-Linearized Model Equations++++++++++++++++++++++
       None
 
 
       %Variance-Covariance Matrix++++++++++++++++++++++++++
       Sigma = [ 0.002704 ];
       
       
       %End Of Model File+++++++++++++++++++++++++++++++++++
       '''
       
       # You could now check if the model also loads with the generated modfile string
       In [7]: rbc_alt = pm.newMOD(modstr, mesg=True)

       
  As you can see, with the power of templating engines [#f1]_ such as ``wheezy.template`` we can generate PyMacLab-conformable DSGE model files
  on-the-fly by passing simple Python data structures to the template and calling its ``render()`` method. In the above script, the DSGE models called
  ``rbc`` and ``rbc_alt`` will be identical save for small numerical discrepancies introduced because of floating-point arithmatics imprecision.
  In the near future PyMacLab will include another template which will allow the automatic generation of Dynare-conformable model files, allowing users
  to compare and contrast results computed in both environments.
  
.. rubric:: Footnotes

.. [#f1] There exist far more popular templating engines than ``wheezy.template``. One of such, perhaps the most popular, is an engine library called
         ``jinja2`` which is often used by programmers to design dynamic webpages. The other candidate is a library called ``cheetah``. In spite of
         ``wheezy.template``'s lesser popularity, it was chosen for PyMacLab because it claims to be the fastest template engine of all of the
         above mentioned candidates.