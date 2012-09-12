.. index:: tutorial; DSGE instance; dynamic updating

.. raw:: latex

   \newpage

Tutorial 3 - The Python DSGE instance updater methods
=====================================================

Introduction
------------

  In the previous tutorial we described the general structure and some of the behaviour of PyMacLab's DSGE model's object-oriented design and
  what kind of advantages model builders can derive from this. In particular, at the end of the last tutorial, we saw how the decision to design
  DSGE model's instances essentially in from of an object-oriented advanced data structure with well-defined data fields and methods allowed us
  to "loop over" a DSGE model instance many times over, each time feeding the instance with slightly a different value for the impatience factor
  :math:`\beta` and obtaining a new value for the corresponding steady state value of physical capital.

  To be able to do something like this `simply` and `intuitively` is one of the many advantages that PyMacLab has over traditional DSGE model
  solving `programs` such as Dynare. While Dynare is a complete `program` allowing users to solve for specific DSGE models with little added
  post-solution programmatic flexibility, PyMacLab is from ground up designed to be a library providing an object-oriented Class-defined advanced
  data structure with data fields and behaviour in form of methods. This makes it incredibly easy to build programs which treat DSGE models as if
  they were simple data structures, such as integers, floats, lists or any other data structure you are familiar with from using a number of
  different programming languages. This makes PyMacLab incredibly flexible for within- and post-solution `program interventions`, where in this
  section we will define more clearly what we meant by `program interventions`. The main point to take away from all of this is that Dynare is a
  canned program to get specific solutions, while PyMacLab is a Python plugin library, whose only purpose is to make available an advanced DSGE
  model class suitable for carrying out a variety of tasks.

  While in programs like Dynare the most important aspect from the point of view of model builders to obtaining the solution of DSGE models lies
  in the specification of so-called model files, the proper way to understand the more flexible operation of PyMacLab is to view the model
  template only as a specific `point of departure` from which to start your analysis from. The only point of making use of model template files
  is to initialize or instantiate a DSGE model instance, the real power in using PyMacLab lies in the methods made available to researchers which
  become available `after` model files have been read in and the DSGE model instances become available inside the Python interactive environment.
  It is this post-instantiation scope for activity and alteration which makes using PyMacLab so much fun for researchers. As already illustrated
  in the previous tutorial, this design decision opens up the possibility of easily and intuitively obtaining a large permutation of possible
  solution outcomes quickly. There are two convenience methods or avenues open to the researcher to intelligently update DSGE model instances
  dynamically at runtime which we will describe in some detail next.


One-off alteration of one specific model property
-------------------------------------------------

  Once a DSGE model is instantiated using a model template files as a point of departure and a specific source of information from which model
  properties get parsed and attached as data fields to the DSGE model instance, all we have done is to initialize or read into memory a specific
  `state` of our DSGE model. You may also recall that the nature of this `state` differs crucially depending upon what kind of `initlev` parameter
  value was passed at DSGE model instantiation. So we recall that at `initlev=0` only information gets parsed in but nothing at all gets solved,
  not even the steady state. Other values for `initlev` gave rise to instantiated model instances which were solved all down to ever deeper levels
  including possibly including the steady state and dynamic elasticities (policy functions). In this section I will briefly describe methods
  allowing a functionality akin to comparative statics analysis in PyMacLab in which models can get easily be re-initialized with different model
  information injected into existing DSGE model instances dynamically at runtime.

  One way of using PyMacLab to carry out a kind of comparative statics analysis has already been touched upon for one specific case in the
  previous tutorial, using the example of the ``rbc1.updaters.paramdic`` attribute. In general you will discover that navigating into the
  ``rbc1.updaters.`` branch of the model exposes a number of data fields which in exact name are also available at the root ``rbc1.`` of the
  model. So how then does ``rbc1.updaters.paramdic`` differ from its counterpart available at the root ``rbc1.paramdic``?

  The answer to this question is that while ``rbc1.paramdic`` is just the parameter dictionary which get populated at model instantiation time,
  ``rbc1.updaters.paramdic`` is a "wrapped" version of the same dictionary with the added behaviour that if a new values gets inserted into it
  using the relevant dictionary methods, such as ``rbc1.updaters.paramdic['R_bar']=1.02`` or ``rbc1.updates.paramdic.update(dic_alt)`` where
  `dic_alt` is some other dictionary containing updated values of one or more keyed item in the original paramdic, then an internal function
  gets triggered at assignment time which re-initializes the model using the updated values for paramdic.

  If the newly assigned value is exactly identical to the value which was already stored in the original paramdic, then the model will not get
  updated as its state has remained unaltered. The following smart one-off updaters are available all of which possess the above described
  behaviour. Notice also that the DSGE model will only get updated down to the level first specified in `initlev` at model instantiation time.
  Also, when any of these wrapped updater objects just gets called in the normal way they behave exactly as their non-wrapped counterparts by
  returning the values stored in them without any additional behaviour updating the model.

  +------------------------------------+----------------------------------------------------------------------------------------------------+
  | Wrapped updater object             |                                  Description                                                       |
  +====================================+====================================================================================================+
  |``model.updaters.paramdic``         | Inserting updated values re-initializes the model with new set of parameter values                 |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters.vardic``           | Inserting new values updates/changes values and attributes defined for variables used in the model |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters.nlsubsdic``        | Inserting new values updates the dictionary of variable substitution items beginning with @        |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters.foceqs``           | Inserting new string equations updates the set of non-linear first-order conditions of optimality  |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters.manss_sys``        | In this list users can insert new equations as strings into the closed form steady state system    |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters.ssys_list``        | In this list users can insert new equations as strings into the numerical steady state system      |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters.sigma``            | A wrapped matrix which updates the model when changed values are inserted in the varcovar-matrix   |
  +------------------------------------+----------------------------------------------------------------------------------------------------+

  Notice how in this branch ``rbc1.updaters.`` changing wrapped objects by assigning new values will trigger automatic updating of the DSGE
  model immediately upon assignment. This behaviour may not always be desirable whenever a series of changes need to be made before updating of
  the model can be considered. Whenever such situations occur an alternative route needs to be taken which we will explore next.


Altering many model properties and queued processing
------------------------------------------------------

  At times researchers may want to load or instantiate a particular DSGE model instance using a corresponding template file but then perhaps
  plan to radically modify the model dynamically at runtime, by combining such actions as introducing new time-subscripted variables, altering
  the deep parameter space and adding new or augmenting existing equations in the system of non-linear FOCs. Whenever such radical alterations
  are considered, they will often have to happen in combindation `before` the model gets updated using the new information passed to it. In this
  case users will use the same wrapped objects already described above but instead use them in the ``rbc1.updaters_queued.`` branch.

  Here, first a number of changes can be made to objects such as ``rbc1.updaters_queued.paramdic`` or ``rbc1.updaters_queued.foceqs``, etc.
  which by themselves will `not` trigger an automatic model updating functionality. Instead all changes will be put into a queue which will
  then have to be processed manually by calling the method ``rbc1.updaters_queued.process_queue()`` after all desired changes have been made.
  This addscenormous flexibility to model builders' options, as they can essentially build a completely new model at runtime dynamically
  starting from a simple model instantiated at the outset of their Python scripts/batch files. Therefore, this functionality allows users to
  dynamically update all information at runtime which was first parsed from the model template file, each time re-computing the DSGE model's
  new state given the changes made after the call to the queue processing method has been made.

  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  | Wrapped updater object                |                                  Description                                                       |
  +=======================================+====================================================================================================+
  |``model.updaters_queued.paramdic``     | Inserting updated values re-initializes the model with new set of parameter values                 |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.vardic``       | Inserting new values updates/changes values and attributes defined for variables used in the model |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.nlsubsdic``    | Inserting new values updates the dictionary of variable substitution items beginning with @        |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.foceqs``       | Inserting new string equations updates the set of non-linear first-order conditions of optimality  |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.manss_sys``    | In this list users can insert new equations as strings into the closed form steady state system    |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.ssys_list``    | In this list users can insert new equations as strings into the numerical steady state system      |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.sigma``        | A wrapped matrix which updates the model when changed values are inserted in the varcovar-matrix   |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.queue``        | The actual queue. Here objects which have been altered will be stored as strings                   |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
  |``model.updaters_queued.process_queue``| The queue processing method which finally updates the queued objects in the right order            |
  +---------------------------------------+----------------------------------------------------------------------------------------------------+
