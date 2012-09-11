.. index:: cloud; sphinx theme, sphinx theme; cloud

===================
PyMacLab Philosophy
===================

  The main motivation or reason for why I decided to write PyMacLab was a relatively simple one. Nowadays in Macroeconomic research,
  PhD students and researchers often want to solve DSGE models many thousand times over, compare models with each other, and experiment
  with new methods which make use of existing ones. I therefore wanted to afford myself the convenience of working in an environment in which
  DSGE models could be treated and worked with in terms of simple but flexible and powerful OOP data structures. In addition, open-source OOP
  coding in a robust and easy to maintain language with all of the advantages discussed would translate into easy extensibility and quick
  adaptability to whatever needs may arise.

  The whole point of PyMacLab is to view the model template files solely as instantiating information,
  with flexible model "moulding" turning into the standard mode of interaction thereafter. To be able to loop over a DSGE model and extend its
  functionality in the future easily due to its modular structure were important motivating aspects driving this project forward. This
  distinguishes PyMacLab from existing `programs` such as Dynare which feel somewhat as inflexible take-as-is routines with a standardized
  output. PyMacLab, in contrast, is `not` a program, but instead a Python library providing an abstract DSGE model class, from which users can
  instantiate as many instances as they like which in turn exhibit an incredible amount of post-instantiation scope for transformability. In
  Dynare the model solving task after the program has run is often considered finished, in PyMacLab the fun starts after models have been loaded
  from template model files when researchers change models' properties dynamically at runtime.
