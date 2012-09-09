.. index:: cloud; sphinx theme, sphinx theme; cloud

=======================
PyMacLab Tutorial
=======================

Steady State Solution Methods - Building Blocks
===============================================

*Introduction*

  In the previous tutorial we discovered the general structure of the PyMacLab DSGE model instance and saw how this programming approach lent
  itself well to the idea of inspecting and exploring the instantiated models' current state, summarized by its data fields and supplied
  instance methods equipping them with functionality. This section finally discusses how the library equips DSGE model instances with methods
  which make use of the models' computed Jacobian and Hessian, which are evaluated at the models' numerical steady state. As a short reminder,
  we may recall here that it is often this step of obtaining a steady state can prove difficult in the case of a number of well-known models.
  Notwithstanding, for the remainder of this section we will assume that a steady state has successfully been attained and that the model's
  Jacobian and Hessian have been computed. Let's first start with our usual setup of lines of code entered into an IPython shell:

