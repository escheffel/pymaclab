

=======================================
Building a Linux Scientific Environment
=======================================

Understanding the wider system
==============================

*Introduction*

  Many researchers wanting to do work within some scientific environment on their computers typically depend on software which handles at least
  two kinds of tasks well. One of them is setting up and solving problems in linear algebra and/or problems requiring non-linear or other
  optimisation techniques while the other one is declaring and manipulating algebraic expressions, usually done with so-called CASs, or computer
  algebra systems. Also, sometimes it is not only important to be able to have access to environments which can handle these tasks, but also to
  handle these tasks as fast and efficiently as possible, conditional on any given hardware platform researchers may work with. In this short
  text I will describe how to set up a Python numerical/scientific environment in Linux which handles the above mentioned tasks - as well as
  others - and focus on the step-by-step installation of various individual packages, which together make up such an environment. Among other,
  this will make clear how many modern interactive scientific environments continue to follows the tried-and-tested approach of often exposing
  their core functions as wrappers to open-source industry standard libraries of which many have been around for decades. PyMacLab which uses
  Numpy and Scipy forms no exception here.

*Starting the Build*

  As just mentioned PyMacLab makes use of the Numpy and Scipy libraries which users often tend to download in tandem. This is because Numpy
  is a library, whose main contribution lies in the provision of Python numerical array and matrix classes, apart from some small number of
  often statistical convenience methods employing these new data types, supplies little else in the domain of advanced scientific computing.
  This is where Scipy then enters which provides many such routines, such as for instance many currently popular optimisation algorithms.
  So Scipy requires Numpy to install and run, but does Numpy and its numerical data types also require some other package or library to run?
  As it turns out, if one wishes to run Numpy as efficiently and as fast as possible, some reference BLAS and LAPACK implementations need to be
  installed on the system in form of dynamically loadable libraries. One some systems and for inexperienced Linux users it may be no easy feat
  to configure their system to correctly and most efficiently make use of these tried and tested industry standard libraries.