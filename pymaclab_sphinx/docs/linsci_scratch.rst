

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
  handle these tasks as fast and efficiently as possible, conditional on any given hardware platform researchers may work with.
  
  In this short text I will describe how to set up a Python numerical/scientific environment in Linux which handles the above mentioned tasks -
  as well as others - and focus on the step-by-step installation of various individual packages, which together make up such an environment. Among other,
  this will make clear how many modern interactive scientific environments continue to follow the tried-and-tested approach of often exposing
  their core functions as wrappers to open-source industry standard libraries of which many have been around for decades.
  
  PyMacLab which uses Numpy and Scipy makes no exception in this regard. What this discussion will make clear is that at their core,
  scientific programming suites such as Matlab, Octave, Scilab, Gauss, Numpy/Scipy or indeed PyMacLab (which builds on Numpy/Scipy) are in some sense
  all identical in that they inherit the functionality of tried-and-tested libraries which have been around for decades and which are permitted to be
  included in free as well as proprietary software packages. What sets many of these packages apart are the added methods they provide and the syntactic
  choices they make.

*Starting the Build*

  As just mentioned PyMacLab makes use of the Numpy and Scipy libraries which users often tend to download in tandem. This is because Numpy
  is a library, whose main contribution lies in the provision of Python numerical array and matrix classes, apart from some small number of
  often statistical convenience methods employing these new data types, supplies little else in the domain of advanced scientific computing.
  This is where Scipy then enters which provides many such routines, such as for instance many currently popular optimisation algorithms.
  
  So Scipy requires Numpy to install and run, but does Numpy and its numerical data types also require some other package or library to run?
  As it turns out, if one wishes to run Numpy as efficiently and as fast as possible, some reference BLAS and LAPACK implementations need to be
  installed on the system in form of dynamically loadable libraries. One some systems and for inexperienced Linux users it may be no easy feat
  to configure their system to correctly and most efficiently make use of these tried and tested industry standard libraries.
  
  *All* matrix manipulation software, whether it is for-pay such as Matlab or free such as Octave or Scilab or indeed Numpy/Scipy, *always* relies
  on the BLAS and LAPACK libraries. This is because they are free, have been around for decades and have been tweaked for performance for such a
  long period of time as to make them close to perfect. They are also coded in Fortran, which makes them super-efficient. It would be impossible
  to re-invent the wheel here (although see comments about GotoBlas further below). That is why all of these packages use them.

  BLAS stands for "Basic Linear Algebra Suite" and LAPACK stands for "Linear Algebra Package". BLAS is the starting point for everything and indeed
  LAPACK needs an existing BLAS installation in order to compile and run, so BLAS is a core dependency for LAPACK. Both are usually installed as Fortran
  libraries and are hence very fast and efficient, but BLAS, which contains the most rudimentary matrix and vector functionality can be heavily optimized
  for different hardware set-ups, i.e. different CPU architectures.
  
  Enter the special BLAS library with the name ATLAS = "Automatically tuned linear algebra
  suite", which is essentially another BLAS implementation, but one which detects the hardware setup when it gets compiled from source and intelligently
  adapts its own source code to fit the existing hardware perfectly, like a glove, so-to-speak. The performance differences between a vanilla
  run-of-the-mill BLAS reference library and the optimized ATLAS version can be staggering, depending on what kind of hardware you are running.
  Professional scientists running on expensive number-crunching hardware will always use ATLAS or something similar.

  And since I am describing everything here in full I might as well not conceal the fact that in 2007 or around that time a Japanese computer
  genius with the surname "Goto" hand-coded using assembly language a BLAS implementation for different hardware platforms which beats even ATLAS in speed
  by some margin and is called GotoBlas or GotoBlas2, depending on what kind of version you prefer. The moment the public heard of this, I think he got
  snatched up by Microsoft or something like that. Also, he stopped working on GotoBlas and it kind of got neglected. Some Chinese researchers took
  GotoBlas' source code and continue to maintain and improve the code basis under the name "OpenBlas". You may also find other people which may have started
  maintaining their own version of GotoBlas, but to my knowledge OpenBlas is the most popular one.

  PyMacLab does *not* work correctly with a Numpy version which was compiled against an OpenBlas library, which is a shame, but tolerable given that the speed
  difference between ATLAS BLAS and GotoBlas is not that big. Actually it does work, but it breaks the use of parallel computation using PP
  (Parallel Python) in PyMacLab, which is one of the great features of PyMacLab helping to speed it up considerably exploiting multi-core CPUs. I found out
  about this bug a while ago but I don't know it's exact cause. The point to take away from this is to only use Numpy compiled against either a standard
  vanilla BLAS implementation or the older faster implementation called ATLAS.
  
  And what was the point of all of what I have written so far? Only to tell you that PyMacLab needs Numpy, Numpy needs LAPACK and LAPACK needs BLAS
  (of which there exist various versions, some of which are tightly optimized for existing hardware architectures). When installing PyMacLab, both some
  BLAS and LAPACK need to be installed system-wide in your OS' library directory, which must be the same in a Macintosh OS as in Linux/BSD, which after
  all is just some variant of BSD Unix which itself is structured similarly to Linux.