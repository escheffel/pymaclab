.. index:: cloud; sphinx theme, sphinx theme; cloud

================================
Macroeconomic Analysis in Python
================================

  PyMacLab, so far, is the only Python library designed with the specific purpose in mind to permit solving DSGE models conveniently. For some
  this may raise the question of why one would want to make available such a library for the Python programming environment in the first place,
  especially in light of the fact that many other alternatives already exist. In this section I will attempt to briefly outline reasons for
  why having and using PyMacLab is of benefit to many potential users and why it fits well into the existing software ecosystem. This section
  may very well read like your usual list of advantages vs. disadvantages. So let's perhaps first start with the advantages of using PyMacLab
  in Python.

Advantages of PyMacLab in Python
================================

*Python's growing scientific user community*

  Perhaps the largest benefit of having access to a library for solving DSGE models programmed in the Python language, is that for various
  easily identifiable reasons, Python is rapidly turning itself into the language best supplied with ready-to-use libraries aimed at the
  requirements of a sophisticated scientific community. An important stepping stone in this larger development was the availability of the
  Numpy/Scipy library suite which has rapidly turned Python into an open-source replacement for proprietary software environments such as
  Matlab.

  But the growing availability of mature scientific libraries has not faltered since and has continued to grow at a dramatic pace. To provide
  an exhaustive list of all production-ready scientific libraries for Python would be a difficult task to achieve, so I will limit myself to a
  few well-known and preferred packages:

  +------------------------------------+----------------------------------------------------------------------------------------------------+
  | Library                            |                                  Description                                                       |
  +====================================+====================================================================================================+
  |Numpy                               | Doing linear algebra and providing an array and matrix data type in Python                         |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |Scipy                               | Library of many scientific routines,such as basic statistics, optimization, filtering, etc.        |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |Statsmodels                         | More advanced and dedicated library for advanced statistics in Python                              |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |Pandas                              | Library providing a data frame and time series data type and a large number of data methods        |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |MDP Toolkit                         | A data processing library with wrappers for unsupervised learning routines, etc.                   |
  +------------------------------------+----------------------------------------------------------------------------------------------------+
  |Matplotlib                          | The Python de facto standard library for all-purpose graphing and plotting                         |
  +------------------------------------+----------------------------------------------------------------------------------------------------+ 
  |Scikits.TimeSeries                  | The first library to provide a convenient library containing an advanced time series data type     |
  +------------------------------------+----------------------------------------------------------------------------------------------------+ 

  This list barely touches the surface of what is currently out there available for free for Python prgrammers wanting to do scientific
  computing. Even if a specific library does not exist directly, it is usually easy to produce wrappers for traditional and mature libraries
  originally written in C, C++ or Fortran. This last point brings me straight to the next advantage users can expect to benefit from when choosing
  to use Python in their scientific work.

*Python glues well into traditional scientific languages*

  Do you come from a Fortran, C or C++ background. Do you still have some old tried and tested routines in source code lying around? Then you
  will be glad to hear that Python has a number of outstanding tools and built-in properties available which allow you to easily link your
  existing source code into Python programs allowing them to be called inside Python scripts as if they were normal Python routines. For C++
  aficionados, there is `Swig <http://swig.org/>`_, a translation tool turning C++ code into Python modules. Fortran users can make use of the
  easy-to-use `Py2f <http://www.scipy.org/F2py>`_ while C users should have no problems whatsoever using the
  `Ctypes <http://docs.python.org/library/ctypes.html>`_ library.

  The latter users are perhaps best supported, as the original implementation
  of Python is actually implemented and written in C itself, which explains why it still does well in terms of execution speed in spite of being
  a dynamically typed and interpreted languages (more on this later). The lesson to take away from this - and this is something Python as a
  language is well-known for - is that Python is a great "gluing" language which allows you to work well with a large number of software
  libraries originally coming from quite disparate software environments/ecosystems. Using tools such as
  `Scikits.mlabwrap <http://mlabwrap.sourceforge.net/>`_ or `RPy2 <http://rpy.sourceforge.net/rpy2.html>`_ you can even interface Python with
  Mathworks Matlab or Gnu R.

*Python is dynamic and interpreted*

  Python by itself is a programming language like any other such as Java, C++ or C and supports pretty much any functionality these languages
  are also capable of. What sets it apart is that it is not compiled and linked, but insteady interpreted and thus belongs to the family of
  scripting languages. This means that compile-time and run-time are woven together in one single environment which makes the programming
  experience much more seemless, interactive and transparent.

  This turns Python into a so-called RAD tool - a rapid application development tool, which dramatically cuts down development time and allows
  developers to design code which is much easier to read and maintain. Also, Python is a mixed language supporting both OOP and procedural code.
  Coupled with the Python-specific interactive shell `IPython <http://ipython.org/>`_ Python programming is just as interactive and dynamic as
  working in a Matlab interactive environment, only much more powerful and flexible and its abilities stretch far beyond matrix algebra and
  scientific computing.. Work for a while in an IPython shell and you will know the difference. Although it has not happened yet, it stands to
  reason to expect that one day a fully-fledged Python compiler may appear, giving developers the choice to compile their programs all the way
  down to machine code.

*PyMacLab in Python encourages learning and extending*

  Many routines aimed at solving DSGE models often feel like canned algorithms which by their very design encourage use of them as simple
  and unreflective input-blackbox-output procedures in which the users are mostly concerned with learning the syntactic rules of the program
  to quickly "get out of it what they need". I feel that for the sake of productivity this is not an entirely wrong or indeed deplorable
  circumstance, quite to the contrary. But it often does imply that users substitute away from learning and understanding under-the-hood
  details of implementation which in themselves would be worthwhile try to come to grasp with as a means of learning. PyMacLab, as a result
  of the language it is written in and the way it chooses to implement DSGE modelling in form of an "intelligent" DSGE model instance,
  encourages students and researchers to look underneath the hood and to use the structure of the DSGE data type as it exists at any given
  point in time in order to develop easy-to-add extensions. The open-source, improved readability and maintainability nature of Python and
  PyMacLab itself further enforce this advantage. Canned routines encourage unreflective use, but does human capital theory not teach us that
  learning-by-doing is an important aspect of stimulating economic growth?



Disadvantages of PyMacLab in Python
===================================

*Python is dynamic and interpreted*

  The previous stated `advantage` of Python is simultaneously also its disadvantage. In many areas of scientific research in which
  heavy-duty `number-crunching` and `brute-force` methods prevail, execution speed is usually perceived as a top priority. Python's dynamism
  comes at the cost of much slower execution speed than comparable source code written in Fortran or C++ compiled all the way down to machine
  code. However, this last point needs to be qualified in light of what has already been pointed out above. Since Python glues in well with
  existing traditional programming languages, it is comparatively easy to design Python programs in which CPU-intensive code is simply
  "outsourced" to a dynamically linked library originally written and compiled in Fortran.

  This last remark is particularly relevant when reference is made to the well-known 20/80 rule of computing, stating that for most computer
  programs 20% of its code uses up 80% of its total execution time. Writing the other 80% of your code in easily maintainable Python source code
  and the remaining 20% in Fortran or another compiled language is a golden recipe which is advocated and applied by many professional users.
  Actually, the execution speed vs. development speed is the only real drawback worth the trouble to mention. And given the above recipe and
  the plausible possibility of one day seeing a real Python compiler, the benefits of Python in scientific computing by far outweigh its
  drawbacks.