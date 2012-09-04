PyMacLab - The Python Macroeconomics Laboratory
===============================================

About
-------
PyMacLab is the Python Macroeconomics Laboratory which currently primarily serves the purpose of providing
a convenience framework written in Python to solve non-linear DSGE models easily. At the time of this writing the library
supports solving DSGE models using 1st and 2nd order perturbation methods which are computed around the steady state.
In particular, the library provides wrapper functions for [Paul Klein's](http://paulklein.ca/newsite/start/start.php)
1st-order accurate method based on the Schur Decomposition as well a more recently published method by the same author
(co-authored with Paul Gomme, see [here](http://ideas.repec.org/a/eee/dyncon/v35y2011i4p604-615.html)) which provides
2nd-order accurate solutions without using Tensor Algebra (using the Magnus and Neudecker 1999 definition of the
Hessian matrix).

The library is extremely user-friendly in the sense of employing a model text file parser similar to that present in
[Dynare](http://www.dynare.org/) which requires users to only write down the original set of non-linear first-order
conditions of optimality. In addition, users are offered a menu of options of how to provide information required for
calculating the steady state of the model. Once the model is parsed and read in, several options of solving it exist
and users are provided with further convenience methods suitable for simulating solved models and investigating dynamic
statistical properties.

It should also be mentioned that because PyMacLab is a convenience library of highly modular nature (using
a object-oriented programming approach) it is very easy to loop over one model several thousand times each time changing
the original set of primitive parameters, such as depreciation rates, impatience factors, etc. in order to compute
solutions of the same model over a large set of conceivable parameterisations. Also, whenever solution methods require
the calculation of the Jacobian or Hessian, this is always done analytically (symbolically) using the Python
symbolic computation library [SympyCore](http://code.google.com/p/sympycore/) and not numerically as in other software
packages. Sympycore is not supplanted by Sympy, but it works well at the moment so we will alter PyMacLab at a later
stage to reflect this.

PyMacLab was authored by [Eric M. Scheffel](http://www.ericscheffel.com) who is currently working as [Assistant Professor
in Economics at Nottingham University China](http://www.nottingham.edu.cn/en/business/people/staffprofile/eric-scheffel.aspx)
and is distributed under the Apache License v2.0. When using this Python library please make sure to cite
the project (using Latex) as:

```latex
@Misc{,
  author =    {Eric M. Scheffel},
  title =     {{PyMacLab}: Open source {Python} Macroeconomics Laboratory},
  year =      {2007--},
  url = "http://github.com/escheffel/pymaclab/"
}
```


Dependencies
-------
[numpy](http://numpy.scipy.org/)  
[scipy](http://www.scipy.org/)  
[sympycore](http://sympy.org/en/index.html) (comes supplied with pymaclab)
[pp](http://www.parallelpython.com/) (optional - but highly recommended for speed)
[ipython](http://ipython.org/) (optional)    
[mlabwrap](http://mlabwrap.sourceforge.net/) (optional)  
[scikits.timeseries](http://pytseries.sourceforge.net/) (optional) 

While numpy, scipy, and sympycore are absolutely necessary to run PyMacLab, Parallel Python (PP)
can autodetect multi-core processors and can potentially speed up the computation of elements such as the analytical
Hessian of a non-linear DSGE model, while mlabwrap have been marked as deprecated and are in the process of
being pruned from the code. Scikits.timeseries is not used in any of the classes themselves, but is used in the example
test file that ships with the library, so it is not essential. To explore the use of PyMacLab interactively a special
IPython environment can be launched with pre-loaded convenience functions, but again, ipython is not a necessary
dependency.

Installation
------------
```python
python setup.py build
sudo python setup.py install
```

Usage
-----
The below code snippet is indicative of how PyMacLab can be used in an environment such as IPython. Since the library
is still undergoing change, it is always advisable to enter the subdirectory containing the test.py test-script and
inspect/run this to better understand PyMacLab's syntax.
```python
# Import the PyMacLab library
import pymaclab as pm
# Import a convenience handler to filepaths for supplied DSGE model files
from pymaclab.modfiles import models

# Model Instantiation with data
rbc1 = pm.newMOD(models.rbc1,db1)

# Model instantiation without data
rbc2 = pm.newMOD(models.rbc2)

# Solve the models using different methods
rbc1.modsolvers.forkleind.solve()
rbc2.modsolvers.pyklein2d.solve()

# Simulate the models after they have been solved
rbc1.modsolvers.forkleind.sim(400)
rbc2.modsolvers.pyklein2d.sim(1000)
```
