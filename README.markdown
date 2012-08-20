PyMacLab README
================================

About
-------
PyMacLab stands for Python Macroeconomics Library which currently primarily serves the purposes of providing
a convenience framework written in Python to solve non-linear DSGE models. It was authored by Eric M. Scheffel
and is distributed under the GNU General Public License v3.0. When using this Python library please make sure to
cite the the project (using Latex) as:

@Misc{,
  author =    {Eric M. Scheffel},
  title =     {{PyMacLab}: Open source {Python} Macroeconomics Laboratory},
  year =      {2007--},
  url = "http://github.com/escheffel/pymaclab/"
}


Dependencies
-------
numpy
scipy
sympy
sympycore
ipython (optional)
pp (optional)
mlabwrap (optional)
scikits.timeseries

Installation
------------
```
python setup.py build
sudo python setup.py install
```

Usage
-----
```
import pymaclab as pm

# Model Instantiation with data
rbc1 = pm.newMOD('modfiles/rbc1.txt',db1)

# Model instantiation without data
rbc2 = pm.newMOD('modfiles/rbc2.txt')

# Solve the models using different methods
rbc1.modsolvers.forkleind.solve()
rbc2.modsolvers.pyklein2d.solve()

# Simulate the models after they have been solved
rbc1.modsolvers.forkleind.sim(400)
rbc2.modsolvers.pyklein2d.sim(1000)
```
