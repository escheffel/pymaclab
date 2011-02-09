TODO: change names of sections for steady-state

First the steady-state is solved for. If the .mod file section Steady State 
Non-Linear System [Manual] is provided then dsge.solvers.steadystate.Fsolve
class.  Else if the Steady States [Closed Form] has an entry then these are 
used to solve for the steady-date. Non-linear equations are solved using 
fsolve, which is called for in the .mod file using ROOT().


Next, the model is solved. There are three solution methods available.

Linear Methods
--------------
Rely on log-linearized model equations in .mod file.
Uhlig's method (MatUhlig and PyUhlig) and Klein's Method (MatKlein and 
ForKlein).

1st-Order Non-Linear Methods
----------------------------
Rely on non-linear focs in .mod file
Woodford's Method (MatWood) and Klein's method (ForKleinD)

2nd-Order Non-Linear Methods
----------------------------
Rely on Variance-Covariance Matrix being specified in .mod file.
Klein's method (MatKlein2D and PyKlein2D)

