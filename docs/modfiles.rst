Sections
========

Model Information
------------------
name=
desc=

Parameters
----------

Variable Vectors
----------------
[#] @vari:viid:varn{vtype}[mod]


# -- the variable number, e.g., 1
@ -- the variable substition marker; optional
vari -- the variable, e.g., z(t)
viid -- the variable id, e.g., eps_z(t); optional (used for exo)
varn -- the variable name, e.g., techshock
vtype -- the type of variable, can be exo, endo, con, or other; optional
            exo --
            endo --
            con --
            other --
mod -- variable modification, can be exp, log, or hp in order to say that
            the variable should be exponentiated, logged, or filtered using
            the Hodrick-Prescott filter. e.g, log,hp; optional

Example Exogenous
[1]  z(t):eps_z(t):techshock{exo}[log,hp]

Boundary Conditions
-------------------

Variable Substitution Non-Linear System
---------------------------------------

Non-Linear First-Order Conditions
---------------------------------

Steady State Non-Linear System [Manual]
---------------------------------------
If you want to input multiline equations, you need to use 
... at the beginning of the following line

Steady States [Closed Form]
---------------------------

Log-Linearized Model Equations
------------------------------

Variance-Covariance Matrix
--------------------------

Minford Model Evaluation
------------------------

