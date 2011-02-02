Pymaclab history
----------------
Written by E.M. Scheffler
Takes place of Uhlig's toolkit plus some other helper stuff

Uhlig's Toolkit
---------------
Specify three kinds of variables
    Engogenous State, x
    Endogenous Other, y
    Exogenous State, z

Uses Method of Undetermined Coefficients to solve for the matrices
    P,Q,R, and S such that
    x_t = Px_{t-1]} + Qz_{t}
    y_{t} = Rx_{t-1} + Sz_{t}
