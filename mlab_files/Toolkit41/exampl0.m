% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL0.M:
% Solving the stochastic neoclassical growth model with the "toolkit"

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

disp('EXAMPLE 0: The stochastic neoclassical growth model');


disp('Hit any key when ready...');
pause;


% Setting parameters:

Z_bar     = 1;    % Normalization
rho       = .36;  % Capital share
delta     = .025; % Depreciation rate for capital
R_bar     = 1.01; % One percent real interest per quarter
eta       = 1.0;  % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
psi       = .95;  % autocorrelation of technology shock
sigma_eps = .712; % Standard deviation of technology shock.  Units: Percent.

% Calculating the steady state:

betta   = 1.0/R_bar;  % Discount factor beta
K_bar   = ((rho*Z_bar)/(R_bar - 1 + delta))^(1.0/(1 - rho));
Y_bar   = Z_bar*K_bar^rho;
C_bar   = Y_bar - delta*K_bar;


% Declaring the matrices. 


VARNAMES = ['capital    ',
            'consumption',
            'return     ',
            'output     ',
            'technology '];

% Translating into coefficient matrices.  
% The loglinearized equations are, conveniently ordered:
% 1) 0 = - (1 - betta(1-delta))(1-rho) k(t-1) + (1 - betta(1-delta)) z(t) - r(t)
% 2) 0 = - K/C k(t) + K/(betta C) k(t-1) + (1 + delta K/C) z(t) - c(t)
% 3) 0 = - y(t) + z(t) + rho * k(t-1)
% 4) 0 = E_t [ - eta c(t+1) + r(t+1) + eta c(t) ]
% 5) z(t+1) = psi z(t) + epsilon(t+1)
% CHECK: 5 equations, 5 variables.
% The third equation is superfluous for the calculation of the
% dynamics, but was added to enable calculations of correlations with
% output: this can only done, if there is an output variable!
%
% Endogenous state variables "x(t)": k(t)
% Endogenous other variables "y(t)": c(t),  r(t), y(t)
% Exogenous state variables  "z(t)": z(t).
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

% DETERMINISTIC EQUATIONS:

% for k(t):
AA = [ 0
       - K_bar/C_bar 
       0             ];

% for k(t-1):
BB = [ - (1 - betta*(1-delta))*(1-rho)
       K_bar/(betta * C_bar) 
       rho                             ];

% For [ c(t), r(t), y(t) ]
CC = [     0,   -1,    0   % Equ. 1)
          -1,    0,    0   % Equ. 2)
           0,    0,   -1 ];% Equ. 3)

% For z(t):
DD = [ (1 - betta*(1-delta))
       (1 + delta * K_bar/C_bar)
       1                        ];

% EXPECTATIONAL EQUATIONS:

% For k(t+1)
FF = [ 0 ];

% For k(t)
GG = [ 0 ];

% For k(t-1)
HH = [ 0 ];

% For [ c(t+1), r(t+1), y(t+1) ]
JJ = [    -eta,     1 ,  0     ];

% For [ c(t), r(t), y(t) ]
KK = [   eta,   0,    0  ];

% For z(t+1)
LL = [ 0 ];

% For z(t)
MM = [ 0 ];

% AUTOREGRESSIVE MATRIX FOR z(t)

NN = [psi];

Sigma = [ sigma_eps^2  ];

% Setting the options:

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);

  
PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 4; % Index of output among the variables selected for HP filter
IMP_SELECT = 1:(m_states+n_endog+k_exog);
   %  a vector containing the indices of the variables to be plotted
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
DO_IMPRESP = 0;
DO_SIMUL   = 0; % Calculates Simulations
DO_MOMENTS = 0; % Calculates Moments

% Starting the calculations:

do_it;
