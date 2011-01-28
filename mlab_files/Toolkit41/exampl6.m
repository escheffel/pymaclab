% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL6.M: A small open economy stochastic neoclassical growth model with 
% adjustment costs for capital.

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

disp('EXAMPLE 6: The small open economy stochastic neoclassical growth model');
disp('           New capital is a CES function of old capital and investment');
disp('           (this introduces ``adjustment costs'').');

disp('Hit any key when ready...');
pause;

% The model is: max E[ sum beta^t (C(t)^(1-eta) - 1)/(1 - eta) ]
% s.t. 
% (1) C(t) + X(t) + A(t) = Y(t) + R(t) A(t-1) (feasibility)
% (2) Y(t) = Z(t) K(t-1)^rho                  (production function)
% (3) K(t) = F(t) - delta * K(t-1)            (production of new capital)
% where the intermediate capital F(t) is given by             
% (4) F(t) = ( K(t-1)^theta + X(t)^theta )^(1/theta)  (``fabricated'', undepreciated new capital)  
% Exogenous processes are: Z(t) and R(t) (the "world return" on assets).
% One finds the first order conditions
% (5) 1 = beta E_t [ (C(t)/C(t+1))^eta R(t+1) ]  ( investing in assets )
% (6) 1 = beta E_t[ (C(t)/C(t+1))^eta (F(t)/X(t))^(1-theta) ...
%           ( rho Y(t+1)/K(t) + (X(t+1)/K(t))^(1-theta) - delta (X(t+1)/F(t+1))^(1-theta) )      ]

% Setting parameters:

Z_bar     = 1;    % Normalization
NPV_frac  = .5;   % Steady state asset holdings = NPV_frac * (NPV of output). 
                  % NOTE: NPV_frac = 0 does not lead to a sensible solution, since we are calculating
                  % PERCENTAGE deviations from steady state, and percentage dev. from zero are all infinity!
rho       = .36;  % Capital share
delta     = .025; % Depreciation rate for capital
R_bar     = 1.01; % One percent real interest per quarter.  FIXED DUE TO "SMALL OPEN ECONOMY"
eta       = 1.0;  % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
theta     = .8;  % CES coefficient in F(t) = ( K(t-1)^theta + X(t)^theta )^(1/theta).
                  % We need 0 < theta <= 1.  For theta = 1, one obtains the usual linear case.
psi_z     = .95;  % autocorrelation of technology shock
sigma_z   = .712; % Standard deviation of technology shock.  Units: Percent.
psi_r     = .95;  % autocorrelation of interest rate shock
sigma_r   = 1.0;  % Standard deviation of world interest rate shock.  Units: Percent
corr_z_r  = 0;    % Correlation of technology shock and world interest rate shock.


% Calculating the steady state:

betta   = 1.0/R_bar;  % Discount factor beta LINED UP WITH WORLD RETURN TO GET A STEADY STATE
XK_bar  = ((1+delta)^theta - 1)^(1.0/theta); % investment-to-capital ratio in steady state
FK_omt  = (1+delta)^(1-theta);     % = (F_bar/K_bar)^(1-theta).  Useful constant. "omt" = "one minus theta"
FX_omt  = FK_omt/XK_bar^(1-theta); % = (F_bar/K_bar)^(1-theta).  Useful constant. "omt" = "one minus theta"
YK_bar  = (R_bar - FK_omt + delta)/(rho*FX_omt);  % output-capital ratio.
          %   Note: we need 0 < theta or, better yet, theta close to unity.
K_bar   = (Z_bar/YK_bar)^(1.0/(1-rho));
Y_bar   = Z_bar*K_bar^rho;
X_bar   = XK_bar*K_bar;
F_bar   = (1+delta)*K_bar;
A_bar   = NPV_frac * (Y_bar/(R_bar - 1));
C_bar   = Y_bar - X_bar + (R_bar - 1)*A_bar;

% Declaring the matrices. 


VARNAMES = ['capital    ',
            'assets     ',
            'consumption',
            'investment ',
            'output     ',
            'interm.cap.',
            'technology ',
            'return     ' ];


% Translating into coefficient matrices. 
% The loglinearized equations are, conveniently ordered:
% 1) 0 = Y y(t) + R A (r(t) + a(t-1)) - A a(t) - C c(t) - X x(t)
% 2) 0 = - y(t) + z(t) + rho k(t-1)
% 3) 0 = -  k(t) + (1+delta) f(t) - delta k(t-1)
% 4) 0 = - F_bar^theta f(t) + K_bar^theta k(t-1) + X_bar^theta x(t)
% 5) 0 = E_t[ eta (c(t) - c(t+1)) + r(t+1) ]
% 6) 0 = E_t[ eta (c(t) - c(t+1)) + (1-theta) (f(t)-x(t)) + ...
%             (1/R)*(F/X)^(1-theta)* { rho Y/K (y(t+1)-k(t))  ...
%                                      + (1-theta)(X/K)^(1-theta) (x(t+1)-k(t)) ...
%                                      - delta (1-theta) (X/F)^(1-theta) (x(t+1)-f(t+1) } ]
% 7) z(t) = psi_z z(t-1) + epsilon_z(t)
% 8) r(t) = psi_r r(t-1) + epsilon_r(t)
% CHECK: 8 equations, 8 variables.
% Endogenous state variables "x(t)": k(t), a(t)
% Endogenous other variables "y(t)": c(t), x(t), y(t), f(t)
% Exogenous state variables  "z(t)": z(t), r(t)
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

% DETERMINISTIC EQUATIONS:
% 1) 0 = Y y(t) + R A (r(t) + a(t-1)) - A a(t) - C c(t) - X x(t)
% 2) 0 = - y(t) + z(t) + rho k(t-1)
% 3) 0 = -  k(t) + (1+delta) f(t) - delta k(t-1)
% 4) 0 = - F_bar^theta f(t) + K_bar^theta k(t-1) + X_bar^theta x(t)

% for k(t) and a(t)
AA = [  0, -A_bar
        0, 0 
       -1, 0    
        0, 0  ];

% for      k(t-1) and a(t-1)
BB = [           0, R_bar*A_bar
               rho, 0
            -delta, 0       
       K_bar^theta, 0           ];

% For [  c(t),       x(t), y(t),          f(t) ]
CC = [ -C_bar,     -X_bar,Y_bar,             0   % Equ. 1)
            0,          0,   -1,             0   % Equ. 2)
            0,          0,    0,     (1+delta)   % Equ. 3)
            0,X_bar^theta,    0,(-F_bar^theta) ];% Equ. 4)

% For z(t) and r(t)
DD = [ 0, R_bar*A_bar
       1, 0
       0, 0          
       0, 0          ];

% EXPECTATIONAL EQUATIONS:
% 5) 0 = E_t[ eta (c(t) - c(t+1)) + r(t+1) ]
% 6) 0 = E_t[ eta (c(t) - c(t+1)) + (1-theta) (f(t)-x(t)) + ...
%             (1/R)*(F/X)^(1-theta)* { rho Y/K (y(t+1)-k(t))  ...
%                                      + (1-theta)(X/K)^(1-theta) (x(t+1)-k(t)) ...
%                                      - delta (1-theta) (X/F)^(1-theta) (x(t+1)-f(t+1) } ]

% For k(t+1) and a(t+1)
FF = [    0,    0 
          0,    0 ];

% For k(t) and a(t)
GG = [                                                         0, 0
       ( -((rho/R_bar)*YK_bar*FX_omt + (1-theta)*FK_omt/R_bar) ), 0 ];

% For k(t-1) and a(t-1)
HH = [ 0, 0
       0, 0 ];

% For [ c(t+1),                        x(t+1),                   y(t+1),               f(t+1)]
JJ = [    -eta,                             0,                        0,                    0     
          -eta,(1-theta)*(FK_omt-delta)/R_bar,(rho/R_bar)*YK_bar*FX_omt,delta*(1-theta)/R_bar];

% For [ c(t),   x(t), y(t),    f(t) ]
KK = [   eta,      0,    0,       0
         eta,theta-1,    0, 1-theta ];

% For z(t+1) and r(t+1)
LL = [ 0, 1
       0, 0 ];

% For z(t) and r(t)
MM = [ 0, 0
       0, 0 ];

% AUTOREGRESSIVE MATRIX FOR z(t)

NN = [psi_z,  0
          0,  psi_r ];

Sigma = [ sigma_z^2,                corr_z_r*sigma_z*sigma_r 
          corr_z_r*sigma_z*sigma_r, sigma_r^2                 ];

% Setting the options:
  
[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);


PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 5; % Index of output among the variables selected for HP filter
IMP_SELECT = [1,2,3,5,7,8];
   %  a vector containing the indices of the variables to be plotted
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.


% Starting the calculations:

do_it;
