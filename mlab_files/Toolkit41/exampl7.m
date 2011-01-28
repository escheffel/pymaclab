% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL7.M computes through a two-country stochastic neoclassical growth model.
% A Cobb-Douglas production function is used for each country, i.e.
% Y(t,i) = Z(t,i) K(t-1,i)^rho(i)  i = 1,2
% Labor is assumed to be immobile.
% where rho(i) is the share parameter for country i, where 
% log Z(t,i) = (1-psi(i)) log Z_bar(i) + psi(i) log Z(t-1,i) + epsilon(t,i), i = 1,2
% The evolution of capital follows
% K(t,i) = (1- delta(i)) K(t-1,i) + X(t,i), i = 1,2
% Resource constraint is:
% C(t,1) + C(t,2) + X(t,1) + X(t,2)  = Y_total(t) = Y(t,1) + Y(t,2) 
% The utility function is supposed to be of the constant intertemp. elast. 
% of substitution type,
% but we allow the elasticities to differ across countries.
% Given some weight 0 < alpha < 1, the social planner maximizes
% alpha U1 + (1-alpha) U2
% subject to the constraints above, where
% Ui = E[ sum beta^t ( c(t,i)^(1-eta_i) - 1)/(1-eta_i) ]
% But rather than using alpha as a parameter, the fraction of country 1 consumed in steady state is
% used as a parameter, and the appropriate alpha is calculated.

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

disp('EXAMPLE 7: The stochastic neoclassical growth model with two sectors');
disp('           The productivities of the sectors are chosen so that ');
disp('           both sectors operate in steady state.  An equal distribution');
disp('           of labor is assumed at the steady state.');

disp('Hit any key when ready...');
pause;


% Setting parameters:

Z_bar_1   = 1;    % Normalization
Z_bar_2   = 1;    % Note, that different countries may operate technologies at different productivities,
                  % since labor is not tradeable.
rho_1     = .36;  % Capital share
rho_2     = .36;  % Capital share
delta_1   = .025; % Depreciation rate for capital
delta_2   = .025; % Depreciation rate for capital
N_bar_1   = 1/2;  % Assumption about steady state distribution of labor
N_bar_2   = 1/2;  % Assumption about steady state distribution of labor
C_1_frac  = .5;   % Fraction of entire consumption to be consumed in country 1 in steady state.
betta     = 0.99; % One percent discounting per quarter
eta_1     = 1.0;  % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
eta_2     = 1.0;  % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)

psi_1     = .5;  % autocorrelation of technology shock
psi_2     = .5;  % autocorrelation of technology shock

sigma_eps = .712; % Standard deviation of technology shock.  Units: Percent.
corr_12   = 0;    % correlation between techn. shocks 1 and 2.

   
% Calculating the steady state:

R_bar_1   = 1.0/betta;  % Steady state return.  Note: no variance correction.
R_bar_2   = 1.0/betta;  % Steady state return.  Note: no variance correction.

YK_bar_1  = (R_bar_1 - 1 + delta_1)/rho_1;
YK_bar_2  = (R_bar_2 - 1 + delta_2)/rho_2;

% This follows from the Cobb-Douglas production function:

K_bar_1 = N_bar_1 * ( Z_bar_1 / YK_bar_1 )^(1.0/(1-rho_1));
K_bar_2 = N_bar_2 * ( Z_bar_2 / YK_bar_2 )^(1.0/(1-rho_2));

Y_bar_1 = YK_bar_1 * K_bar_1;
Y_bar_2 = YK_bar_2 * K_bar_2;

X_bar_1 = delta_1 * K_bar_1;
X_bar_2 = delta_2 * K_bar_2;


Y_bar = Y_bar_1 + Y_bar_2;
X_bar = X_bar_1 + X_bar_2;

C_bar   = Y_bar - X_bar;
C_bar_1 = C_1_frac * C_bar;
C_bar_2 = (1-C_1_frac)*C_bar;

% Note: In the social planners problem, it must be the case that
% alpha C(t,1)^(-eta_1) = (1 - alpha) C(t,2)^(-eta_2).  Thus, in steady state,
% C_bar_1^(-eta_1)/C_bar_2^(-eta2) = (1/alpha - 1)
% This is exploited now:
alpha  = 1.0/(C_bar_1^(-eta_1)/C_bar_2^(-eta_2)  + 1);



% Declaring the matrices. 


VARNAMES = ['capital 1   ',
            'capital 2   ',
            'total cons. ',
            'total output',
            'tot.investm.',
            'consumpt. 1 ',
            'consumpt. 2 ',
            'output 1    ',
            'output 2    ',
            'investment 1',
            'investment 2',
            'return 1    ',
            'return 2    ',
            'technology 1',
            'technology 2',];

% Translating into coefficient matrices.  
% The loglinearized equations are, conveniently ordered:
% 1) 0 = C c(t) + X x(t) - Y y(t)
% 2) 0 = - X x(t) + X1 x(t,1) + X2 x(t,2) 
% 3) 0 = - Y y(t) + Y1 y(t,1) + Y2 y(t,2)
% 4) 0 = - C c(t) + C1 c(t,1) + C2 c(t,2)
% 5) 0 = eta_1 c(t,1) - eta_2 c(t,2)
%   (Note: here, marginal utility is equated period by period)
% 6,7)  0 = - y(t,i) + z(t,i) + rho_i k(t-1,i) 
% 8,9)  0 = - k(t,i) + (1-delta_i) k(t-1,i) + delta_i x(t,i)
% 10,11) 0 = - R_i r(t,i) + rho_i YK_bar_i (y(t,i) -  k(t-1,i) )
% 12,13) 0 = E_t [ - eta_i c(t+1,i) + r(t+1,i) + eta_i c(t,i) ]
%   (Note: in principle, there are four Euler conditions.  But since marginal utilities
%          are equated across countries, see (5), two suffice, provided, we use both r(t+1,i)  )
% 14,15) z(t+1,i) = psi_i z(t,i) + epsilon(t+1,i)
% CHECK: 15 equations, 15 variables.

% Endogenous state variables "x(t)": k(t,i), i = 1,2
% Endogenous other variables "y(t)": 
%                 c(t), y(t), x(t), n(t,i), y(t,i), x(t,i), r(t,i), i = 1,2
% Exogenous state variables  "z(t)": z(t,i), i = 1,2.
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

% for k(t,i):
AA = [  0,  0,  
        0,  0,
        0,  0,
        0,  0,
        0,  0, 
        0,  0,
        0,  0,
       -1,  0,
        0, -1,
        0,  0
        0,  0 ];

% for k(t-1,i):
BB = [  0,  0, 
        0,  0, 
        0,  0,  
        0,  0,  
        0,  0, 
    rho_1,  0, 
        0,rho_2,
  (1-delta_1),0,
      0,(1-delta_2),
 -rho_1*YK_bar_1,0,
     0,-rho_2*YK_bar_2 ];

% For   c(t),  y(t),  x(t),   c(t,i),   y(t,i),   x(t,i),    r(t,i)
CC = [ C_bar,-Y_bar, X_bar,    0,  0,    0,  0,    0,  0,    0,  0        % 1)
           0, 0,-X_bar,        0,  0,    0,0, X_bar_1,X_bar_2, 0,0        % 2)
           0,-Y_bar, 0,        0,0, Y_bar_1,Y_bar_2, 0,0,    0,  0        % 3)
      -C_bar,0,0,       C_bar_1,C_bar_2,   0,0,    0,  0,    0,  0        % 4)
           0,0,0,       eta_1,-eta_2,    0,  0,    0,  0,    0,  0        % 5)
           0,0,0,              0,  0,   -1,  0,    0,  0,    0,  0        % 6)
           0,0,0,              0,  0,    0, -1,    0,  0,    0,  0        % 7)
           0,0,0,              0,  0,    0,  0,  delta_1,0,  0,  0        % 8)
           0,0,0,              0,  0,    0,  0,  0,delta_2,  0,  0        % 9)
           0,0,0,              0,0,  rho_1*YK_bar_1,0,  0,0, -R_bar_1,0   % 10)
           0,0,0,              0,0,  0,rho_2*YK_bar_2,  0,0, 0,-R_bar_2 ];% 11)

% For z(t,i):
DD = [ 0, 0, 
       0, 0, 
       0, 0, 
       0, 0, 
       0, 0, 
       1, 0, 
       0, 1, 
       0, 0
       0, 0
       0, 0
       0, 0 ];


% EXPECTATIONAL EQUATIONS:

% For k(t+1,i)
FF = [ 0, 0
       0, 0 ];

% For k(t)
GG = [ 0, 0
       0, 0 ];

% For k(t-1)
HH = [ 0, 0
       0, 0 ];

% For   c(t+1),y(t+1),x(t+1),c(t+1,i),y(t+1,i),x(t+1,i),r(t+1,i)
JJ = [       0,     0,     0,-eta_1,0,    0, 0,    0, 0,   1, 0   
             0,     0,     0,0,-eta_2,    0, 0,    0, 0,   0, 1  ];
 
% For     c(t),  y(t),  x(t),  c(t,i),  y(t,i),  x(t,i),  r(t,i)
KK = [       0,     0,     0, eta_1,0,    0, 0,    0, 0,    0, 0  
             0,     0,     0, 0,eta_2,    0, 0,    0, 0,    0, 0  ];  

% For z(t+1,i)
LL = [ 0, 0
       0, 0 ];

% For z(t)
MM = [ 0, 0
       0, 0 ];

% AUTOREGRESSIVE MATRIX FOR z(t)

NN = [ psi_1,     0,  
           0, psi_2  ];

Sigma = sigma_eps^2 * ...
         [ 1,       corr_12  
           corr_12,    1     ];

% Setting the options:
  
[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);


PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 4; % Index of output among the variables selected for HP filter
IMP_SELECT = [1:9,12:15];
   %  a vector containing the indices of the variables to be plotted
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.


% Starting the calculations:

do_it;
