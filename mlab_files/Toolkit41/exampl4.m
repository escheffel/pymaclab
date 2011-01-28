% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL4.M calculates through the models in
% Farmer, R.E.A. and J.T. Guo, "Real Business Cycles and the 
% Animal Spirits Hypothesis," Journal of Economic Theory 63, 42-72 (1994).
% You can choose the model by setting FARMER_GUO = 1, 2 or 3.
% The default is FARMER_GUO = 3 to show the failure.
% 
% First, parameters are set and the steady state is calculated. Next, the matrices are
% declared.  In the last line, the model is solved and analyzed by calling DO_IT.M

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

disp('EXAMPLE 4: this model is taken from Farmer, R.E.A.');
disp('           and J.T. Guo, "Real Business Cycles and the Animal');
disp('           Spirits Hypothesis," Journal of Economic Theory');
disp('           63, 42-72 (1994).');
disp('           The model by Farmer and Guo is a discrete time version');
disp('           of the model by Benhabib and Farmer, same JET issue.');
disp('           You can choose the model by setting FARMER_GUO = 1, 2 or 3');
disp('           The default is FARMER_GUO = 3.');
disp('   ==== NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE ===');
disp('   Example 4 is meant to be instructive: it is set up wrong for');
disp('   the third model, although not for the other two models!');
disp('   In this version of the example, the solution gets stuck with model 3:');
disp('   there are two imaginary roots, but only one state variable.');
disp('   Thus, the method cannot proceed further, although the');
disp('   algorithm attempts to do so buy cutting off the imaginary parts.');
disp('   To do this example right, check EXAMPL5.M!!!!');
disp('   === END OF NOTE ==========================================');

disp('Hit any key when ready...');
pause;

% Setting parameters:

if exist('FARMER_GUO')~=1,
   FARMER_GUO = 3;
end;
if FARMER_GUO == 1,  % Model 1
   lambda    = 1.00; % see Table I in Farmer-Guo:
                     % Model 1: = 1, Model 2: = 0.7, Model 3: = 0.58
   b         = 0.64; % Labor share, see Table I in Farmer-Guo:
                     % Model 1: = 0.64, Model 2: = 0.63, Model 3: = 0.7
   a         = 0.36; % Capital share, see table II in Farmer-Guo:
                     % Model 1: = 0.36, Model 2: = 0.3, Model 3: = 0.23                  
   psi       =1.0/1.05; % autocorrelation of technology shock
                     % = 1/theta in Farmer-Guo, see table III.  
                     % Model 1 and 2: = 1/1.05, model 3: N/A,
                     % For model 3, we leave it in for the heck of it anyhow,
                     % but we make the aver. shock size sigma_eps really tiny.              
   sigma_eps = 0.7;  % Percentage standard deviation of technology shock.
                     % Model 1: = 0.7, Model 2: = 0.465, model 3: N/A.
   psi_nu    = 0;    % autocorrelation of sunspot.  Only in model 3.
   sigma_nu  = 0.001;% Standard deviation of technology shock.  
                     % Units: Percent. Check Table IV in Farmer-Guo:
                     % Model 1: N/A, Model 2: N/A, model 3: = 0.217.
elseif FARMER_GUO == 2,   % Model 2
   lambda    = 0.7;
   b         = 0.63;
   a         = 0.3;           
   psi       = 1.0/1.05;            
   sigma_eps = 0.465;
   psi_nu    = 0;    
   sigma_nu  = 0.001;
else   % Model 3
   lambda    = 0.58; 
   b         = 0.7;  
   a         = 0.23;              
   psi       =1.0/1.05; % autocorrelation of technology shock
                     % = 1/theta in Farmer-Guo, see table III.  
                     % Model 1 and 2: = 1/1.05, model 3: N/A,
                     % For model 3, we leave it in for the heck of it anyhow,
                     % but we make the aver. shock size sigma_eps really tiny.              
   sigma_eps = .001; 
   psi_nu    = 0;    % autocorrelation of sunspot.  Only in model 3.
   sigma_nu  = 0.217;
end;

L_bar     = 1.0/3;% Steady state employment is a third of total time endowment
Z_bar     = 1;    % Normalization
delta     = 0.025;% Depreciation rate for capital
rho       = 0.99; % Discount rate per quarter.  
                  % NOTATION: FARMER-GUO.  So, this is NOT
                  % the capital share, as in most other models.
eta       = 1.0;  % constant of relative risk aversion 
                  % = 1/(coeff. of intertemporal substitution)
                  % Farmer-Guo only consider the case eta = 1;
                  

% Calculating the steady state:

alpha   = a/lambda;
betta   = b/lambda;  % NOTATION: FARMER-GUO.
                     % THIS IS THE ELASTICITY OF OUTPUT WITH RESPECT
                     % TO LABOR, NOT THE DISCOUNT FACTOR, AS USUAL.

R_bar   = 1.0/rho;  % Steady state return
YK_bar  = (R_bar + delta - 1)/a;  % = Y_bar / K_bar
K_bar   = (YK_bar / (Z_bar*L_bar^betta) )^(1.0/(alpha-1));
I_bar   = delta * K_bar;
Y_bar   = YK_bar * K_bar;
C_bar   = Y_bar - delta*K_bar;
A       = C_bar^(-eta) * b * Y_bar/L_bar; % Parameter in utility function

% Declaring the matrices. 


VARNAMES = ['capital    ',
            'consumption',
            'output     ',
            'labor      ',
            'interest   ',
            'investment ',
            'technology ',
            'sunspot    '];

% Translating into coefficient matrices.  
% The equations are, conveniently ordered:
% 1) 0 = - I i(t) - C c(t) + Y y(t)
% 2) 0 = I i(t) - K k(t) + (1-delta) K k(t-1)
% 3) 0 = alpha k(t-1) - y(t) + betta l(t) + z(t)
% 4) 0 = -eta c(t) + y(t) - l(t)
% 5) 0 = - a Y/K k(t-1) + a Y/K y(t) - R r(t)
% 6) 0 = E_t [ - eta c(t+1) + r(t+1) + eta c(t) ]
% 7) z(t+1) = psi z(t) + epsilon(t+1)
% 8) s(t+1) = eps_s(t+1): iid sunspot
% CHECK: 8 equations, 8 variables.
%
% Endogenous state variables "x(t)": k(t)
% Endogenous other variables "y(t)": c(t), y(t), l(t), r(t), i(t)
% Exogenous state variables  "z(t)": z(t).
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

% for k(t):
AA = [ 0
       - K_bar
       0
       0
       0 ];

% for k(t-1):
BB = [ 0
       (1-delta)*K_bar
       alpha
       0
       - a * YK_bar ];

%Order:   consumption  output      labor     interest  investment
CC = [    -C_bar,      Y_bar,      0,        0,        -I_bar % Equ. 1)
          0,           0,          0,        0,        I_bar  % Equ. 2)
          0,           -1,         betta,    0,        0      % Equ. 3)      
          -eta,        1,          -1,       0,        0      % Equ. 4)
          0,           a*YK_bar,   0,        - R_bar,  0 ];   % Equ. 5)

% DD = [ 0, 0
%        0, 0
%        1, 0
%        0, 0
%        0, 0 ];

DD = [ 0
       0
       1
       0
       0 ];

FF = [ 0 ];

GG = [ 0 ];

HH = [ 0 ];

JJ = [ -eta,  0,  0,  1,  0];

KK = [ eta,   0,  0,  0,  0];

% LL = [ 0, 0 ];
LL = [ 0 ];

% MM = [ 0, 0 ];
MM = [ 0 ];

% NN = [psi    0
%        0    psi_nu ];
NN = [ psi ];

% Sigma = [ sigma_eps^2,    0
%               0,       sigma_nu^2      ];
Sigma = [ sigma_eps^2 ];

% Setting the options:
  
[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);


PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 3; % Index of output among the variables selected for HP filter
IMP_SELECT = [1:5,7];
   %  a vector containing the indices of the variables to be plotted
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.



% Starting the calculations:

do_it;


if FARMER_GUO == 1,
  disp('EXAMPL4.M: Model 1 of Farmer-Guo has been calculated.');
elseif FARMER_GUO == 2,
  disp('EXAMPL4.M: Model 2 of Farmer-Guo has been calculated.');
else
  disp('EXAMPL4.M: Model 3 of Farmer-Guo has been calculated.');
end;

 
       

