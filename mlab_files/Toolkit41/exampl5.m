% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL5.M calculates through the models in
% Farmer, R.E.A. and J.T. Guo, "Real Business Cycles and the 
% Animal Spirits Hypothesis," Journal of Economic Theory 63, 42-72 (1994).
% Here the problem is set up right also for model 3
% by simply making c(t) into an additional state variable as well.
% You can choose the model by setting FARMER_GUO = 1, 2 or 3.
% The default is FARMER_GUO = 3, since that is the "weirdest" case.
% 
% First, parameters are set and the steady state is calculated. Next, the matrices are
% declared.  In the last line, the model is solved and analyzed by calling DO_IT.M

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

disp('EXAMPLE 5: this model is taken from Farmer, R.E.A.');
disp('           and J.T. Guo, "Real Business Cycles and the Animal');
disp('           Spirits Hypothesis," Journal of Economic Theory');
disp('           63, 42-72 (1994).');
disp('           The model by Farmer and Guo is a discrete time version');
disp('           of the model by Benhabib and Farmer, same JET issue.');
disp('           You can choose the model by setting FARMER_GUO = 1, 2 or 3');
disp('           The default is FARMER_GUO = 3.');
disp('   ==== NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE ===');
disp('   Example 5 is set up correctly in contrast to exampl4.m:');
disp('   we have simply introduced an additional state variable c(t)');
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
   sigma_eps = 0.7;  % Percentage standard deviation of technology shock.
                     % Model 1: = 0.7, Model 2: = 0.465, model 3: N/A.
elseif FARMER_GUO == 2,   % Model 2
   lambda    = 0.7;
   b         = 0.63;
   a         = 0.3;           
   psi       = 1.0/1.05;            
   sigma_eps = 0.465;
else   % Model 3
   lambda    = 0.58; 
   b         = 0.7;  
   a         = 0.23;              
   psi       =1.0/1.05; % autocorrelation of technology shock
                     % = 1/theta in Farmer-Guo, see table III.  
                     % For Model 3: N/A, but we leave it in for the heck
                     % of it to see what it would do in model 3.
   sigma_eps = 0.001;% Percentage standard deviation of technology shock.
                     % we make it really small, because it should be = 0 for model 3.
   psi_nu    = 0;    % autocorrelation of sunspot.  Only in model 3.
   sigma_nu  = 0.217;% Standard deviation of technology shock.  
                     % Units: Percent. Check Table IV in Farmer-Guo:
                     % Model 1: N/A, Model 2: N/A, model 3: = 0.217.
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

SUN_SCALE = 2.0969;  % This was computed from a previous run of this example
        % to make sure that a one-percent move in the sunspot variable
        % results in a one-percent move of consumption in model 3.
        % This is easy to do: just run the program first with FARMER_GUO = 3,
        % and SUN_SCALE = 1,and then set SUN_SCALE =  1/QQ(2,2)

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
% 1) 0 = - I i(t) - C c(t) + Y y(t),            for models 1 and 2
% 1) 0 = - I i(t) - C c(t) + Y y(t) + SUN_SCALE s(t),   for model 3
%    In model 3, c(t) is moved by the sunspot variable s(t) "above and
%    beyond" what would be predicted by past state variables.
% 2) 0 = I i(t) - K k(t) + (1-delta) K k(t-1)
% 3) 0 = alpha k(t-1) - y(t) + betta l(t) + z(t)
% 4) 0 = -eta c(t) + y(t) - l(t)
% 5) 0 = - a Y/K k(t-1) + a Y/K y(t) - R r(t)
% 6) 0 = E_t [ - eta c(t+1) + r(t+1) + eta c(t) ]
% 7) z(t+1) = psi z(t) + epsilon(t+1)
% 8) s(t+1) = eps_s(t+1): iid sunspot, only for model 3
% CHECK: 8 equations, 8 variables.
%
% Endogenous state variables "x(t)": k(t), c(t)
% Endogenous other variables "y(t)": y(t), l(t), r(t), i(t)
% Exogenous state variables  "z(t)": z(t), s(t)
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

% for k(t) and c(t)
AA = [ 0,       -C_bar      
       - K_bar, 0
       0,       0
       0,       -eta
       0,       0   ];

% for k(t-1) and c(t-1)
BB = [ 0,                  0
       (1-delta)*K_bar,    0
       alpha,              0
       0,                  0
       - a * YK_bar,       0 ];

%Order:   output      labor     interest  investment
CC = [     Y_bar,      0,        0,        -I_bar % Equ. 1)
           0,          0,        0,        I_bar  % Equ. 2)
           -1,         betta,    0,        0      % Equ. 3)      
           1,          -1,       0,        0      % Equ. 4)
           a*YK_bar,   0,        - R_bar,  0 ];   % Equ. 5)

if FARMER_GUO < 3,
   DD = [ 0
          0
          1
          0
          0 ];
else
   DD = [ 0, SUN_SCALE
          0, 0
          1, 0
          0, 0
          0, 0 ];
end;

FF = [ 0, -eta];

GG = [ 0, eta ];

HH = [ 0, 0 ];

JJ = [  0,  0,  1,  0];

KK = [  0,  0,  0,  0];

 
if FARMER_GUO < 3,
   LL = [ 0 ] ;   
   MM = [ 0 ];
   NN = [ psi ];
   Sigma = [sigma_eps^2 ];
else
   LL = [ 0, 0 ];
   MM = [ 0, 0 ];
   NN = [psi    0
          0    psi_nu ];
   Sigma = [ sigma_eps^2,    0
             0,       sigma_nu^2      ];
end;

% Setting the options:
PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 3; % Index of output among the variables selected for HP filter
if FARMER_GUO < 3,
   IMP_SELECT = [1:5,7];
   HP_SELECT = [1:7];
else
   IMP_SELECT = [1:5,7,8];
   HP_SELECT = [1:8];
end;


% Starting the calculations:

do_it;

if FARMER_GUO == 1,
  disp('EXAMPL5.M: Model 1 of Farmer-Guo has been calculated.');
elseif FARMER_GUO == 2,
  disp('EXAMPL5.M: Model 2 of Farmer-Guo has been calculated.');
else
  disp('EXAMPL5.M: Model 3 of Farmer-Guo has been calculated.');
end;
 
       

