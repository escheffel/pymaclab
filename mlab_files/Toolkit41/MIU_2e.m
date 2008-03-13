clear all

% MIU_e2.M: used to solve MIU function model in C.E.Walsh, "Monetary 
% Theory and Policy," second edition."
disp(' Money in the utility function model,');
disp(' used in Chapter 2 of Monetary Theory and Policy 2nd ed.');
disp(' June 2001');

% Setting parameters;

    % "Technology" parameters
N_bar = 1/3; % steady state employment is a third of total time endowment
alpha = .36; % Capital share
delta = .019; % Depreciation rate for capital

    %"Taste" parameters
beta = .989; 	%discount parameters
PHI = 3; 		% coefficient of relative risk aversion
eta = 1; 		%leisure parameter
a = .95; 		% weight on consumption in composite good definition
					% 1-a is weight on real money balances in utility function.
b =  1/.39; 	% interest elasticity, from Chari, Kehoe, and McGratten

THETA = 1.0125;% 1 plus quarterly rate of nominal money growth

%Exogenous stochastic processes
% z is the technology shock
rho_z = .95;  %autocorrelation of technology shock
sigma_z = .7; %Standard deviation of technology shock innovation. Units: percent

%Exogenous policy specification
% money growth g(t) = rho_g*g(t-1) + phi*z(t-1) + v(t)
sigma_m = .89; %Standard deviation of money growth rate shock. Units: Percent
rho_g = 0.5; 		%Autocorrelation of money shock
phi = 0; 		%Correlation of technology shock and money shock
sigma_v = ((1-rho_g^2) * (sigma_m)^2 - (phi^2/(1-rho_z^2))*sigma_z^2)^.5;


%Calculating the steady state:
R_bar = 1/beta; %one percent real interest per quarter
YK_bar = (R_bar - 1 + delta)/alpha;
CK_bar = YK_bar - delta;
NK_bar = (YK_bar)^ (1/(1-alpha));
inf_bar = THETA - 1;
NL_bar = N_bar/(1.0 - N_bar);
YN_bar = YK_bar /(1.0 - N_bar);
I_bar = R_bar*THETA;
i_bar = I_bar - 1;
mK_bar = (CK_bar) * (a*i_bar/((1-a) * (1+i_bar)))^ (-1/b);
mC_bar = mK_bar /CK_bar;
X = a* (CK_bar^ (1-b)) + (1-a) * (mK_bar^(1-b));
XX = (PHI * a * ((CK_bar) ^ (1-b)) + (1-a) * b*((mK_bar) ^ (1-b)))/X;
XM = (1-a) * (mK_bar^(1-b))/X;
gamma = 1/(1 + ((1-a)/a)*((mC_bar)^(1-b)));
A = gamma* PHI + (1-gamma)*b;
B = (b - PHI)*(1-gamma);

Parameters = ['alpha     ',
              'delta     ',
              'beta      ',
              'PHI       ',
              'eta       ',
              'a         ',
              'b         ',
              'THETA     ',
              'rho_z     ',
              'sigma_z   ',
              'rho_g     ',
              'phi       ',
              'sigma_m   ',
              'sigma_v   '];
           
values = [alpha, delta, beta, PHI, eta, a, b, THETA, rho_z, sigma_z, rho_g, phi, sigma_m, sigma_v];
% printing parameter values
Parameters
values'
%printing steady_state value
R_bar
YK_bar
CK_bar
mK_bar
mC_bar
NK_bar
inf_bar
sigma_v

%Declaring the matries.

VARNAMES = ['capital        ',
            'money balances ',
            'output         ',
            'consumption    ',
            'employment     ',
            'lambda         ',
            'nominal rate   ',
            'real rate      ',
            'inflation      ',
            'investment     ',
            'technology     ',
            'money growth   '];
    
%Translating into coefficient matries.

%Check: 12 variables.
%Endogenous state variables "x(t)": k(t), m(t)
%Endonenous other variables "y(t)": y(t), c(t), n(t), lambda(t), i(t), r(t), inf(t), x(t).
%Exogenous state variables "z(t)": z(t), g(t).
%Switch to that notation. Find matrices for format
%0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) ++ KK y(t) = LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon (t+1) with E_t [ epsilon(t+1) ] = 0,

%equation order
%1). resource constraint
%2). production function
%3). money demand equation
%4). evolution of m
%5). labor market
%6). euler equation 1 
%7). euler equation 2
%8). investment evolvement equation
%8). fisher equation
%9). ex-post real interest

% for k(t), m(t):
k_7 = alpha *(alpha-1)* YK_bar*(eta*NL_bar);

AA = [1, 0
    0, 0
    0, 1
    0, 1
    0, 0
    0, -B
    k_7, 0
    1, 0];

%for k(t-1), m(t-1):
BB = [ (delta -1),  0
    alpha, 0
    0, 0
    0, -1
    0, 0
    0, 0
    0, 0
    (delta-1), 0];

%order: y(t), c(t), n(t), lambda(t), i(t), r(t), inf(t), x(t)
CC = zeros(8, 8);
CC(1, 1) = -YK_bar;
CC(1, 2) = CK_bar;
CC(2, 1) = -1;
CC(2, 3) = 1-alpha;
CC(3, 2) = -1;
CC(3, 5) = 1/b;
CC(4, 7) = 1;
CC(5, 1) = 1;
CC(5, 3) = -(1 + eta*NL_bar);
CC(5, 4) = 1;
CC(6, 2) = A;
CC(6, 4) = 1;
CC(7, 4) = alpha*(1-alpha)*YK_bar;
CC(7, 6) = -(alpha + eta *NL_bar + alpha*(1-alpha)*YK_bar);
CC(8, 8) = -delta;

DD = zeros(8, 2);
DD(2, 1) = 1;
DD(4, 2) = -1;
DD(7, 1) = alpha*YK_bar*(1+eta*NL_bar)*rho_z;


FF = zeros(2, 2);

GG = zeros(2, 2);

HH = zeros(2, 2);

JJ = zeros(2, 8);
JJ(1, 7) = 1;
JJ(2, 4) =1;

KK = [0, 0, 0, 0, -1, 1, 0, 0
    0, 0, 0, -1, 0, 1, 0, 0];

LL = [0, 0
    0, 0];

MM = [0, 0
    0, 0];

NN = [rho_z, 0
    phi, rho_g];

Sigma = [sigma_z^2 0
    0 sigma_v^2];

%setting the options:
[l_equ, m_states] = size(AA);
[l_equ, n_endog] = size(CC);
[l_equ, k_exog] = size(DD);

HORIZON = 32; %how far out should the inpulse response be calculated
PERIOD = 4; %number of periods per year, i.e. 12 for monthly, 4 for quarterly
DO_PLOTS = 1; %if impulse response plots should be made, = 0, if not
% IMP_SELECT is a vector containing the indices of the variables to be plotted
IMP_SELECT = [1, 2, 3, 5, 7, 8, 9]; % for money growth shock
%IMP_SELECT = [1, 3, 4, 5, 8]; % for money growth shock -- real variables
%IMP_SELECT = [1, 3, 5, 10];  % for technology shock
%IMP_SELECT = 1:(m_states+n_endog+k_exog); % for plotting all variables
HP_SELECT = 1: (m_states + n_endog + k_exog); %selecting the variables for the HP filter calcs.
GNP_INDEX = 3; %Index of output among the variables selected for HP filter
DO_MOMENTS = 1;
DISPLAY_AT_THE_END = 0;
%Start the calculations;
do_it;