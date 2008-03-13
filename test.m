%PARAMETERS

Z_bar     = 1;
alpha     = 0.36;
delta     = 0.025;
betta     = 0.99;
Pi_star   = 1.05;
theta     = 1.0;
rho       = 0.95;
phi       = 0.15;
gamma     = 0.5;
sigma_eps = 0.712;



%STEADY STATES

x0 = zeros(3,1);
x0(1,1) = 10.0;
x0(2,1) = 2.0;
x0(3,1) = 0.3;

k_star = x(1,1);
c_star = x(2,1);
x_star = x(3,1);
n_star = 1-x_star;

y_star = k_star^alpha*(1-x_star)^(1-alpha);
c_star = y_star - delta*k_star;
i_star = 1/betta;
R_star = i_star + Pi_star -1;



%DECLARING VARNAMES 

VARNAMES = ['capital      ',
            'money        ',
            'output       ',
            'consumption  ',
            'investment   ',
            'labour       ',
            'inflation    ',
            'nom_interest ',
            'nom_interest '];



%DETERMINISTIC EQUATIONS:

AA = zeros(8,2);
AA(1,2) = -theta;
AA(3,1) = -alpha*y_bar/k_bar;
AA(5,1) = (betta*alpha)+(phi);
AA(6,1) = 1;
AA(7,2) = 1;
AA(8,2) = -1;

BB = zeros(8,2);
BB(3,1) = (alpha*y_bar/k_bar)*(alpha*(1+eta*n_bar/l_bar))/((alpha+eta*n_bar/l_bar));
BB(4,1) = alpha;
BB(6,1) = -(1-delta);
BB(8,2) = gamma;

CC = zeros(8,8);
CC(1,8) = -1;
CC(1,6) = -(1-theta);
CC(1,7) = (1-theta);
CC(2,4) = -(1+eta*n_bar/(1-n_bar));
CC(2,1) = 1;
CC(2,8) = -1;
CC(3,8) = (alpha*y_bar/k_bar)*((1-alpha))/((alpha+eta*n_bar/l_bar));
CC(3,7) = (-1)+((alpha*y_bar/k_bar)*((1-alpha))/((alpha+eta*n_bar/l_bar)));
CC(4,1) = -1;
CC(4,4) = (1-alpha);
CC(5,1) = -1;
CC(5,3) = (quasi)+(1);
CC(6,3) = -delta;
CC(7,2) = -1;
CC(8,5) = -1;

DD = zeros(8,2);
DD(1,2) = -theta*gamma;
DD(1,1) = -theta*phi;
DD(3,1) = (alpha*y_bar/k_bar)*(rho*(1+eta*n_bar/l_bar))/((alpha+eta*n_bar/l_bar));
DD(4,1) = 1;
DD(8,2) = 1;
DD(9,2) = -1;
DD(10,1) = -1;



%EXPECTATIONAL EQUATIONS:

FF = zeros(2,2);

GG = zeros(2,2);

HH = zeros(2,2);

JJ = zeros(2,8);
JJ(1,8) = -1;
JJ(2,5) = 1;

KK = zeros(2,8);
KK(1,7) = -1;
KK(1,8) = 1;
KK(2,6) = -1;
KK(2,7) = 1;

LL = zeros(2,2);

MM = zeros(2,2);



%AUTOREGRESSIVE MATRIX FOR z(t)

NN = zeros(2,2);
NN(1,2) = -1;
NN(2,1) = -1;



%OPTIONS:

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);

  
PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 4; % Index of output among the variables selected for HP filter
IMP_SELECT = 1:(m_states+n_endog+k_exog);
   %  a vector containing the indices of the variables to be plotted
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
DO_SIMUL   = 1; % Calculates Simulations
DO_MOMENTS = 1; % Calculates Moments
DISPLAY_IMMEDIATELY = 1;
% Starting the calculations:

do_it;
