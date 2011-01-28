% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL3.M calculates through an extensions of Hansens benchmark real business
% cycle model with a simple time-to-decay feature, leading to echo effects.  
% The change to Hansens model consists in assuming that capital
% does not depreciate at all for p periods, and then depreciates 100 percent.  We set p to 4.
% This example was prepared to discuss a paper by Boucekkine, Germain and Licandro

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.


% Setting parameters:

N_bar     = 1.0/3;  % Steady state employment is a third of total time endowment
Z_bar     = 1; % Normalization
rho       = .36; % Capital share
R_bar     = 1.01; % One percent real interest per quarter
eta       = 1.0; % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
psi       = .95; % autocorrelation of technology shock
sigma_eps = .712; % Standard deviation of technology shock.  Units: Percent.
p_echo    = 4;

disp('EXAMPLE 3: Hansen benchmark real business cycle model ');
disp('           with a time-to-decay feature, leading to echo effects.');
disp(sprintf('          New capital does not depreciate for %2.0f periods',p_echo));
disp('          but completely depreciates after that');

disp('Hit any key when ready...');
pause;





% Calculating the steady state:

betta   = 1.0/R_bar;  % Discount factor beta
YK_bar  = (1- betta)/((1 - betta^p_echo)*betta*rho); % = Y_bar / K_bar
K_bar   = (YK_bar / Z_bar)^(1.0/(rho-1)) * N_bar;
I_bar   = K_bar / p_echo ;
Y_bar   = YK_bar * K_bar;
C_bar   = Y_bar - I_bar;
Lam_bar = C_bar^(- eta); % Lambda bar, steady state marg. utility of cons.
Mu_bar  = rho*Lam_bar*YK_bar;  % Lagrangian on capital = sum(j=1..p_echo) of lagged investments
A       = Lam_bar * (1 - rho) * Y_bar/N_bar; % Parameter in utility function

% Declaring the matrices. 


VARNAMES = ['investment      ',  % 1
            'investment(t-1) ',  % 2
            'investment(t-2) ',  % 3
            'investment(t-3) ',  % 4
            'E_t[mu(t+2)]    ',  % 5
            'E_t[mu(t+3)]    ',  % 6
            'E_t[mu(t+4)]    ',  % 7
            'consumption     ',  % 8
            'output          ',  % 9
            'capital         ',  % 10
            'labor           ',  % 11
            'marginal utility',  % 12
            'mu              ',  % 13
            'E_t[mu(t+1)]    ',  % 14
            'Solow parameter ']; % 15

% Translating into coefficient matrices.  Let i(t,j)=i(t-j),j=0,..,p-1, and mu_e(t,j) = E_t[mu(t+j)],j=0,..,p
% The equations are, conveniently ordered:
% 1)  0 = - I i(t,0) - C c(t) + Y y(t)
% 2)  0 =  i(t-1,0) + i(t-1,1) + i(t-1,2) + i(t-1,3) - n k(t) 
% 3)  0 = rho k(t) - y(t) + (1-rho) n(t) + z(t)
% 4)  0 = lambda(t) + y(t) - n(t)
% 5)  0 = (-Lam_bar/Mu_bar) lambda(t) + betta mu_e(t,1) + betta^2 mu_e(t,2) + betta^3 mu_e(t,3) + betta^4 mu_e(t,4)
% 6)  0 = lambda(t) + y(t) - k(t) - mu_e(t,0)
% 7)  0 = lambda(t) + eta c(t)
% 8)  0 = i(t,1)-i(t-1,0) = 0
% 9)  0 = i(t,2)-i(t-1,1) = 0
% 10) 0 = i(t,3)-i(t-1,2) = 0
% 11) 0 = E_t [ mu_e(t,1) - mu_e(t+1,0) ]
% 12) 0 = E_t [ mu_e(t,2) - mu_e(t+1,1) ]
% 13) 0 = E_t [ mu_e(t,3) - mu_e(t+1,2) ]
% 14) 0 = E_t [ mu_e(t,4) - mu_e(t+1,3) ]
% 15) z(t+1) = psi z(t) + epsilon(t+1)
% CHECK: 15 equations, 15 variables: c(t), y(t), k(t), n(t), lambda(t), mu_e(t,j), j=0,..,p, i(t,j), j=0,..,p-1, z(t)
%
% Endogenous state variables "x(t)": i(t,j), j=0,..,p-1, mu_e(t,j), j = 2,..,p
%    Here, the mu_e(t,j) are "artifically" turned into state variables to avoid
%    that the matrix CC is of too low rank.
% Endogenous other variables "y(t)": c(t), y(t), k(t), n(t), lambda(t), mu_e(t,j), j=0,1,
% Exogenous state variables  "z(t)": z(t).
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,


% EQUATIONS (1) THROUGH (10):

% for i(t,j), j = 0, .., p-1,  mu_e(t,j), j = 2,..,p
AA = [ -I_bar, 0, 0, 0,             0,        0,        0   % Equ. 1)
            0, 0, 0, 0,             0,        0,        0   % Equ. 2)
            0, 0, 0, 0,             0,        0,        0   % Equ. 3)
            0, 0, 0, 0,             0,        0,        0   % Equ. 4)
            0, 0, 0, 0,       betta^2,  betta^3,  betta^4   % Equ. 5)
            0, 0, 0, 0,             0,        0,        0   % Equ. 6)
            0, 0, 0, 0,             0,        0,        0   % Equ. 7)
            0, 1, 0, 0,             0,        0,        0   % Equ. 8)
            0, 0, 1, 0,             0,        0,        0   % Equ. 9)
            0, 0, 0, 1,             0,        0,        0 ];% Equ. 10)

% for i(t-1,j), j = 0, .., p-1,  mu_e(t-1,j), j = 2,..,p
BB = [      0, 0, 0, 0,             0,        0,        0   % Equ. 1)
            1, 1, 1, 1,             0,        0,        0   % Equ. 2)
            0, 0, 0, 0,             0,        0,        0   % Equ. 3)
            0, 0, 0, 0,             0,        0,        0   % Equ. 4)
            0, 0, 0, 0,             0,        0,        0   % Equ. 5)
            0, 0, 0, 0,             0,        0,        0   % Equ. 6)
            0, 0, 0, 0,             0,        0,        0   % Equ. 7)
           -1, 0, 0, 0,             0,        0,        0   % Equ. 8)
            0,-1, 0. 0,             0,        0,        0   % Equ. 9)
            0, 0,-1, 0,             0,        0,        0 ];% Equ. 10)

% for endogeneous other variables
%Order:   consumption,output,capital,labor,lambda,mu_e(t,0),mu_e(t,1)
CC = [         -C_bar, Y_bar,      0,    0,     0,        0,        0  % Equ. 1)
                    0,     0,-p_echo,    0,     0,        0,        0  % Equ  2) 
                    0,    -1,    rho,(1-rho),   0,        0,        0  % Equ  3)            
                    0,     1,      0,   -1,     1,        0,        0  % Equ  4)    
                    0,     0,      0,0,(-Lam_bar/Mu_bar), 0,    betta  % Equ  5)           
                    0,     1,     -1,    0,     1,       -1,        0  % Equ  6) 
                  eta,     0,      0,    0,     1,        0,        0  % Equ  7)     
                    0,     0,      0,    0,     0,        0,        0  % Equ  8)           
                    0,     0,      0,    0,     0,        0,        0  % Equ  9)           
                    0,     0,      0,    0,     0,        0,        0];% Equ  10)           
 
DD = [ 0
       0
       1
       0
       0
       0
       0
       0
       0
       0 ];

% EQUATIONS (11) THROUGH (14):


FF = [ 0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,   -1, 0, 0
       0, 0, 0, 0,    0,-1, 0 ] ;

GG = [ 0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    1, 0, 0
       0, 0, 0, 0,    0, 1, 0
       0, 0, 0, 0,    0, 0, 1 ] ;

HH = [ 0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    0, 0, 0 ] ;




%Order:   consumption,output,capital,labor,lambda,mu_e(t,0),mu_e(t,1)

JJ = [              0,     0,      0,    0,     0,       -1,        0  % Equ  11)    
                    0,     0,      0,    0,     0,        0,       -1  % Equ  12)           
                    0,     0,      0,    0,     0,        0,        0  % Equ  13)           
                    0,     0,      0,    0,     0,        0,        0];% Equ  14)           

KK = [              0,     0,      0,    0,     0,        0,        1  % Equ  11)    
                    0,     0,      0,    0,     0,        0,        0  % Equ  12)           
                    0,     0,      0,    0,     0,        0,        0  % Equ  13)           
                    0,     0,      0,    0,     0,        0,        0];% Equ  14)           

LL = [ 0
       0
       0
       0 ];

MM = [ 0
       0
       0
       0 ];

NN = [psi];

Sigma = [ sigma_eps^2  ];

% Setting the options:

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);

 
PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 2; 
% IMP_SELECT = 1:(m_states+n_endog+k_exog);
   %  a vector containing the indices of the variables to be plotted
IMP_SELECT = [ 1, 8, 9, 10, 11, 15];  % excluding investment in the impulse-response plots
% HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
HP_SELECT = [1, 9, 15 ];


% Starting the calculations:

do_it;


 
       

