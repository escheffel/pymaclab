from mako.template import Template

template = """\
%Model Description+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This is just a dynare++ to pymaclab translated file which does not possess any model desc.

%Model Information+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Name = dynarepp-to-pymaclab translated model;


%Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rho       = 0.36;
delta     = 0.025;
betta     = 1.0/1.01;
eta	  = 2.0; 
psi	  = 0.95;
z_bar     = 1.0;
sigma_eps = 0.052;

%Variable Vectors+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1]  k(t):capital{endo}[log,bk]
[2]  c(t):consumption{con}[log,bk]
[4]  y(t):output{con}[log,bk]
[4]  R(t):rrate{con}[log,bk]
[5]  z(t):eps(t):productivity{exo}[log,bk]

%Boundary Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
None


%Variable Substitution Non-Linear System++++++++++++++++++++++++++++++++++++++++++++++++
None


%Non-Linear First-Order Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Insert here the non-linear FOCs in format g(x)=0

[1]   y(t)-@inv(t)-c(t) = 0;
[2]   betta*(@MU(t+1)/@MU(t))*E(t)|R(t+1)-1 = 0;
[3]   @F(t)-y(t) = 0;
[4]   R(t) - (1+@Fk(t)-delta) = 0;
[5]   LOG(z(t))-psi*LOG(z(t-1))-eps(t) = 0;


%Steady States [Closed Form]+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
None


%Steady State Non-Linear System [Manual]+++++++++++++++++++++++++++++++++++++++++++++++++

[1]   c_bar = 1.0;
[2]   y_bar = 1.0;
[2]   k_bar = 1.0;
[3]   R_bar = 1.01;


%Log-Linearized Model Equations++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
None


%Variance-Covariance Matrix++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Sigma = [sigma_eps**2];


%End Of Model File+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




"""

mako_pml_template = Template(template)