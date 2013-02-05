from mako.template import Template

template = """\
%Model Description+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This is just a dynare++ to pymaclab translated file which does not possess any model desc.

%Model Information+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Name = dynarepp-to-pymaclab translated model;


%Parameters++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% for lino in paramvals[0]:
  ${lino}
% endfor

%Variable Vectors+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Build the necessary list of variables we need to iter through.
<%
  varli = []
  counter = 1
  for keyo in vardic.keys():
    for i1,itemo in enumerate(vardic[keyo]['var']):
      varli.append([counter, itemo[0], itemo[1], itemo['mod'][i1]])
      counter += 1
%>
[1]  k(t):capital{endo}[log,bk]
[2]  c(t):consumption{con}[log,bk]
[4]  y(t):output{con}[log,bk]
[4]  R(t):rrate{con}[log,bk]
[5]  z(t):eps(t):productivity{exo}[log,bk]
% for keyo in vardic.keys():
  % for itemo in vardic[keyo]['var']
  [${loop.index+1}]  ${vardic[keyo]['var'][0]}:${vardic[keyo]['var'][1]}{keyo}
  % endfor
% endfor

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