%
% Function: solab2
%
% Purpose: Solves for the second-order accurate approximation 
%          of the solution to a dynamic model
%
% The model is E_t[f(x_{t+1},y_{t+1},x_t,y_t)]=0
%          
% Inputs: gra, an (nx+ny) by 2*(nx+ny) matrix (the gradient of f)
%         hes, an 2*(nx+ny)^2 by 2*(nx+ny) matrix (the Hessian of f)
%         ssigma, an nx by nx matrix; the variance of eps_{t+1}:=x_{t+1}-E_t[x_{t+1}]
%         nx, a natural number
%
% For the definition of the gradient and the Hessian, see
% Magnus and Neudecker: Matrix Differential Calculus with 
% Applications in Statistics and Econometrics.
%
% Outputs: the matrices ff, pp, ee and gg and vectors kx and ky
%          whose roles are described below
%
% x_{t+1} = kx + Px_t + (1/2)*kron(eye(nx),x'_t)Gx_t + eps_{t+1}
%
% y_{t}   = ky + Fx_t + (1/2)*kron(eye(ny),x'_t)Ex_t
%
% Calls: solab, tracem
%
% Code written by Paul Klein in May 2005
% Code debugged and improved in July 2005 as suggested by Paul Gomme
% Code debugged and comments modified in November 2006
% New modifications in November 2006
%
% For details see Gomme and Klein: Second-order approximations of dynamic
%                                  models without the use of tensors
% 

function [ff,pp,ee,gg,kx,ky] = solab2(gra,hes,ssigma,nx);

m = size(gra,1);   % Number of variables
ny = m-nx;			 % Number of non-state variables

aa = gra(:,1:m);
bb = -gra(:,m+1:end);

[ff,pp] = solab(aa,bb,nx);  % Linear approximation

f1 = gra(:,1:nx);
f2 = gra(:,nx+1:m);
f4 = gra(:,m+nx+1:2*m);

mm = [pp;ff*pp;eye(nx);ff];
aa1 = kron(eye(m),mm')*hes*mm;
bb1 = kron(f1,eye(nx));
bb2 = kron(f2,eye(nx));
bb4 = kron(f4,eye(nx));
cc1 = kron(eye(ny),pp');
cc2 = kron(ff,eye(nx));

aa = [kron(eye(nx),bb4)+kron(pp',bb2*cc1) kron(eye(nx),bb1+bb2*cc2)];
sol = -aa\aa1(:);
ee = reshape(sol(1:nx^2*ny),nx*ny,nx);
gg = reshape(sol(nx^2*ny+1:end),nx^2,nx);

ma = 2*[f1+f2*ff f2+f4];
eyeff = [eye(nx);ff;zeros(m,nx)];
ve = f2*tracem(kron(eye(ny),ssigma)*ee)+...
     tracem(kron(eye(m),eyeff')*hes*eyeff*ssigma);
kxy = -ma\ve;
kx = kxy(1:nx);
ky = kxy(nx+1:end);
