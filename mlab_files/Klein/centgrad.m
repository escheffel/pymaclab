function grdd = centgrad(f,x0);
%
% Purpose: To calculate the gradient of the function f at the point x0
%          using central differences
% Format: 
%
% The format is
%
% gradient = grad('f',x0)
%
% where f is a function and f(x) is an mx1 vector
%
% x0 is a nx1 matrix
%
% Output: the gradient, an mxn matrix defined as in Magnus & Neudecker (1988).
%
%

f0 = feval(f,x0);
m = size(f0,1);
n = size(x0,1);
grdd = zeros(m,n);
eps = 5e-6;

for i=1:n
   h = zeros(n,1);
   h(i)=1;
   grdd(:,i) = (feval(f,x0+eps*h)-feval(f,x0-eps*h))/(2*eps);
end


