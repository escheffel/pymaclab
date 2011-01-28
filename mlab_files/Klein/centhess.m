function hess = centhess(f,x0);
%
% Purpose: To calculate the Hessian of the function f at the point x0
%          using central differences
% Format: 
%
% hess = centhess('f',x0)
%
% where f is a function and f(x) is an mx1 vector
%
% x0 is a nx1 vector
%
% Output: the Hessian, an mnxn matrix defined as in Magnus & Neudecker (1988).
%
%

f0 = feval(f,x0);

m = size(f0,1);
n = size(x0,1);
hess = zeros(m*n,n);
eps = 1e-4;

for i=1:n
   for j = i:n
      h1 = zeros(n,1);
      h1(i)=eps;
      h2 = zeros(n,1);
      h2(j)=eps;
      y = ((feval(f,x0+h1+h2)-feval(f,x0+h1-h2))-(feval(f,x0-h1+h2)-feval(f,x0-h1-h2)))/(4*eps^2);
      for k=1:m
         hess(n*(k-1)+j,i) = y(k);
         if i~=j
            hess(n*(k-1)+i,j) = y(k);
         end
      end
   end
end