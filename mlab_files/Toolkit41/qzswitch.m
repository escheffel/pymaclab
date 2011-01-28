function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, interchanges
% diagonal elements i and i+1 of both A and B, while maintaining 
% Q'AZ' and Q'BZ' unchanged.  Does nothing if ratios of diagonal elements
% in A and B at i and i+1 are the same.  Aborts if diagonal elements of
% both A and B are zero at either position.
%

% Copyright: C.A. Sims, 1996, Yale University.

a = A(i,i); d = B(i,i); b = A(i,i+1); e = B(i,i+1);
c = A(i+1,i+1); f = B(i+1,i+1); 
wz = [c*e-f*b, (c*d-f*a)'];
xy = [(b*d-e*a)', (c*d-f*a)'];
n = sqrt(wz*wz');
m = sqrt(xy*xy');
if n == 0
   return
else
   wz = n\wz;
   xy = m\xy;
   wz = [wz; -wz(2)', wz(1)'];
   xy = [xy;-xy(2)', xy(1)'];
   A(i:i+1,:) = xy*A(i:i+1,:);
   B(i:i+1,:) = xy*B(i:i+1,:);
   A(:,i:i+1) = A(:,i:i+1)*wz;
   B(:,i:i+1) = B(:,i:i+1)*wz;
   Z(:,i:i+1) = Z(:,i:i+1)*wz;
   Q(i:i+1,:) = xy*Q(i:i+1,:);
end
