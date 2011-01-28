%
% Function: solab
%
% Written by Paul Klein
%
% Rewritten in November 2007 after a suggestion by Alexander Meyer-Gohde
%
% Format: [f,p] = solab(a,b,nk);
%
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations.
%
% Inputs: Two square matrices a and b and a natural number nk
%
% a and b are the coefficient matrices of the difference equation
%
% a*x(t+1) = b*x(t)
% 
% where x(t) is arranged so that the state variables come first, and
%
% nk is the number of state variables.
%
% Outputs: the decision rule f and the law of motion p. If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)   = f*k(t) and
%
% k(t+1) = p*k(t).
%
% Calls: qz, ordqz
%

function [f,p,z11] = solab(a,b,nk);

[s,t,q,z] = qz(a,b);                % upper triangular factorization of the matrix pencil b-za
[s,t,q,z] = ordqz(s,t,q,z,'udo');   % reordering of generalized eigenvalues with the block inside the unit circle in the upper left
z21 = z(nk+1:end,1:nk);
z11 = z(1:nk,1:nk);

if rank(z11)<nk;
	error('Invertibility condition violated')
end

z11i = z11\eye(nk);
s11 = s(1:nk,1:nk);
t11 = t(1:nk,1:nk);

if abs(t(nk,nk))>abs(s(nk,nk)) | abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1));
   warning('Wrong number of stable eigenvalues.');
end

dyn = s11\t11;

f = real(z21*z11i);
p = real(z11*dyn*z11i);
