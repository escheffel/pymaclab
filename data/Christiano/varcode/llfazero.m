function [llf,cnsts]=llfazeroa(x,sigma,FUNA);   
%LLFAZEROA
%computes the likelihood function llf and a vector cnsts that is
%used by constr to constrain the matrix A0 to have positive diagonal elements. 
%x is a vector whose values are the matrix A0.
%sigma is the variance covariance matrix of the data.
%FUNA user supplied procedure to calculate the A0 matrix.
%returns the matrix a0 that it makes out of the vector x

%Note that I haven't scaled this
%by the number of observations or included any constant terms. 
%this however does produce a simple function

[a0]=feval(FUNA,x);

%a few calculations
bt = a0*sigma*a0';
tra = sum(diag(bt));
llf = -log(abs(det(a0))) + 0.5*tra;       

%this constraint the diagonal elements to be positive
%Additional  constraints can be added. See 
cnsts=[-1*diag(a0)];    
   