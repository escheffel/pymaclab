function VCVK = vcv(VCVX, G, H);
%
% ****************************************************************
%                 VARIANCE-COVARIANCE MATRIX OF K
% ****************************************************************
%
% ****************************************************************
% This function computes the stationary variance-covariance matrix
% of k(t), VCVK, given a variance covariance matrix VCVX of x(t),
% since, according to the model solution:
%
% k(t+1) = G*k(t) + H*x(t)
%
% The computation follows Hamilton (1994), pp. 264-6.
% ****************************************************************

NK = size(G,1);

HVCVXH = H*VCVX*H';

vecHVCVXH = HVCVXH(:);

vecVCVK = inv(eye(NK*NK) - kron(G, G))*vecHVCVXH;

VCVK = reshape(vecVCVK, NK, NK);

