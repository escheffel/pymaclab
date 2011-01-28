function ACM = acm(LAG, VCVX, VCVK, D, F, G, H);
%
% ****************************************************************
%                 AUTOCOVARIANCE MATRIX OF Y
% ****************************************************************
%
% ****************************************************************
% This function computes ACM = E[y(t)*y(t-LAG)'], LAG >=0, given
% the assumed VCVX and the corresponding VCVK (computed by vcv.m),
% for the model solution:
%
% y(t) = D*k(t) + F*x(t)
%
% k(t+1) = G*k(t) + H*x(t)
%
% We use the facts that x(t) is serially uncorrelated and that
% k(t), depending only on past x, is uncorrelated with the current
% and all future x. Then we find that:
%
%       If LAG = 0: ACM = D*VCVK*D' + F*VCVX*F'
%
%       If LAG > 0: ACM = D*G^LAG*VCVK*D' + D*G^(LAG-1)*H*VCVX*F'
%
% We also check for the magnitude of the imaginary component
% of the resulting ACM, and supress this component to issue
% report.
%
% Note that, for LAG = 0, ACM is simply the stationary variance-
% covariance matrix of y(t).
% ****************************************************************

% ****************************************************************
% [1] Compute the true ACM
% ****************************************************************

if LAG == 0
        ACM = D*VCVK*D' + F*VCVX*F';
        else
        ACM = D*G^LAG*VCVK*D' + D*G^(LAG-1)*H*VCVX*F';
        end

% ****************************************************************
% [2] Check for magnitude of imaginary component of ACM and
% suppress it.
% ****************************************************************

IMACM=imag(ACM);

if max(max(abs(IMACM))) > .0000001
        disp('Maximum imaginary component of D is larger than')
        disp('10 to the -8.  The magnitude of this part is:')
        disp(max(max(abs(IMACM))))
        pause
        end

ACM = real(ACM);

