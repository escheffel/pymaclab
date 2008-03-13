function ACF = acf(NACF, VCVX, VCVK, D, F, G, H);
%
% ****************************************************************
%                 AUTOCOVARIANCE FUNCTIONS OF Y
% ****************************************************************
%
% ****************************************************************
% This function computes the NACF first autocovariances of each
% entry in y(t) (not including the variance itself), given the
% assumed VCVX and the corresponding VCVK computed by vcv.m, for
% the model solution:
%
% y(t) = D*k(t) + F*x(t)
%
% k(t+1) = G*k(t) + H*x(t)
%
% As argued in the notes to acm.m, we know that:
%
% E[y(t)*y(t-J)'] = D*G^J*VCVK*D' + D*G^(J-1)*H*VCVX*F'
%
% for any J > 0. We compute such autocovariance matrix for
% J = 1, ..., NACF, and for each J we take its main diagonal as
% the J-th autocovariances of each entry in y(t). The output is a
% NYxNACF matrix, the rows of which contain the autocovariance
% functions for the corresponding entries in y(t), for lags 1 to
% NACF.
%
% We again check for the magnitude of the imaginary component
% of the resulting ACF, and supress this component to issue
% report.
% ****************************************************************

% ****************************************************************
% [1] Compute the true ACF
% ****************************************************************

NY = size(D, 1);

ACF = zeros(NY, NACF);

GVCVK = VCVK*D';
GVCVX = H*VCVX*F';

for J = 1:NACF
        GVCVK = G*GVCVK;
        ACF(:, J) = diag(D*(GVCVK + GVCVX));
        GVCVX = G*GVCVX;
        end

% ****************************************************************
% [2] Check for magnitude of imaginary component of ACF and
% suppress it.
% ****************************************************************

IMACF=imag(ACF);

if max(max(abs(IMACF))) > .0000001
        disp('Maximum imaginary component of D is larger than')
        disp('10 to the -8.  The magnitude of this part is:')
        disp(max(max(abs(IMACF))))
        pause
        end

ACF = real(ACF);

