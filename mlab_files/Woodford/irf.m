function Y = irf(SHOCK, NIR, D, F, G, H);
%
% ****************************************************************
%                    IMPULSE RESPONSE FUNCTIONS
% ****************************************************************
%
% ****************************************************************
% This function computes the path of y(t) in response to a unit
% impulse in the SHOCK-th component of x(0), with all other
% entries zero, k(0) = 0 as initial condition, and x(t) = 0 for
% t > 0, in:
%
% y(t) = D*k(t) + F*x(t)
%
% k(t+1) = G*k(t) + H*x(t)
%
% The paths are computed with length NIR, as a matrix NYxNIR.
% ****************************************************************

% ****************************************************************
% [1] Compute the true impulse responses
% ****************************************************************

NY = size(D, 1);
NX = size(F, 2);
NK = size(D, 2);

X = zeros(NX, 1);
K = zeros(NK, NIR);
Y = zeros(NY, NIR);

X(SHOCK) = 1;

K(:, 2) = G*K(:, 1) + H*X;

for J = 3:NIR

        K(:, J) = G*K(:, J-1);
        end

Y(:, 1) = D*K(:, 1) + F*X;

for J = 2:NIR

        Y(:, J) = D*K(:, J);
        end

% ****************************************************************
% [2] Check the magnitude of imaginary components of responses
% and suppress these components to issue report
% ****************************************************************

IMY=imag(Y);

if max(max(abs(IMY))) > .0000001
        disp('Maximum imaginary component is larger than')
        disp('10 to the -8.  The magnitude of this part is:')
        disp(max(max(abs(IMY))))
        pause
        end

Y=real(Y);

