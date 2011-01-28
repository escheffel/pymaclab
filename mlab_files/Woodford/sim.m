function Y = sim(K0, X, D, F, G, H);
%
% ****************************************************************
%                    SIMULATION OF SOLUTION
% ****************************************************************
%
% ****************************************************************
% This function simulates the path of {y(t)} = Y given an initial
% condition k(0) = K0 and a sequence of shocks {x(t)} = X, from:
%
% y(t) = D*k(t) + F*x(t)
%
% k(t+1) = G*k(t) + H*x(t)
%
% Note that X and Y are written as matrices with NX and NY rows,
% respectively, and as many columns as the length of the desired
% simulation, NSIM.
% ****************************************************************

% ****************************************************************
% [1] Compute the true simulated series
% ****************************************************************

NY = size(D, 1);
NX = size(F, 2);
NK = size(D, 2);
NSIM = size(X, 2);

K = zeros(NK, NSIM);
Y = zeros(NY, NSIM);

K(:, 1) = K0;

for J = 2:NSIM

        K(:, J) = G*K(:, J-1) + H*X(:, J-1);
        end

for J = 1:NSIM

        Y(:, J) = D*K(:, J) + F*X(:, J);
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

