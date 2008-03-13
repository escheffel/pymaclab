function [D,F,G,H] = soldsf(Br,Cr,Lr,NY,NX,NK,NF)
% soldsf.m
%
% Syntax:  [D,F,G,H] = soldsf(Br,Cr,Lr,NY,NX,NK,NF)
%
% ****************************************************************
%                  SOLUTION OF DYNAMIC SYSTEMS
% ****************************************************************
%
% ****************************************************************
% The starting point now is the output of redsf.m, i.e., the
% already reduced system. This system is represented by the
% current Br and Cr (recall that A must have been reduced to zeros
% all over, except for a southeastern identity block). We have
% also kept track of Lr, the reordering matrix.
%
% The dynamic solution is found for the "d(t)" equations of the
% reduced system, which are:
%
% Ed(t+1|t) = Bse*d(t) + Cs*x(t)
%
% where Bse is the southeastern block of Br.
%
% Once that solution is found, we obtain the complete solution
% including the flow part of the system.
%
% The final solution will be cast in the form:
%
% y(t) = D*k(t) + F*x(t)
%
% k(t+1) = G*k(t) + H*x(t)
%
% where k(t) is the vector of predetermined variables.
% ****************************************************************

% ****************************************************************
% [1] Obtain eigenvalues of Bse and associated left eigenvectors
%
% Denote by LE the matrix of left eigenvectors of Bse, i.e., if
% MU is the diagonal matrix of eigenvalues, then LE*Bse = MU*LE.
% ****************************************************************

[LE1 MU1] = eig(Br(NF+1:NY, NF+1:NY)');
LE1 = LE1';
MU1 = MU1';

[LE2 MU2] = eig(Br(NF+1:NY, NF+1:NY)', 'nobalance');
LE2 = LE2';
MU2 = MU2';

if sum(sum(abs(LE1*Br(NF+1:NY, NF+1:NY) - MU1*LE1))) > ...
                sum(sum(abs(LE2*Br(NF+1:NY, NF+1:NY) - MU2*LE2)))
        LE = LE2;
        MU = MU2;
        else
        LE = LE1;
        MU = MU1;
        end

% ****************************************************************
% [2] Check standard conditions for solvability
%
% First, we locate and count the unstable eigenvalues - those with
% absolute value greater than 1 plus a tolerance level (set to
% 0.000001).
%
% If there are unstable eigenvalues at all, we then turn LE and MU
% into the submatrices of the original LE and MU corresponding to
% these unstable eigenvalues. We take the first columns of LE (as
% many as it has rows, i.e., as many as the number of unstable
% eigenvalues) to form a square matrix with which the Boyd-Dotsey
% rank condition is checked.
%
% Otherwise, we set LE and MU to [0].
%
% In any case, we also check whether the number of unstable
% eigenvalues equals NY-NF-NK, i.e., the number of non-
% predetermined variables left in d(t).
% ****************************************************************

LUE = abs(diag(MU)) > 1.000001;

NUE = sum(LUE);

if NUE > 0

        LE = LE(LUE, :);
        MU = MU(LUE, LUE);
        
        if NUE > rank(LE(:, 1:NUE), 10^(-5))
                error('Rank condition violated')
                end
        
        else
        
        LE = [0];
        MU = [0];
                
        end

if NUE > NY-NF-NK

        error('Too many unstable eigenvalues')
        end

if NUE < NY-NF-NK

        error('Too few unstable eigenvalues')
        end

% ****************************************************************
% [3] Solution in the absence of unstable roots
%
% Without unstable eigenvalues, the solution strategy can be very
% simple because there are no variables to be solved "forward".
% In this case, d(t) and k(t) are the same, Bse is simply the
% feedback matrix and C governs the response to x(t).
%
% So, the NY-NK equations are for flows, having the form:
%
% f(t) + Bne*k(t) + Cn*x(t) = 0
%
% and the last NK are in the form:
%
% k(t+1) = Bse*k(t) + Cs*x(t)
%
% because Ek(t+1|t)=k(t+1) by predeterminateness.
%
% We have nothing further ado in this case.
% ***************************************************************

% ***************************************************************
% [5] Solution in the presence of unstable eigenvalues
%
% Consider again the d(t) equations:
%
% Ed(t+1|t) = Bse*d(t) + Cs*x(t)
%
% Then, premultiplying by LE:
%
% LE*Ed(t+1|t) = LE*Bse*d(t) + LE*Cs*x(t)
%
% But LE*Bse = MU*LE by definition (this holds for any subset of
% eigenvectors), and we can rewrite the latter equation as:
%
% Eu(t+1|t) =  MU*u(t) + RHO*z(t),
%
% where u(t) = LE*d(t) and RHO = LE*Cs. But MU is diagonal, and
% so we have arrived at "separated" dynamic equations for u(t).
%
% We can iterate each of these NUE equations "forward" to obtain,
% imposing stability in every case:
%
% u(i)(t) = -[RHO(i, :)/(MU(i, i)*I)]*x(t)
%
% We denote:
%
% THETA(i, :) = -RHO(i, :)/(MU(i, i)*I)
%
% and form the matrix THETA by stacking all the NUE rows
% THETA(i, :), each of length cols(Cr) = cols(Cs).
%
% But now we have:
%
% THETA*x(t) = u(t) = LE*d(t)
%
% which we want to solve for d(t) in terms of x(t).
%
% Let h(t) be the first NUE of elements of d(t), and k(t) be
% the remaining NY-NF-NUE. Partition LE conformably into
% [LEH LEK]. We can then solve the above equality for h(t):
%
% h(t) = HK*k(t) + HX*x(t)
%
% where HK = -LEH\LEK and HX = LEH\THETA.
%
% We want to impose this solution on the dynamic system, i.e.,
% "substitute out" for h(t) using the above rule.
%
% The system being solved included equations of the form:
%
% B2*h(t) + B3*k(t) + Cr*x(t) = 0.
%
% Premultiply the previous equation by B2:
%
% -B2*h(t) + B2*HK*k(t) + B2*HX*x(t) = 0,
%
% and add to the latter to get:
%
% (B3 + B2*HK)*k(t) + (Cr + B2*HX)*x(t) = 0.
% ****************************************************************

if NUE > 0

        RHO = LE*Cr(NF+1:NY, :);

        THETA = zeros(NUE, NX);

        for J = 1:NUE;
                THETA(J, :) = ...
                        -RHO(J, :)/(MU(J,J)*eye(NX));
                end

        B2 = Br(:, NF+1:NF+NUE);
        B3 = Br(:, NF+NUE+1:NY);

        if NK > 0
                HK = - LE(:, 1:NUE)\LE(:, NUE+1:NY-NF);
                Br(:, NF+NUE+1:NY) = B3+B2*HK;
                end

        HX = LE(:, 1:NUE)\THETA;
        Cr = Cr+B2*HX;

        Br(:, NF+1:NF+NUE) = zeros(size(B2));
        Br(NF+1:NF+NUE, NF+1:NF+NUE) = eye(NUE);

        if NK > 0
                Br(NF+1:NF+NUE, NF+NUE+1:NY) = -HK;
                end

        Cr(NF+1:NF+NUE, :) = -HX;

        end

% ****************************************************************
% [6] The complete solution
%
% We obtain the D and F matrices for the new ordering and then
% transform them in D and F matrices for the old ordering. We also
% obtain G and H.
% ****************************************************************

D = [-Br(1:NY-NK, NY-NK+1:NY)
             eye(NK)         ];

F = [-Cr(1:NY-NK, :)
     zeros(NK, NX)];

D = Lr'*D;
F = Lr'*F;

G = Br(NY-NK+1:NY, NY-NK+1:NY);
H = Cr(NY-NK+1:NY, :);


