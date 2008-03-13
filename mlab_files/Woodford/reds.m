% reds.m
%
% ****************************************************************
%                    REDUCTION OF DYNAMIC SYSTEMS
% ****************************************************************
%
% ****************************************************************
% Consider a singular dynamic system of the form:
%
% A*Ey(t+1|t) = B*y(t) + C*x(t),
%
% where y(t) is ordered so that the non-predetermined variables
% appear first and the predetermined variables appear last, and
% x(t) is a martingale difference sequence.
%
% This program reduces it to a dynamic system of the form:
%
% | 0  0 | |Ef(t+1|t)|   | I   Bne | |f(t)|    | Cn |
% |      | |         | = |         | |    |  + |    | x(t)
% | 0  I | |Ed(t+1|t)|   | 0   Bse | |d(t)|    | Cs |
%
% where:
%
%            | f(t) |
% L*y(t) =   |      |
%            | d(t) |
%
% with L being a permutation matrix.
%
% Throughout, we denote NY = dim(y), NX = dim(x), NK = dim(k)
% (where k(t) is the vector of predetermined variables) and
% NF = dim(f).
%
% Matrices A, B and C and parameters NY, NX, and NK are the
% inputs.
%
% As we progress with the reduction process, we will have at each
% iteration systems of the form:
%
% | 0   0 | |Ef(t+1|t)|   | I     Bne| |f(t)|    | Cn |
% |       | |         | = |          | |    |  + |    | x(t)
% | 0  Ase| |Ed(t+1|t)|   | 0     Bse| |d(t)|    | Cs |
%
% At each step we "remove a flow" (i.e., transfer a variable from
% d to f). We iterate until either det(Ase) becomes sufficiently
% large (in absolute value) or no additional flows are removed at
% the prior iteration.
%
% Each iteration contains three consecutive blocks: singular
% value decomposition, QR factorization and classical state
% reduction.
% ****************************************************************

% ****************************************************************
% [1] Preliminary check of solvability
%
% We first check whether the system satisfies a necessary
% condition for solvability: |Az-B| nonzero.
% ****************************************************************

J=1;

for Z = -100:1/4:100
        POLYN(J) = det(A*Z-B);
        J = J+1;
        end

if POLYN == 0
        error('Model unsolvable: |Az-B| identically null')
        end

% ****************************************************************
% [2] INITIALIZATION
%
% Initialize L = I and NF = 0.
%
% Define:
%
% FRST = 1 <=> it is the first iteration
% SING = 1 <=> det(Ase) is close enough to zero
%
% NNF will be defined below as the number of new flows "removed"
% in the previous iteration. We initialize with NNF = 1 just to
% let the first iteration occur.
% ****************************************************************

L = eye(NY);

NF = 0;
NNF = 1;

FRST = 1;
SING = (abs(det(A)) < .000001);

while ((SING == 1)&(NNF > 0))
        
% ****************************************************************
% [3] THE SINGULAR VALUE DECOMPOSITION BLOCK
%
% The first problem is to get some more zero rows in Ase, turning
% it into:
%
% Ase= |  0     0  |
%      | Ase1  Ase2|
%
% Our standard way of doing this is to compute the singular value
% decomposition of Ase.  But this may not be necessary if we are
% on the first iteration and there are lots of rows of zeros in
% A: it will then suffice to interchange rows of A, placing zeros
% first, and perform the same reordering of B and C.
%
% If FRST=0 (it is not the first iteration) or if it is the
% first iteration and but NZEROS=0 (no rows of zeros found), then
% we perform the singular value decomposition of Ase = U*S*V',
% where U and V are unitary matrices (U*U'=I and V*V'=I) and S is:
%
% S = | Snw   0 |
%     |  0    0 |
%
% with Snw being a diagonal matrix with rank(Ase) positive
% numbers (the "singular values" of Ase) on its diagonal.
%
% (The rank of a matrix is theoretically the number of non-zero
% singular values.  We approximate this here as in the MATLAB
% rank command, where the computation is based on the number of
% singular values greater than a given tolerance level.)
%
% So, we can rewrite our d(t) equations as:
%
% S*V'*Ed(t+1|t) = U'*Bse*d(t) + U'*Cs*x(t)
%
% while leaving the f(t) equations unaltered.
%
% We also want to rearrange the ordering of the equations of the
% system so that those with zeros appear before the non-
% degenerate ones. This involves multiplying by a matrix:
%
%   | I  0  0|
%   | 0  0  I|
%   | 0  I  0|
%
% with suitably chosen dimensions of the identity matrices. Note
% that the upper I keeps the f(t) equations unchanged, and the
% and the southeastern 2x2 block turns S upside down.
% ****************************************************************

        if FRST == 1

                [SORT1 SORT2] = sort(max(abs(A')));

                A = A(SORT2, :);
                B = B(SORT2, :);
                C = C(SORT2, :);

                NZEROS = sum(SORT1 < .000001);

                RA = NY-NZEROS;

                FRST = 0;

                end

        if NZEROS+FRST == 0

                [U S V] = svd(A(NF+1:NY, NF+1:NY));

                RA = diag(S);
                RA = sum(RA > .00001);

                TR = eye(NY);
                TR(NF+1:NY, NF+1:NY) = U';
                TR2 = eye(NY);
                TR2 = [TR2(:, 1:NF) TR2(:, NY-RA+1:NY) ...
                                         TR2(:, NF+1:NY-RA)];
                TR = TR2*TR;

                A = TR*A;
                B = TR*B;
                C = TR*C;

                end

% ****************************************************************
% [4] THE QR FACTORIZATION BLOCK
%
% The "d(t) equations" of the system are now in the form:
%
% | 0     0  | |Eh(t+1|t)|
% |          | |         | = ...
% |Ase1  Ase2| |Ek(t+1|t)|
%
%                  |Bse11   Bse12| |h(t)|   | Cs1 |
%            ... = |             | |    | + |     | x(t)
%                  |Bse21   Bse22| |k(t)|   | Cs2 |
%
%
% where h(t) are the non-predetermined variables still in d(t)
% and k(t) are the predetermined ones.
%
% Notice that the first row of the system above:
%
% Bse11*h(t) + Bse12*k(t) + Cs1*x(t) = 0
%
% suggests that there are "candidate flows" in h(t), i.e.,
% elements of h(t) attached to behavioral equations that contain
% no leads. We want to solve for as many of these as possible.
%
% Our Bse11 above is a square matrix of size NY-NK-NF. In order
% not to try to solve for more candidate flows than possible, we
% will QR factor only its upper block of NY-max(RA,NK)-NF rows.
%
% If that block is [ ], QR factorization involves finding matrices
% P, Q and R such that [ ]*P = Q*R, P is a permutation matrix,
% Q is unitary and R is upper triangular. If the row rank of
% [ ] is less than full, then the last rows of R are zero.
%
% We then use P to reorder the system placing the candidate flows
% first within d(t). Next, we multiply the equations of these
% candidate flows by Q'. Note that the former operation does not
% affect C, while the latter does not affect L.
%
% We determine how many linearly independent rows of the equations
% there are by finding the total number of elements in each row
% that are close to zero; computing the number of rows with all
% zeros and subtracting from the number of rows.
% ****************************************************************

        [Q,R,P] = qr(B(NF+1:NY-max([RA NK]), NF+1:NY-NK));

        TR = eye(NY);

        TR(NF+1:NY-NK, NF+1:NY-NK) = P;

        A = A*TR;
        B = B*TR;
        L = TR'*L;

        TR = eye(NY);

        TR(NF+1:NY-max([RA NK]), NF+1:NY-max([RA NK])) = Q';

        A = TR*A;
        B = TR*B;
        C = TR*C;

        NNF = (abs(R) < .00001);
        NNF = sum(NNF');
        NNF = NNF/size(R, 2);
        NNF = sum(NNF == 1);
        NNF = size(R, 1)-NNF;

        if (NNF > 0)

                if (NNF < NY-max([RA NK])-NF)

                        if (max(max(abs(C(NF+NNF+1: ...
                            NY-max([RA NK]), :)))) > .000000001)
                                error('Failed to reduce')
                                end
                        
                        if (max(max(abs(B(NY-NK+1:NY, :)))) ...
                            > .000000001)
                                error('Failed to reduce')
                                end

                        end

                NF = NF+NNF;

% ****************************************************************
% [5] THE CLASSICAL STATE REDUCTION BLOCK:
%
% Then we undertake the classical state reduction problem for
% the dynamic system in the form:
%
% A*Ey(t+1|t) = B*y(t) + C*x(t)
%
% already with the posited structure
%
% | 0     0 | |Ef(t+1|t)|   |Bff   Bfd| |f(t)|    | Cf |
% |         | |         | = |         | |    |  + |    | x(t)
% |Adf   Add| |Ed(t+1|t)|   |Bdf   Bdd| |d(t)|    | Cd |
%
% where Bff is a nonsingular (NF x NF) matrix; f(t) are the first
% NF elements of y(t) and d(t) are the last NY-NF elements of
% y(t).
%
% We compute the state reduction using a transformation
% T(F) = T1*F+T2 such that the transformed system is in the form:
%
% |0   0 | |Ef(t+1|t)|   |I    Bne| |f(t)|    | Cn |
% |      | |         | = |        | |    |  + |    | x(t)
% |0  Ase| |Ed(t+1|t)|   |0    Bse| |d(t)|    | Cs |
%
% That is, the transformed system has specified control variables
% and state variables.
%
% The transformation matrices are:
%
% T1 = | 0      0 |         T2 = | I/Bff     0  |
%      |Adf/Bff 0 |              |-Bdf/Bff   I  |
%
% ****************************************************************


                BFF = B(1:NF, 1:NF);

                T1 = zeros(NY, NY);
                T1(NF+1:NY, 1:NF) = A(NF+1:NY, 1:NF)/BFF;

                T2 = eye(NY, NY);
                T2(1:NF, 1:NF) = eye(NF, NF)/BFF;
                T2(NF+1:NY, 1:NF) = -B(NF+1:NY, 1:NF)/BFF;

                A = T2*A-T1*B;
                B = T2*B;
                C = T2*C;
                                
                NZEROS = 0;

                end

        SING = (abs(det(A(NF+1:NY, NF+1:NY))) < .000001);

        end

% ****************************************************************
% [4] FINAL OUTPUT:
%
% If we don't find an irreducible singularity (i.e., iterations
% end and Ase is still singular), we put the system in "standard
% form".
% ****************************************************************

if (SING == 1)

        error('System contains irreducible singularity in A')

        else
        TR = eye(NY);
        TR(NF+1:NY, NF+1:NY) = ...
                A(NF+1:NY, NF+1:NY)\TR(NF+1:NY, NF+1:NY);
        A = TR*A;
        B = TR*B;
        C = TR*C;

        end



