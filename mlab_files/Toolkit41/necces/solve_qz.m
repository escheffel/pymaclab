% VERSION 3.0, September 2002, COPYRIGHT H. UHLIG.
% SOLVE_QZ.M is called by SOLVE.M, provided the option DO_QZ = 1 has been set.
% SOLVE_QZ.M employs the QZ-method, due to C.A. Sims (1989, 1996)
% and P. Klein (1997), adapted to the method of undetermined coefficients
% The QZ-method is perhaps numerically more stable than the generalized
% eigenvalue method employed in the standard solve algorithm.
% More importantly, the QZ-method works, even if PP is
% not diagonalizable.
% SOLVE_QZ.M in turn makes use of the routines
% QZDIV.M and QZSWITCH.M, written by C.A. Sims (1996).
%
% Once the PP matrix has been obtained, SOLVE_QZ proceeds
% to solve for the other matrices in the same way as SOLVE.M.

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.


% Theory: with the matrix definitions in the paper,
% find the QZ-decomposition for Delta and Xi, i.e.
% unitary matrices U and V and upper triangular matrices
% Delta_up and Xi_up so that
% U * Delta * V = Delta_up
% U * Xi * V = Xi_up
% and such that the ratios of the diagonal entries
% | Xi_up(i,i) / Delta_up(i,i) | are ordered ascendingly.
% Let U[a,b] be the m x m submatrices of U, where a,b = 1,2
% and likewise for V', where V' is the complex conjugate transpose.
% If V'[2,1] and U[2,1] are invertible
% then it can be shown that P = - inv(V'[2,1])*V'[2,2]
% is a solution to the matrix quadratic equation with
% the most stable roots selected.
% A problem that this version of the program does not
% yet solve is to drop a low real eigenvalue, if otherwise
% two conjugate complex roots are separated.

[Delta_up,Xi_up,UUU,VVV,Xi_eigvec]=qz(Delta_mat,Xi_mat);
Xi_eigval = diag( diag(Xi_up)./max(diag(Delta_up),TOL) );
[Xi_sortabs,Xi_sortindex] = sort(abs(diag(Xi_eigval)));
Xi_sortval = diag(Xi_eigval(Xi_sortindex,Xi_sortindex));
Xi_select = 1 : m_states;
stake = max(abs(Xi_sortval(Xi_select))) + TOL;
  % sorting the eigenvalues using the Sims (1996) code:
[Delta_up,Xi_up,UUU,VVV]=qzdiv(stake,Delta_up,Xi_up,UUU,VVV);


% Some warning messages:
if imag(Xi_sortval(m_states))~=0,
  if (abs( Xi_sortval(m_states) - conj(Xi_sortval(m_states+1)) ) < TOL),
  % NOTE: THIS LAST LINE MIGHT CREATE PROBLEMS, IF THIS EIGENVALUE OCCURS MORE THAN ONCE!!
  % IF YOU HAVE THAT PROBLEM, PLEASE TRY DO_QZ = 0.   
%  message = ['SOLVE_QZ.M: You may be in trouble. You have complex eigenvalues, and   '
%             '   I am not yet programmed to find a real eigenvalue to drop to only   '
%             '   have conjugate-complex pairs.  You should probably try DO_QZ = 0.   '
%             '   Put differently: now, your PP matrix will contain complex numbers.  '];
%  if DISPLAY_IMMEDIATELY, disp(message); end;
  warnings = [warnings;message];
  drop_index = 1;
  while (abs(imag(Xi_sortval(drop_index)))>TOL) & (drop_index < m_states),
     drop_index = drop_index + 1;
  end;
  if drop_index >= m_states,
     message = ['SOLVE_QZ.M: You are in trouble. You have complex eigenvalues. I cannot '
                '   find a real eigenvalue to drop to only have conjugate-complex pairs.'
                '   Put differently: your PP matrix will contain complex numbers. Sorry!'
                '   Try increasing the dimension of your state space. You may then get  '
                '   sunspots, too.                                                      '];
     if DISPLAY_IMMEDIATELY, disp(message); end;
     warnings = [warnings;message];
  else
     message = ['SOLVE_QZ.M: I will drop the lowest real eigenvalue to get real PP.     '
                '         I hope that is ok. You may have sunspots.                     ']; 
     if DISPLAY_IMMEDIATELY, disp(message); end;
     warnings = [warnings;message];
     for i = drop_index:m_states,
        [Delta_up,Xi_up,UUU,VVV]=qzswitch(i,Delta_up,Xi_up,UUU,VVV);
     end;
     % ... thus, successively moving the low real eigenvalue into the m+1-st 
     % position.  There may be a more elegant way, but this works.
     Xi_select = [ 1: (drop_index-1), (drop_index+1):(m_states+1)];
   end; % if drop_index >= m_states,
end;
if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
   message = ['SOLVE_QZ.M: You may be in trouble. There are stable roots NOT used for '
              '         PP. I have used the smallest roots: I hope that is ok.        '  
              '         Your alternative: DO_QZ = 0 and manual roots selection.       '
              '         For manual root selection, take a look at the file solve.m    '
              '         Watch out for sunspot solutions.                              '
              '         Better yet: move the time index of some endogenous variables  '
              '         back by one and turn them into (predetermined) state variables'];
   if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
   end; 
end;
if max(abs(Xi_sortval(Xi_select)))  > 1 + TOL,
   message = ['SOLVE_QZ.M: You may be in trouble.  There are unstable roots used for  '
              '         PP. Keep your fingers crossed or change your model.           '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;
if abs( max(abs(Xi_sortval(Xi_select))) - 1  ) < TOL,
   message = ['SOLVE_QZ.M: Your matrix PP contains a unit root. Your steady state is  '
              '         is probably not unique, is it?  Should not be a problem, but  '
              '         you do not have convergence back to steady state after a shock'
              '         and you should better not trust long simulations.             '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;
Lambda_mat = diag(Xi_sortval(Xi_select)); % to help sol_out.m

% Proceeding with the calculations:

VVV = VVV'; % complex conjugate transpose
VVV_2_1 = VVV(m_states +1 : 2*m_states, 1 : m_states);
VVV_2_2 = VVV(m_states +1 : 2*m_states, m_states+1 : 2*m_states);
UUU_2_1 = UUU(m_states+1 : 2*m_states,1:m_states);
VVV = VVV'; % back into original form, just in case.

if abs(det(UUU_2_1)) < TOL,
   message = ['SOLVE_QZ.M: A necessary condition for the validity of the formula for  '
              '         the matrix PP is not satisfied.  I will proceed anyhow, but   '
              '         keep your fingers crossed.                                    '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;
if abs(det(VVV_2_1)) < TOL,
   message = ['SOLVE_QZ.M: A matrix which I need to invert for calculating PP is not  '
              '         invertible.  I will proceed anyhow, but you are probably in   '
              '         trouble, and PP contains infinite entries. Try DO_QZ = 0.     '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;

PP = - VVV_2_1 \ VVV_2_2;
PP_imag = imag(PP);
PP = real(PP);

if sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001,
   message = ['SOLVE_QZ.M: PP is complex.  I proceed with the real part only.         '  
              '         Hope that is ok, but you are probably really in trouble!!     '
              '         You should better check everything carefully and be           '
              '         distrustful of all results which follow now.                  '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;

% End of calculating the PP matrix.  Now comes the rest.
calc_qrs;
