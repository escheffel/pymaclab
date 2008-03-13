% VERSION 2.1, APRIL 1997, COPYRIGHT H. UHLIG.
% SOLVE.M solves for the decision rules in a linear system,
% which is assumed to be of the form
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,
% where it is assumed that x(t) is the endogenous state vector,
% y(t) the other endogenous variables and z(t) the exogenous state
% vector.  It is assumed that the row dimension of AA is at least as large as
% the dimensionality of the endogenous state vector x(t).  
% The program solves for the equilibrium law of motion
% x(t) = PP x(t-1) + QQ z(t)
% y(t) = RR x(t-1) + SS z(t).
% To use this program, define the matrices AA, BB, .., NN.
% SOLVE.M then calculates PP.   
% 
% A few additional variables are used
% overwriting variables with the same names
% that might have been used before.  They are:
% sumimag, sumabs, message, warnings, CC_plus, CC_0, Psi_mat, Gamma_mat, Theta_mat, Xi_mat,
% Delta_mat, Xi_eigvec, Xi_eigval, Xi_sortabs, Xi_sortindex, Xi_sortvec, Xi_sortval,
% Xi_select, drop_index, Omega_mat, Lambda_mat, PP_imag, VV, LLNN_plus_MM, QQSS_vec
% 
% Source: H. Uhlig (1995) "A Toolkit for Solving Nonlinear Dynamic
% Stochastic Models Easily," Discussion Paper, Institute for
% Empirical Macroeconomis, Federal Reserve Bank of Minneapolis #101 or
% Tilburg University, CentER DP 9597.
%
% This update includes the suggestion by Andrew Atkeson to use generalized
% eigenvalues to perform the computations.  IN PARTICULAR, THE DEFINITION OF
% PSI_MAT, GAMMA_MAT and THETA_MAT have been changed!
% Furthermore, there is an option to use the QZ-method, due to Sims (1989, 1996)
% instead.  You can activate it by setting DO_QZ = 1;
%
% You can also select roots manually.  For the manual selection procedure,
% see the instructions upon inspecting this file (filename: solve.m)

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %  MANUAL SELECTION OF ROOTS PROCEDURE.  INSTRUCTIONS:     %
      % For manual selection, set MANUAL_ROOTS = 1 and %
      % define Xi_manual somewhere earlier in your calculations. %
      % Xi_manual should be a vector of length m_states with     %
      % distinct integer entries between 1 and 2*m_states. The   %
      % program then uses the roots Xi_sortval(Xi_manual) and    %
      % the corresponding eigenvectors Xi_sortvec(Xi_manual).    %
      % Thus, to choose the desired roots, run the program once  %
      % with automatic root selection, take a look at Xi_sortval,%
      % and Xi_sortvec and write down the indices of the desired %  
      % roots.                                                   %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[q_expectational_equ,m_states] = size(FF);
[l_equ,n_endog ] = size(CC);
% [l_equ,m_states] = size(AA);
% [l_equ,n_endog ] = size(CC);
k_exog = min(size(NN));
sumimag = sum(sum(abs(imag(AA))))+sum(sum(abs(imag(BB))))+sum(sum(abs(imag(CC))));
sumabs  = sum(sum(abs(AA)))      +sum(sum(abs(BB)))      +sum(sum(abs(CC)));
if sumimag / sumabs > .000001,
  message = ['SOLVE.M: I hate to point this out to you, but some of your matrices    '  
             '         contain complex numbers, which does not make much sense. You  '
             '         should check your steady state parameters and calculations.   '
             '         I will proceed anyhow, but you will probably get nonsense.    '];
  if DISPLAY_IMMEDIATELY, disp(message); end;
  warnings = [warnings;message];
end;
if rank(CC)<n_endog,
  error('SOLVE.M: Sorry!  Rank(CC) needs to be at least n! Cannot solve for PP. ')
  message = 'SOLVE.M: Sorry!  Rank(CC) needs to be at least n! Cannot solve for PP. ';
  if DISPLAY_IMMEDIATELY, disp(message); end;
  warnings = [warnings;message];
else
    if l_equ == 0,
        CC_plus = pinv(CC);
        CC_0 = (null(CC'))';
        Psi_mat   = FF ;
        Gamma_mat = - GG;
        Theta_mat =  - HH;
        Xi_mat    = [ Gamma_mat,     Theta_mat
            eye(m_states), zeros(m_states) ];
        Delta_mat = [ Psi_mat,       zeros(m_states)
            zeros(m_states), eye(m_states) ];
    else
        CC_plus = pinv(CC);
        CC_0 = (null(CC'))';
        Psi_mat   = [ zeros(l_equ-n_endog,m_states)
            FF - JJ*CC_plus*AA           ];
        Gamma_mat = [ CC_0 * AA
            JJ*CC_plus*BB - GG + KK*CC_plus*AA ];
        Theta_mat = [ CC_0 * BB
            KK*CC_plus*BB - HH                 ];
        Xi_mat    = [ Gamma_mat,     Theta_mat
            eye(m_states), zeros(m_states) ];
        Delta_mat = [ Psi_mat,       zeros(m_states)
            zeros(m_states), eye(m_states) ];
 end;
  if DO_QZ,
     solve_qz;
  else
     [Xi_eigvec,Xi_eigval] = eig(Xi_mat,Delta_mat);
     if rank(Xi_eigvec)<m_states,
        error('SOLVE.M: Sorry! Xi is not diagonalizable! Cannot solve for PP')
        message = ['SOLVE.M: Sorry! Xi is not diagonalizable! Cannot solve for PP.         '
                   '         Try to run your program again with DO_QZ = 1.                 '];
        if DISPLAY_IMMEDIATELY, disp(message); end;
        warnings = [warnings;message];
     else
       [Xi_sortabs,Xi_sortindex] = sort(abs(diag(Xi_eigval)));
       Xi_sortvec = Xi_eigvec(1:2*m_states,Xi_sortindex);
       Xi_sortval = diag(Xi_eigval(Xi_sortindex,Xi_sortindex));
       Xi_select = 1 : m_states;
       if imag(Xi_sortval(m_states))~=0,
         if (abs( Xi_sortval(m_states) - conj(Xi_sortval(m_states+1)) ) < TOL),
         % NOTE: THIS LAST LINE MIGHT CREATE PROBLEMS, IF THIS EIGENVALUE OCCURS MORE THAN ONCE!!
         % IF YOU HAVE THAT PROBLEM, PLEASE TRY MANUAL ROOT SELECTION.  
           drop_index = 1;
           while (abs(imag(Xi_sortval(drop_index)))>TOL) & (drop_index < m_states),
             drop_index = drop_index + 1;
           end;
           if drop_index >= m_states,
             message = ['SOLVE.M: You are in trouble. You have complex eigenvalues, and I cannot'
                        '   find a real eigenvalue to drop to only have conjugate-complex pairs.'
                        '   Put differently: your PP matrix will contain complex numbers. Sorry!'
                        '   Try increasing the dimension of your state space. You may then get  '
                        '   sunspots, too.                                                      '];
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
           else
             message = ['SOLVE.M: I will drop the lowest real eigenvalue to get real PP.        '
                        '         I hope that is ok. You may have sunspots.                     ']; 
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
             Xi_select = [ 1: (drop_index-1), (drop_index+1):(m_states+1)];
           end; % if drop_index >= m_states,
         end; % if (abs( Xi_sortval(m_states) - ...
       end; % if imag(Xi_sortval(m_states))~=0,
       if MANUAL_ROOTS,
         message = ['SOLVE.M: You have chosen to select roots manually.  I am crossing my   '
                    '         fingers that you are doing it correctly.  In particular,      '
                    '         you should have defined Xi_manual.  Type help solve           '
                    '         and inspect SOLVE.M to get further information on how to do it'];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
         if exist('Xi_manual'),
            Xi_select = Xi_manual;
         else
            message = ['SOLVE.M: You have not defined Xi_manual.  Either define it or turn off '
                       '         the manual roots selection procedure with                     '
                       '         MANUAL_ROOTS = 0                                              '
                       '         Right now, I better let your calculations crash - sorry!      '
                       '         If you get results, they are based on previous calculations.  '];
            disp(message);
            warnings = [warnings;message];
         end; % if exist('Xi_manual'),
       else
         if max(Xi_select) < 2*m_states,
           if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
             message = ['SOLVE.M: You may be in trouble. There are stable roots NOT used for PP.'
                        '         I have used the smallest roots: I hope that is ok.            '  
                        '         If not, try manually selecting your favourite roots.          '
                        '         For manual root selection, take a look at the file solve.m    '
                        '         Watch out for sunspot solutions.                              '
                        '         Better yet: move the time index of some endogenous variables  '
                        '         back by one and turn them into (predetermined) state variables'];
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
           end; % if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
         end; % if max(Xi_select) < 2*m_states,
       end; % if MANUAL_ROOTS,
       if max(abs(Xi_sortval(Xi_select)))  > 1 + TOL,
         message = ['SOLVE.M: You may be in trouble.  There are unstable roots used for PP. '
                    '         Keep your fingers crossed or change your model.               '];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       end; % if max(abs(Xi_sortval(Xi_select))) ... 
       if abs( max(abs(Xi_sortval(Xi_select))) - 1  ) < TOL,
         message = ['SOLVE.M: Your matrix PP contains a unit root. You probably do not have '
                    '         a unique steady state, do you?  Should not be a problem, but  '
                    '         you do not have convergence back to steady state after a shock'
                    '         and you should better not trust long simulations.             '];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       end; % if abs( max(abs(Xi_sortval(Xi_select))) - 1 ... 
       Lambda_mat = diag(Xi_sortval(Xi_select));
       Omega_mat  = [Xi_sortvec((m_states+1):(2*m_states),Xi_select)];
       if rank(Omega_mat)<m_states,
         message = 'SOLVE.M: Sorry! Omega is not invertible. Cannot solve for PP.          ';
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       else
         PP = Omega_mat*Lambda_mat/Omega_mat;
         PP_imag = imag(PP);
         PP = real(PP);
         if sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001,
           message = ['SOLVE.M: PP is complex.  I proceed with the real part only.            '  
                      '         Hope that is ok, but you are probably really in trouble!!     '
                      '         You should better check everything carefully and be           '
                      '         distrustful of all results which follow now.                  '];
           if DISPLAY_IMMEDIATELY, disp(message); end;
           warnings = [warnings;message];
         end; % if sum(sum(abs(PP_imag)))
      end; % if rank(Omega_mat)<m_states,
      % End of calculating the PP matrix.  Now comes the rest.
      calc_qrs;
    end; % if rank(Xi_eigvec)<m_states,
  end; % if DO_QZ,
end; % if rank(CC)<n_endog,


       
     
      

    
  
      
      

     
