% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% SOL_OUT.M prints the coefficients of the decision rules,
% delivered by SOLVE.M.
% It is assumed, that VARNAMES, a matrix with m+n+k rows has
% been set, containing the names of all the variables.
% This program overwrites m_states, k_exog and n_endog.


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.


[m_states,k_exog] = size(QQ);
[n_endog,k_exog] = size(SS);
disp('Exogenous states z(t):');
disp(VARNAMES((m_states+n_endog+1):(m_states+n_endog+k_exog),:));
disp(' ');
disp('Endogenous states x(t):');
disp(VARNAMES(1:m_states,:));
disp(' ');
if DISPLAY_ROOTS,
   disp('All the roots are:');
   disp('    root         abs(root)   ');
   disp([diag(Xi_eigval(Xi_sortindex,Xi_sortindex)),...
   abs(diag(Xi_eigval(Xi_sortindex,Xi_sortindex)))] );
   disp('The chosen roots are:');
   disp('    root         abs(root)   ');
   disp([diag(Lambda_mat),abs(diag(Lambda_mat))]);
   disp(' ');
end;
disp('PP: Recursive equilibrium law of motion for x(t) on x(t-1):');
disp(PP);
disp('QQ: Recursive equilibrium law of motion for x(t) on z(t):');
disp(QQ);
disp(' ');
disp('Other endogenous variables y(t):');
disp(VARNAMES((m_states+1):(m_states+n_endog),:));
disp(' ');
disp('RR: Recursive equilibrium law of motion for y(t) on x(t-1):');
disp(RR);
disp('SS: Recursive equilibrium law of motion for y(t) on z(t):');
disp(SS);
disp(' ');
