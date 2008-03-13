% VERSION 4.0, November 2002, COPYRIGHT H. UHLIG.
% CALC_QRS.M calculates the matrices Q, R and S as well
% as some other matrices given the matrix P. I.e.W with the property [x(t)',y(t)',z(t)']=W [x(t)',z(t)'].
% (The program uses the name PP for the matrix P, etc.)

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

if l_equ == 0,
	RR = zeros(0,m_states);
	VV = [ kron(NN',FF)+kron(eye(k_exog),(FF*PP+GG)), ...
                             kron(NN',JJ)+kron(eye(k_exog),KK) ];
else
	RR = - CC_plus*(AA*PP+BB);
	VV = [ kron(eye(k_exog),AA),   kron(eye(k_exog),CC)
           kron(NN',FF)+kron(eye(k_exog),(FF*PP+JJ*RR+GG)), ...
                             kron(NN',JJ)+kron(eye(k_exog),KK) ];
end;
if ( (rank(VV) < k_exog*(m_states+n_endog)) & ...
                                (~IGNORE_VV_SING) ),
   message = ['SOLVE.M: Sorry! V is not invertible.  Cannot solve for QQ and SS. You  '
              '         can try setting IGNORE_VV_SING = 1 and wish for the best...   '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
else
   if ( (rank(VV) < k_exog*(m_states+n_endog)) ),
     message = ['SOLVE.M: Warning! V is not invertible.  However, you have set          '
                '         IGNORE_VV_SING = 1, and thus, since you have told me to       '
                '         ignore this, I will proceed.  Keep your fingers crossed...    ']
     if DISPLAY_IMMEDIATELY, disp(message); end;
     warnings = [warnings;message];
   end;
   LLNN_plus_MM = LL*NN + MM;
   QQSS_vec = - VV \ [ DD(:)
                       LLNN_plus_MM(:) ];
   if max(abs(QQSS_vec)) == Inf,
      message = ['SOLVE.M: You probably are in trouble!  QQ or SS contain undefined      '
                 '         entries! Most likely, the matrix VV is not invertible.        '];
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
   end;           
   QQ = reshape(QQSS_vec(1:m_states*k_exog),m_states,k_exog);
   SS = reshape(QQSS_vec((m_states*k_exog+1):((m_states+n_endog)*k_exog)),n_endog,k_exog);
   WW = [ eye(m_states)         , zeros(m_states,k_exog)
          RR*pinv(PP)           , (SS-RR*pinv(PP)*QQ) 
          zeros(k_exog,m_states), eye(k_exog)            ];
end;
