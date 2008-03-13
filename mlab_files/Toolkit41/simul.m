% VERSION 4.1, May 2003, COPYRIGHT H. UHLIG.
%   minor correction in August 2002, correcting a MATLAB grammar error,
%   for simulations with given epsilons
%   May 2003: major correction: autocorrelation calculation interchanged lags and leads before,
%      this has now been corrected.
%      Thanks goes to Mathias Trabandt for pointing out the error.
%      
%
% SIMUL.M simulates the model.  It is controlled by several options to be described
% below.  All these options are set in OPTIONS.M
% It is assumed,
% that SOLVE.M has been executed before, so that the matrices
% NN, PP, QQ, RR and SS are available, describing the law of motion
%   x(t) = PP x(t-1) + QQ z(t)
%   y(t) = RR x(t-1) + SS z(t)
%   z(t) = NN z(t-1) + epsilon(t)
%
% SIMUL.M can be used in two modes, set by SIM_MODE.  
%
% In mode 1, one  time series only of length SIM_LENGTH is calculated.  To do so, one
% can either set SIM_RANDOM_START = 0 and 
% start the initial state variables at SIM_X_START for time t=0 and SIM_Z_START for t=1,
% with the simulations starting at date t=1,
% or set SIM_RANDOM_START = 1 and simulate first SIM_DISCARD 
% observations which are discarded to make
% sure that the series starts from the ergodic distribution.
% The simulated series are in sim_xyz, which is of size (m_states+n_endog+k_exog) x SIM_LENGTH
%    In this mode, you can also choose the draws for the shocks 
% beforehand  rather than generating them randomly:
% this is useful for date-by-date comparison with data.  To do so,
% set SIM_GIVEN_EPS = 1.  You also need to declare the matrix of shocks, called
%    given_eps 
% which needs to be a matrix of size k_exog x SIM_LENGTH.
%    given_eps(:,1) will not be used (as it is assumed to be "part of" SIM_Z_START),
%    given_eps(:,2) will be used for calculating z(2) = NN * SIM_Z_START + given_eps(:,2), etc.
% It is also possible to choose the initially discarded epsilons, although there are probably
% only very few uses for that. 
%
% In mode 2, SIM_N_SERIES time
% series are created and their autocorrelation tables complete with standard
% errors are calculated.  Again, the series are of length SIM_LENGTH after discarding
% SIM_START initial situations.  These calculations are an alternative to the calculations
% in MOMENTS.M.  In that case, you need to set SIM_N_LEAD_LAGS to some number, say, 6.
%
% Set DO_HP_FILTER to unity, if you want to apply the HP-Filter. 
%
% This warning only applies in case of DO_HP_FILTER = 1:
% WARNING: BECAUSE APPLYING THE HP_FILTER IMPLIES SOLVING A LINEAR SYSTEM OF EQUATIONS OF
% SIZE SIM_LENGTH, IT IS PROBABLY NOT A GOOD IDEA TO HAVE SIM_LENGTH GREATER THAN, SAY, 1000.
% OTHERWISE YOU EASILY CAN ENCOUNTER MEMORY PROBLEMS, AND THE PROCEDURE TAKES A LOT OF TIME.
% It is numerically more sensible to use, say, SIM_LENGTH = 100 and SIM_N_SERIES = 10
% and SIM_MODE = 2 instead.
% 
% output is produced by SIM_OUT.M.


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

% Copied from IMPRESP.M:
II_lag = [ PP, zeros(m_states,n_endog),zeros(m_states,k_exog)
           RR, zeros(n_endog, n_endog),zeros(n_endog, k_exog)
           zeros(k_exog,(m_states+n_endog)), NN                ];
II_contemp = eye(m_states+n_endog+k_exog) + ...
     [ zeros(m_states,(m_states+n_endog)), QQ
       zeros(n_endog, (m_states+n_endog)), SS
       zeros(k_exog,  (m_states+n_endog)), zeros(k_exog,k_exog) ];
% describing [x(t)',y(t)',z(t)']'= II_contemp*II_lag*[x(t-1)',y(t-1)',z(t-1)']';



% The following piece is due to Gerard A. Pfann
if DO_HP_FILTER,
   HP_mat = [1+HP_LAMBDA, -2*HP_LAMBDA, HP_LAMBDA,              zeros(1,SIM_LENGTH-3);
             -2*HP_LAMBDA,1+5*HP_LAMBDA,-4*HP_LAMBDA,HP_LAMBDA, zeros(1,SIM_LENGTH-4);
                           zeros(SIM_LENGTH-4,SIM_LENGTH);
              zeros(1,SIM_LENGTH-4),HP_LAMBDA,-4*HP_LAMBDA,1+5*HP_LAMBDA,-2*HP_LAMBDA;     
              zeros(1,SIM_LENGTH-3),          HP_LAMBDA,   -2*HP_LAMBDA, 1+HP_LAMBDA  ];
   for i=3:SIM_LENGTH-2;
     HP_mat(i,i-2)=HP_LAMBDA;
     HP_mat(i,i-1)=-4*HP_LAMBDA;
     HP_mat(i,i)=1+6*HP_LAMBDA;
     HP_mat(i,i+1)=-4*HP_LAMBDA;
     HP_mat(i,i+2)=HP_LAMBDA;
   end;
end;
% If x is a column vector of length SIM_LENGTH,
% xtr=HP_mat\x; delivers the HP-trend and
% xhp=x-xtr; delivers the HP-filtered series


if SIM_TRACK_N & (SIM_MODE == 2),
   message = ['SIMUL.M: I will inform you about starting each new simulation.         '
              '         Turn this feature off with SIM_TRACK_N = 0.                   '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;
if SIM_TRACK_LENGTH & (SIM_LENGTH >= SIM_TRACK_FAC),
   message = ['SIMUL.M: I will inform you about my progress as I calculate each       '
              '         simulation.  Turn this feature off with SIM_TRACK_LENGTH = 0. '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;

n_select = max(size(SIM_SELECT));
Sig_fac= chol(Sigma);  % upper triangular Cholesky-factor
Sig_fac = Sig_fac';     % to get the lower triangular Cholesky factor
if SIM_MODE == 2,
   autcor_sum = zeros(n_select+1,(2*N_LEADS_LAGS-1));
   autcor_sqr = zeros(n_select+1,(2*N_LEADS_LAGS-1));
   covmat_sum = zeros(n_select,n_select);
   covmat_sqr = zeros(n_select,n_select);
   stdvec_sum = zeros(n_select,1);
   stdvec_sqr = zeros(n_select,1);
   n_sim = SIM_N_SERIES;
   message = ['SIMUL.M: I will now calculate simulation-based moments.                '
              '         The advantage: you get small-sample standard errors.          '
              '         The disadvantage: it may take a LONG time!!                   '
              '         If you hate this, you should turn me off by either reducing   '
              '         SIM_N_SERIES or SIM_LENGTH or by setting SIM_MODE = 1 or,     '
              '         more brutaly, by setting DO_SIMUL = 0.                        '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
elseif SIM_MODE == 1,
   n_sim = 1;
else
   message = ['SIMUL.M: You have set SIM_MODE to a nonsensical values: allowed are    '
              '         only the values 1 or 2.  Type help simul for information.     '
              '         I will proceed as if you had set SIM_MODE = 1.                '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];
end;

for sim_index = 1 : n_sim, 
   if (SIM_MODE == 2) & SIM_TRACK_N,
      disp(sprintf('SIMUL.M: starting simulation %5d out of %5d...',sim_index,n_sim));
   end;
   if SIM_RANDOM_START,
      xyz = zeros(m_states+n_endog+k_exog,1);
      if SIM_GIVEN_DISCARD and (SIM_MODE ~= 2),
         if sim_index == 1,
            message = ['SIMUL.M: I am using your given shocks epsilon.  You should have defined'
                       '         discard_eps as a matrix of size k_exog x SIM_DISCARD.         '
                       '         You can turn this feature off with SIM_GIVEN_DISCARD = 0.     '];
            if DISPLAY_IMMEDIATELY, disp(message); end;
            warnings = [warnings;message];
         end;
         sim_discard_eps  = discard_eps;
      else
         sim_discard_eps = Sig_fac*randn(k_exog,SIM_DISCARD);
      end;
      for t = 1 : SIM_DISCARD,
         xyz =II_contemp*(II_lag*xyz + ...
                             [zeros(m_states+n_endog,1);sim_discard_eps(:,t)]);
      end;
   else
      xyz = II_contemp*(II_lag*[SIM_X_START;zeros(n_endog+k_exog,1)]+...
                             [zeros(m_states+n_endog,1);SIM_Z_START]);
   end;
   sim_xyz = zeros(m_states+n_endog+k_exog,SIM_LENGTH);
   if SIM_GIVEN_EPS & (SIM_MODE ~= 2), % CHANGE IN 2002: Replaced 'and' with '&'
      if sim_index == 1,
         message = ['SIMUL.M: I am using your given shocks epsilon.  You should have defined'
                    '         given_eps as a matrix of size k_exog x SIM_LENGTH.            '
                    '         You can turn this feature off with SIM_GIVEN_EPS = 0.         '];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
      end;
      sim_eps  = given_eps;
   else
      sim_eps = Sig_fac*randn(k_exog,SIM_LENGTH);
   end;
   sim_xyz(:,1) = xyz;
   for t = 2 : SIM_LENGTH,
      if SIM_TRACK_LENGTH & (t == floor(t/SIM_TRACK_FAC)*SIM_TRACK_FAC),
         disp(sprintf('SIMUL.M:      Simulation date t = %10d out of %10d...',t,SIM_LENGTH));
      end;
      xyz =II_contemp*(II_lag*xyz + ...
                          [zeros(m_states+n_endog,1);sim_eps(:,t)]);
      sim_xyz(:,t) = xyz;
   end;
   sim_raw = sim_xyz;
   if DO_HP_FILTER,
         if SIM_TRACK_LENGTH | (SIM_TRACK_N & (SIM_MODE == 2)),
            disp('SIMUL.M:      Doing the HP-Filter...');
         end;
      sim_xyz = sim_raw' - HP_mat\sim_raw';
      sim_xyz = sim_xyz';
   end;  
   covmat_one = zeros(n_select,n_select);
   for j = 1 : n_select,
      sim_ser = ones(n_select,1)*sim_xyz(SIM_SELECT(j),:);
      covmat_j = sum( sim_ser'.*sim_xyz(SIM_SELECT,:)' );
      covmat_one(:,j) = covmat_j';
   end;
   covmat_one = covmat_one/SIM_LENGTH;
   % Note that the expected mean of sim_xyz is zero by construction!
   % Of course, the realized mean of sim_xyz may be different from zero.
   stdvec_one = sqrt(diag(covmat_one)); 
   sim_gnp = [zeros(n_select,N_LEADS_LAGS-1),...
              ones(n_select,1)*sim_xyz(GNP_INDEX,:),...
              zeros(n_select,N_LEADS_LAGS-1)];
   autcor_one = zeros(n_select+1,(2*N_LEADS_LAGS-1));
   for j = 1 : 2*N_LEADS_LAGS-1,
      autcor_j = sum( sim_gnp(:,j:(j-1+SIM_LENGTH))'.*sim_xyz(SIM_SELECT,:)')/...
                    (SIM_LENGTH - abs(N_LEADS_LAGS - j))./...
                        (stdvec_one' * stdvec_one(GNP_INDEX));
      % corresponding to corr(ser(t),gnp(t+j-N_LEAD_LAGS))
      % which is the same as corr(ser(t-j+N_LEAD_LAGS,gnp(t))
      % The latter form is used in the tables, so:
      %
      % OLD VERSION, before May 2003:
      % autcor_one(:,j) = [autcor_j';(j-N_LEADS_LAGS)];
      % NEW VERSION, after May 2003:
      autcor_one(:,2*N_LEADS_LAGS-j) = [autcor_j';(N_LEADS_LAGS-j)];
   end;
   if SIM_MODE == 2,
      autcor_sum = autcor_sum + autcor_one;
      covmat_sum = covmat_sum + covmat_one;
      stdvec_sum = stdvec_sum + stdvec_one;
      autcor_sqr = autcor_sqr + autcor_one.^2;
      covmat_sqr = covmat_sqr + covmat_one.^2;
      stdvec_sqr = stdvec_sqr + stdvec_one.^2;
   end;               
end;
if SIM_MODE == 2,
   autcor_sim = autcor_sum/n_sim;
   covmat_sim = covmat_sum/n_sim;
   stdvec_sim = stdvec_sum/n_sim;
   autcor_std = sqrt( (autcor_sqr - n_sim* autcor_sim.^2)/(n_sim-1) );
   covmat_std = sqrt( (covmat_sqr - n_sim* covmat_sim.^2)/(n_sim-1) );
   stdvec_std = sqrt( (stdvec_sqr - n_sim* stdvec_sim.^2)/(n_sim-1) );
else %  SIM_MODE == 1,
   autcor_sim = autcor_one;
   covmat_sim = covmat_one;
   stdvec_sim = stdvec_one;
end;
  
      
      