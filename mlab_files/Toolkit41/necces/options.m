% VERSION 4.1, MAY 2003, COPYRIGHT H. UHLIG.
%   May 2003: QZ-algorithm is know the default solution method.
%
% OPTIONS.M sets the options to be used in all the subroutines.
% It is called by DO_IT.M.  If you use any of the routines
% "by hand", i.e. without calling them via DO_IT.M, you also
% have to execute OPTIONS.M first.
% If you want to change some options, change them here.
%
% Some options can also set directly inside the programs: they are then not overwritten
% here.  They are only set here to their default values, if they have not been set previously.
% Of course, if you didn't set them inside the program, but they have been set by
% a previous program, which you have executed before, then the values set by this
% previous program will not be overwritten either and thus used.  I.e. IT IS GENERALLY A
% GOOD IDEA TO TYPE CLEAR BEFORE RUNNING YOUR PROGRAM  SO THAT THE OPTIONS WILL ALWAYS BE
% SET THE WAY YOU WANT THEM TO.  Alternatively, if you have run your program once,
% just change some options, and rerun the calculations by typing do_it
%
% The options, which can be set within a program or by hand, are, 
% together with their default values:
%
%       IT IS BAD STYLE NOT TO SET THESE OPTIONS WITHIN YOUR PROGRAM:
%   PERIOD     = 4;    % number of periods per year, i.e. 12 for monthly, 4 for quarterly
%   GNP_INDEX  = 2;    % Index of output among the variables selected for HP filter
%
%       IT IS OK TO TRUST THE DEFAULT SETTINGS FOR THE FOLLOWING OPTIONS:
%   DISPLAY_IMMEDIATELY = 0; % = 1: displays warning messages immediately.
%                                     % DO THIS, IF YOUR CALCULATIONS MYSTERIOUSLY CRASH!
%   DISPLAY_LATER = 1;                % = 1: displays warning messages after SOLVE.M
%                                     % DO THIS, IF YOUR CALCULATIONS MYSTERIOUSLY CRASH!
%   DISPLAY_AT_THE_END = 1;           % = 1: displays warning messages at the end of all computations.
%   MANUAL_ROOTS = 0;  % = 1, if you want to choose your own
%                      % roots, otherwise set = 0. 
%                      % See SOLVE.M for instructions.
%   DO_QZ = 1;         % Set = 0, if you want to employ the generalized eigenvector/eigenvalue 
%                      % decomposition to find the matrix P (called PP in the programs).
%   IGNORE_VV_SING = 1;% =1: Ignores, if VV is singular.  
%                      % Sometimes useful for sunspots.  Cross your fingers...  
%   DO_IMPRESP = 1;    % = 1, if you want to calculate impulse responses.
%   DO_MOMENTS = 0;    % Turn on, if you want to do
%                      % fourier-transforms-based calculations of moments
%   DO_SIMUL = 0;      % Turn on, if you want to do simulations or 
%                      % simulation-based calculations of moments
%
%       Only relevant, if DO_MOMENTS = 1 or DO_SIMUL = 1:
%   DO_HP_FILTER= 1;   % Set to = 1, if you want moments and simulations with HP-filter.
%                      % set to = 0, if you want moments and simulations without HP-filter.
%   HP_LAMBDA = 1600*(PERIOD/4)^4;    % The lambda parameter for the HP-Filter
%   HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calculations.
%   N_LEADS_LAGS = 6;  % calculate autocorr. table up to N_LEADS_LAGS-1 leads and lags.
%   SIM_MODE = 1;      % Set to = 1, if you just want one simulated time series.
%                      % Set to = 2, if you want to do simulation-based calculation
%                      % of moments. 
%   SIM_LENGTH  = 100; % Length of time series to be simulated. 
%   SIM_RANDOM_START = 1; 
%                      % = 1, if you want a random start, = 0, if you want predetermined start
%   SIM_X_START = zeros(m_states,1);
%                      % if predetermined start, SIM_X_START should be a vector of length
%                      % m, setting the initial values for the endogenous state variables x.
%   SIM_Z_START = zeros(k_exog,1); 
%                      % if predetermined start, SIM_Z_START should be a vector of length
%                      % k, setting the initial values for the exogenous state variables z.
%   SIM_GIVEN_EPS = 0; % set = 1, if you want to choose the shocks epsilon yourself.
%                      % In that case, declare them beforehand as
%                      % given_eps 
%                      % which needs to be a matrix of size k_exog x SIM_LENGTH
%
%                      % ONLY RELEVANT FOR SIM_MODE = 2.
%   SIM_N_SERIES = 50; % number of simulated series for the calculations of moments.
%                      % This should be set to a value of at least 4 to get 
%                      % small sample standard errors for the computed moments.
%                      % The moments computed from the data should then lie within a
%                      % chosen confidence interval to be calculated from the moments
%                      % and their small-sample standard errors from the simulations.
%                      % ONLY RELEVANT FOR SIM_MODE = 2.
%   SIM_TRACK_N = 1;   % Prints a little message at the start of each simulation.
%   SIM_TRACK_LENGTH = 1; % Prints a little message every SIM_TRACK_FAC date of
%                      % a simulation.  SIM_TRACK_FAC = 1000 is set in OPTIONS.M.
%                      % It also prints a message when applying the HP-Filter.
%
%       The following are options for modifying graphs:
%   DO_SHOCK_RESP = 1; % = 1, if impulse responses to shocks shall be
%                      % calculated.
%   SELECT_SHOCKS = 1 : k_exog; 
%                      % select the shocks to which impulse responses should be plotted.
%   DO_STATE_RESP = 1; % = 1, if impulse responses to deviations of state
%                      % variables should be calculated.
%   SELECT_STATES = 1 : m_states; 
%                      % select the states to which impulse responses should be plotted.
%   INIT_DATE = 4;     % Number of periods prior to shock.  Default = 4.
%                      % For INIT_DATE = 1, one sees the inital value for the state variable.
%                      % For INIT_DATE > 0, variables are initially zero, before responding.
%   HORIZON      = 32; % how far out should the impulse responses be calculated
%   DO_PLOTS     = 1;  % If impulse response plots should be made, = 0, if not.   
%   IMP_SUBPLOT  = 0;  % Set =1, if impulse respones should be plotted in a subplot   
%   IMP_SUB_FONT = 6;  % Size of the title and labels in the subplots
%   IMP_JOINT    = 0;  % Set =1, if all responses should be plotted in one graph, =0, if not 
%   IMP_SINGLE   = 1;  % If each response should be plotted in a separate graph, =0 if not.
%   IMP_SELECT = 1:(m_states+n_endog+k_exog);
%                      %  a vector containing the indices of the variables to be plotted
%   DO_NO_ZERO_RESPONSE =  0; 
%                      % Set = 1 to avoid plotting responses near zero
%   ZERO_RESPONSE_LEVEL =  0.0001;   
%                      % Defining near zero to be the interval 
%                      % [-ZERO_RESPONSE_LEVEL,ZERO_RESPONSE_LEVEL]
%   DO_SCALE_AXIS = 0; % if the axis should be rescaled. 
%                      % Choices: according to a given scale
%                      %       or according to the max. response among AXIS_SELECT
%   DO_GIVEN_AXIS_SCALE = 0; % Scale according to a given axis scale, see next variable.
%                      % Only use = 1, if you are familiar with scaling axis in MATLAB.
%                      % see help axis for more information.
%   AXIS_SCALE = [0,HORIZON,-1,1];   % The argument for the axis(.) command, 
%                      % normally a vector of length 4 = [xmin,xmax,ymin,ymax].
%   AXIS_SELECT= 1:(m_states+n_endog+k_exog);
%                      %  for rescaling the axis automatically according to the
%                      %  max. response in AXIS_SELECT.
%                      % It is a good idea to make this at least as large as IMP_SELECT.
%   DO_ENLARGE = 0;    % =1, if you need large letter sizes for overheads, say.
%                      % If set to 1, letter size will be manipulated by ENLARGE.M
%   DO_THICK_LINES = 1;% If you set DO_ENLARGE = 1, then you will also get thick lines with DO_THICK_LINES = 1;
%   PRINT_FIG  = 0;    % Set = 1, if you want to sent figures directly to printer
%   DO_COLOR_PRINT = 0;% Set = 1 together with PRINT_FIG = 1, if you want to sent figures
%                      % directly to the currently installed windows color printing device
%   SAVE_FIG   = 0;    % Set = 1, if you want to save figures as encapsulated postscript files.
%                      % The graphic names are automatically chosen.  E.g. for impulse responses,
%                      % the file names are imprespN.eps, where N = 1,2,...  .  
%                      % Note: you need to set PRINT_FIG = 0 to be able to use SAVE_FIG.
%   COLOR_FIG  = 1;    % Set = 0, if you want to save as black-and-white .eps files (old default)
%                      % Set = 1, if you want to save as color .eps files (new default)
%                      % Note: you need to set SAVE_FIG = 1 or SAVE_FIG_JPG = 1 and PRINT_FIG = 0.
%   SAVE_FIG_JPG   = 0;% Set = 1, if you want to save figures as color jpg-files.
%                      % The graphic names are automatically chosen.  E.g. for impulse responses,
%                      % the file names are imprespN.jpg, where N = 1,2,...  .  
%                      % Note: you need to set PRINT_FIG = 0 to be able to use SAVE_FIG_JPG.
%   FILENAMEPIECE = '';% To change the naming of the impulse response files from e.g.
%                      % impresp1.eps to impresp_new_1.eps,
%                      % assign FILENAMEPIECE='_new_'.
%   DO_FIXED_TEXT_POS  = 0; 
%                      % If text labels should be put at TXT_MARKER, set = 1.
%                      % Alternative: use point of maximal response between zero and MAX_TXT_MARKER.
%   TXT_MARKER  = 3;   % Where to put the labels on the graph, a number between 1 and HORIZON/2.
%   MAX_TXT_MARKER  = round(HORIZON/2); 
%                      % Maximal position for text and search for max. response.
%         Only relevant if DO_SIMUL = 1:
%   SIM_GRAPH = 1;     % Set to = 1 to see plots of the simulated series. 
%   SIM_SUBPLOT=0;     % Set =1, if simulated series should be plotted in a subplot.
%   SIM_JOINT=1;       % If all simulated series should be plotted in one graph, =0, if not 
%   SIM_SINGLE=0;      % set =1, if each simulated series should be plotted in one graph.
%   SIM_SUB_FONT=6;    % Size of the title and labels in the subplots
%   SIM_MAX = 200;     % Maximal number of dates to be used for the plots.
%   SIM_SELECT = 1:(m_states+n_endog+k_exog); 
%                      % Selecting the variables for the graphics of simulations.
%   SIM_DATE0  = 0;    % Set to the initial date of the simulated series.
%         Only relevant if DO_MOMENTS = 1:
%   DO_HP_GRAPH=1;     % To plot spectral densities
%   MOM_SUBPLOT=0;     % Set =1, if spectral densities should be plotted in a subplot
%   MOM_JOINT=1;       % If all spectral densities should be plotted in one figure, =0, if not
%   MOM_SINGLE=0;      % Set =1, if spectral densities should be plotted in separate figures
%   MOM_SUB_FONT=6;    % Size of the title and labels in the subplots
%   MOM_PLOT_RAW=0;    % If you like non-HP-filtered spectral density plots

% CHANGING ANY OPTION OTHER THAN THE ONES LISTED ABOVE REQUIRES CHANGING THE FILE OPTIONS.M


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

if exist('AA')==1,
   [l_equ,m_states] = size(AA);
   [l_equ,n_endog ] = size(CC);
   [l_equ,k_exog  ] = size(DD);
end;

%============> GENERAL OPTIONS:
if exist('PERIOD')~=1,
   PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
end;
if exist('GNP_INDEX')~=1,
   GNP_INDEX = 2,   % Index of output among the variables selected for HP filter
end;
if exist('DISPLAY_IMMEDIATELY')~=1,
   DISPLAY_IMMEDIATELY = 0; % = 1: displays warning messages immediately.
end;
if DISPLAY_IMMEDIATELY,
   message = ['OPTIONS.M: You will receive warning messages always immediately.       '
              '           You can turn that off with DISPLAY_IMMEDIATELY = 0          '
              '           Warning messages are stored in the array "warnings"         '];
   disp(message);
   warnings = [warnings;message];   
else
   message = ['OPTIONS.M: You will not receive warning messages immediately.          '
              '           To turn on immediate warnings, set DISPLAY_IMMEDIATELY = 1  '
              '           Warning messages are stored in the array "warnings"         '];
   disp(message);
   warnings = [warnings;message];   
end;
message = ['OPTIONS.M: There are lots of options to choose from.  Type help options'
           '           to get more information about them!                         '];
if DISPLAY_IMMEDIATELY, disp(message); end;
warnings = [warnings;message];   
if exist('DISPLAY_LATER')~=1,
   DISPLAY_LATER = 0;                % = 1: displays warning messages after SOLVE.M
end;
if exist('DISPLAY_AT_THE_END')~=1,
   DISPLAY_AT_THE_END = 1;           % = 1: displays warning messages at the end of all computations.
end;
if exist('DO_HP_FILTER')~=1,
   DO_HP_FILTER= 1;    % Set to = 1, if you want moments and simulations with HP-filter.
                       % set to = 0, if you want moments and simulations without HP-filter.
end;
if exist('HP_LAMBDA')~=1,            % The lambda parameter for the HP-Filter:
   HP_LAMBDA = 1600*(PERIOD/4)^4;    % This is not the literature standard, 
                                     % but I believe it is correct.
   %  HP_LAMBDA = 1600*(PERIOD/4)^2; % This is the literature standard
end;
if exist('HP_SELECT')~=1,
   HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
end;
if exist('N_LEADS_LAGS')~=1,
   N_LEADS_LAGS = 6;       % calculate autocorr. table up to N_LEADS_LAGS-1 leads and lags.
end;
if exist('DO_ENLARGE')~=1, 
   DO_ENLARGE = 0; % =1, if you need large letter sizes for overheads, say.
                   % If set to 1, letter size will be manipulated by
                   % the program ENLARGE.M
end;
if exist('DO_THICK_LINES')~=1, 
   DO_THICK_LINES = 1; % If you set DO_ENLARGE = 1, then you will also get thick lines with DO_THICK_LINES = 1;
end;
if exist('PRINT_FIG')~=1, 
   PRINT_FIG  = 0; % Set = 1, if you want to sent figures directly to printer
end;
if exist('DO_COLOR_PRINT')~=1,
   DO_COLOR_PRINT = 0; % Set = 1, and set PRINT_FIG = 1, if you want to send
   % your figures to the currently installed windows color printing device.
end;
if exist('SAVE_FIG')~=1, 
   SAVE_FIG   = 0; % Set = 1, if you want to save figures as encapsulated postscript files.
                   % The graphic names are automatically chosen.  E.g. for impulse responses,
                   % the file names are imprespN.eps, where N = 1,2,...  .  
                   % Note: you need to set PRINT_FIG = 0 to be able to use SAVE_FIG.
end;
if exist('SAVE_FIG_JPG')~=1, 
   SAVE_FIG_JPG   = 0; % Set = 1, if you want to save figures as .jpg-files
                   % The graphic names are automatically chosen.  E.g. for impulse responses,
                   % the file names are imprespN.eps, where N = 1,2,...  . 
                   % Note: you need to set PRINT_FIG = 0 to be able to use SAVE_FIG_JPG.
end;
if exist('COLOR_FIG')~=1, 
   COLOR_FIG   = 1;  % Set = 0, if you want to save as black-and-white .eps files (old default)
                     % Set = 1, if you want to save as color .eps files (new default)
                     % Note: you need to set SAVE_FIG = 1 or SAVE_FIG_JPG = 1 and PRINT_FIG = 0.
end;


%============> options for DO_IT.M
if exist('DO_IMPRESP')~=1,
   DO_IMPRESP = 1; % =1, if impulse responses should be calculated.
                   % To get their plots, also DO_PLOTS = 1 and
                   % at least one of DO_SHOCK_RESP = 1 or 
                   % DO_STATE_RESP = 1 needs to be set, see below.
end;
if exist('DO_MOMENTS')~=1,
   DO_MOMENTS = 0; % Turn on, if you want to do
                   % fourier-transforms-based calculations of moments
end;
if DO_MOMENTS==1,
   message = ['OPTIONS.M: Note, that I calculate moments, which takes time. If you    '
              '           do not want to have them calculated, set DO_MOMENTS = 0     '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
else
   message = ['OPTIONS.M: Note, that I will not calculate moments.                    '
              '           If you want to have them calculated, set DO_MOMENTS = 1     '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
end;

if exist('DO_SIMUL')~=1,
  DO_SIMUL = 0;    % Turn on, if you want to do simulations or 
                   % simulation-based calculations of moments
end;
if DO_SIMUL,
   message = ['OPTIONS.M: Note, that I calculate simulations, which takes time.       '
              '           If you do not want to have them calculated, set DO_SIMUL = 0'];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
else
   message = ['OPTIONS.M: Note, that I will not calculate simulations.                '
              '           If you want to have them calculated, set DO_SIMUL = 1       '];
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
end;
if DO_SIMUL | DO_MOMENTS,
   if DO_HP_FILTER,
      if HP_LAMBDA >= 1e+11
         line1 = sprintf('OPTIONS.M: I will use the HP-Filter with lambda=%12.4e.          ',HP_LAMBDA);
      elseif HP_LAMBDA >= 1000,
         line1 = sprintf('OPTIONS.M: I will use the HP-Filter with lambda=%12.0f.          ',HP_LAMBDA);
      elseif HP_LAMBDA >= 1,
         line1 = sprintf('OPTIONS.M: I will use the HP-Filter with lambda=%12.4f.          ',HP_LAMBDA);
      elseif HP_LAMBDA >= .0001,
         line1 = sprintf('OPTIONS.M: I will use the HP-Filter with lambda=%12.8f.          ',HP_LAMBDA);
      else
         line1 = sprintf('OPTIONS.M: I will use the HP-Filter with lambda=%12.4e.          ',HP_LAMBDA);
      end;
      message = [line1
                 '           If you want to change that, modify DO_HP_FILTER and         '
                 '           HP_LAMBDA in OPTIONS.M or in your program or by changing    '
                 '           e.g.HP_LAMBDA now and restarting the calculations with do_it'];
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
   else
      message = ['OPTIONS.M: I will NOT use the HP-Filter.                               '
                 '           If you want to change that, modify DO_HP_FILTER and         '
                 '           HP_LAMBDA in OPTIONS.M or in your program or by modifying   '
                 '           them now and restarting your calculations by typing do_it   '];
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
   end;
end;
   


%============> options for SOLVE.M
if exist('TOL')~=1,
   TOL = .000001; 
end;           % Roots smaller than TOL are regarded as zero.
               % Complex numbers with distance less than TOL are regarded as equal.
if exist('MANUAL_ROOTS')~=1,
   MANUAL_ROOTS = 0; % = 1, if you want to choose your own
                               % roots, otherwise set = 0. 
                               % See SOLVE.M for instructions.
end;
if exist('DO_QZ')~=1,
   DO_QZ = 1; % = 0, if you want to employ the generalized eigenvector/eigenvalue 
%             % decomposition to find the matrix P (called PP in the programs).
end;
if exist('IGNORE_VV_SING')~=1,
   IGNORE_VV_SING = 1; % =1: Ignores, if VV is singular.  
                       % Sometimes useful for sunspots.  Cross your fingers...  
end;


%============> options for SOL_OUT.M
DISPLAY_ROOTS = 1;  % Set = 1, if you want to see the roots.



%============> options for IMPRESP.M
if exist('DO_PLOTS')~=1,
   DO_PLOTS   = 1; % if impulse response plots should be made, = 0, if not.
end;
if DO_PLOTS,
   message = ['OPTIONS.M: Note, that I plot impulse responses, which you may dislike. '
              '           If you do not want to have them plotted, set DO_PLOTS = 0   '];
   warnings = [warnings;message];
   message = ['OPTIONS.M: For imp.responses you can choose between joint-, sub-, and  ' 
              '           single plots. Try, SIM_JOINT=1, SIM_SUBPLOT=1, SIM_SINGLE=1 ']; 
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
else
   message = ['OPTIONS.M: Note, that I will not plot impulse responses.               '
              '           If you want to have them plotted, set DO_PLOTS = 1          '];
   warnings = [warnings;message];
   message = ['OPTIONS.M: For imp.responses you can choose between joint-, sub-, and  ' 
              '           single plots. Try, SIM_JOINT=1, SIM_SUBPLOT=1, SIM_SINGLE=1 '];  
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
end;
if exist('IMP_SUBPLOT')~=1, % if impulse respones should be plotted in a subplot 
    IMP_SUBPLOT=0;
end;
if exist('IMP_JOINT')~=1, % if all responses should be plotted in one graph, =0, if not 
    IMP_JOINT=1;
end; 
 
if exist('IMP_SINGLE')~=1, % if each response should be plotted in a single graph, =0, if not 
    IMP_SINGLE=0;
end; 
if exist('IMP_SUB_FONT')~=1,  % Size of the title and labels in the subplots
    IMP_SUB_FONT=6;
end;
if exist('INIT_DATE')~=1,
   INIT_DATE    = 4; % number of periods prior to shock.
end;
if exist('HORIZON')~=1,
   HORIZON    = 32;% how far out should the impulse responses be calculated\
end;
if exist('IMP_SELECT')~=1,
   IMP_SELECT = 1:(m_states+n_endog+k_exog);
                %  a vector containing the indices of the variables to be plotted
end;
if exist('DO_NO_ZERO_RESPONSE')~=1,
   DO_NO_ZERO_RESPONSE =  0;
                % Set = 1 to avoid plotting responses near zero
end;
if exist('ZERO_RESPONSE_LEVEL')~=1,
   ZERO_RESPONSE_LEVEL =  0.0001;
                % Defining near zero to be the interval 
                % [-ZERO_RESPONSE_LEVEL,ZERO_RESPONSE_LEVEL]
end;
if exist('DO_SCALE_AXIS')~=1,
   DO_SCALE_AXIS =  0;
                %  if the y-axis should be scaled according to different set of variables
end;
if exist('DO_GIVEN_AXIS_SCALE')~=1,
   DO_GIVEN_AXIS_SCALE =  0;
                % Scale according to a given axis scale, see next variable.
                % Only use = 1, if you are familiar with scaling axis in MATLAB.
                % see help axis for more information.
end;
if exist('AXIS_SCALE')~=1,
   AXIS_SCALE =  [0,HORIZON,-1,1];
                % The argument for the axis(.) command, 
                % normally a vector of length 4 = [xmin,xmax,ymin,ymax].
end;
if exist('AXIS_SELECT')~=1,
   AXIS_SELECT = 1:(m_states+n_endog+k_exog);
                %  a vector containing the indices of the variables to which to scale the axis
end;
% To change the naming of the impulse response files from e.g.
% impresp1.eps to impresp_new_1.eps,
% assign FILENAMEPIECE='_new_'.
% Default is FILENAMEPIECE = '';
if ~exist('FILENAMEPIECE'),
    FILENAMEPIECE = ''; 
end;
if exist('DO_SHOCK_RESP')~=1,
   DO_SHOCK_RESP = 1; % = 1, if impulse responses to shocks should be plotted.
end;
if exist('SELECT_SHOCKS')~=1,
   SELECT_SHOCKS = 1 : k_exog; 
                      % select the shocks to which impulse responses should be plotted.
end;
if exist('DO_STATE_RESP')~=1,
   DO_STATE_RESP = 1; % = 1, if impulse responses to shocks should be plotted.
end;
if exist('SELECT_STATES')~=1,
   SELECT_STATES = 1 : m_states; 
                      % select the states to which impulse responses should be plotted.
end;
if exist('DO_FIXED_TEXT_POS')~=1,
   DO_FIXED_TEXT_POS  = 0; % If text labels should be put at TXT_MARKER, set = 1.
       % Alternative: use point of maximal response between zero and MAX_TXT_MARKER.
end; 
if exist('TXT_MARKER')~=1,
   TXT_MARKER  = 3; % Where to put the labels on the graph, a number between 1 and HORIZON/2.
end; 
if exist('MAX_TXT_MARKER')~=1,
   MAX_TXT_MARKER  = round(HORIZON/2); % Maximal position for text and search for max. response.
end; 
%============> options for ENLARGE.M
LET_SIZE = 18;  % Size of letters in text on axes and title.
LAB_SIZE = 24;  % Size of letters in labels
% OLD_LET_SIZE = 12; % Default letter size of MATLAB. - NO LONGER NEEDED, H. Uhlig, May 2001
% OLD_LAB_SIZE = 12; % Default letter size of MATLAB. - NO LONGER NEEDED, H. Uhlig, May 2001


%============> options for SIMUL.M
if exist('SIM_MODE')~=1,
   SIM_MODE = 1;    % Set to = 1, if you just want one simulated time series.
                    % Set to = 2, if you want to do simulation-based calculation
                    % of moments.
end;
if exist('SIM_LENGTH')~=1,
   SIM_LENGTH  = 100;  % Length of time series to be simulated.
end;  
if exist('SIM_RANDOM_START')~=1,
   SIM_RANDOM_START = 1; 
                    % = 1, if you want a random start, = 0, if you want predetermined start
end;
if exist('SIM_X_START')~=1,
   SIM_X_START = zeros(m_states,1);
                    % if predetermined start, SIM_X_START should be a vector of length
                    % m, setting the initial values for the endogenous state variables x.
end;
if exist('SIM_Z_START')~=1,
   SIM_Z_START = zeros(k_exog,1); 
                    % if predetermined start, SIM_X_START should be a vector of length
                    % k, setting the initial values for the exogenous state variables z.
end;
SIM_DISCARD = 100;  % Length of initial discarded observations, if random start.
if exist('SIM_GIVEN_EPS')~=1,
   SIM_GIVEN_EPS = 0; % set = 1, if you want to choose the shocks epsilon yourself.
                      % In that case, declare them beforehand as
                      % given_eps 
                      % which needs to be a matrix of size k_exog x SIM_LENGTH
end;
SIM_GIVEN_DISCARD = 0;% set = 1, if you want to choose the shocks epsilon 
                      % for the initially discarded observations (SIM_RANDOM_START = 1)
                      % yourself.
                      % In that case, declare them beforehand as
                      % discard_eps 
                      % which needs to be a matrix of size k_exog x SIM_DISCARD
if exist('SIM_N_SERIES')~=1,
   SIM_N_SERIES = 50;  % number of simulated series for the calculations of moments.
                       % This should be set to a value of at least 4 to get 
                       % small sample standard errors for the computed moments.
                       % The moments computed from the data should then lie within a
                       % chosen confidence interval to be calculated from the moments
                       % and their small-sample standard errors from the simulations.
                       % ONLY RELEVANT FOR SIM_MODE = 2.
end;
if exist('SIM_TRACK_N')~=1,
   SIM_TRACK_N = 1;   % Prints a little message at the start of each simulation.
                      % ONLY RELEVANT FOR SIM_MODE = 2.
end;
if exist('SIM_TRACK_LENGTH')~=1,
   SIM_TRACK_LENGTH = 1; % Prints a little message every SIM_TRACK_FAC date of
                      % a simulation.  
end;
SIM_TRACK_FAC = 1000;


%============> options for SIM_OUT.M
if exist('SIM_GRAPH')~=1,
   SIM_GRAPH = 1;   % Set to = 1 to see plots of the simulated series.
end;
if (SIM_GRAPH & DO_SIMUL),
   message = ['OPTIONS.M: Note, that I plot simulations, which you may dislike.       '
              '           If you do not want to have them plotted, set SIM_GRAPH=0    '];
   warnings = [warnings;message];
   message = ['OPTIONS.M: For simulations you can choose between joint-, sub-, and    ' 
              '           single plots. Try, SIM_JOINT=1, SIM_SUBPLOT=1, SIM_SINGLE=1 ']; 
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
else   
   message = ['OPTIONS.M: Note, that I will not plot simulations.                     '
              '           If you want to have them plotted, set SIM_GRAPH=1           '];
   warnings = [warnings;message];
   message = ['OPTIONS.M: For simulations you can choose between joint-, sub-, and    ' 
              '           single plots. Try, SIM_JOINT=1, SIM_SUBPLOT=1, SIM_SINGLE=1 ']; 
   if DISPLAY_IMMEDIATELY, 
       disp(message); 
   end;
   warnings = [warnings;message];   
end;

if exist('SIM_SUBPLOT')~=1, % if simulated series should be plotted in a subplot, =0, if not
    SIM_SUBPLOT=0;
end;
if exist('SIM_JOINT')~=1, % if all simulated series should be plotted in one graph, =0, if not 
    SIM_JOINT=1;
end; 
if exist('SIM_SINGLE')~=1, % if each simulated series should be plotted in one graph, =0, if not 
    SIM_SINGLE=0;
end; 
if exist('SIM_SUB_FONT')~=1,  % Size of the title and labels in the subplots
    SIM_SUB_FONT=6;
end;
if exist('SIM_SELECT')~=1,
   SIM_SELECT = 1:(m_states+n_endog+k_exog); 
                      % Selecting the variables for the graphics of simulations.
end;
if exist('SIM_DATE0')~=1,
   SIM_DATE0  = 0;    % Set to the initial date of the simulated series.
end;
if exist('SIM_MAX')~=1,
   SIM_MAX         = 200; % Maximal amount of dates to be plotted of a simulation.
end;
SIM_DO_DISP1    = 1;% Set to = 1 to see printout of the autocorrelation matrix. 
SIM_DO_DISP2    = 0;% Set to = 1 to see printout of the variance-covariance matrix.
SIM_DO_DISP3    = 1;% Set to = 1 to see printout of the vector of variances.
   % For the graphics output in SIM_OUT.M:
% SIM_TXT_MARKER = min(SIM_LENGTH,SIM_MAX)/10; 
                    % Where to put the labels on the graph, 
                    % a number between 1 and SIM_LENGTH.
                    % This is currently deactivated.  Instead, label
                    % positions are chosen "intelligently".
SIM_CUT = 0.6;      % This controls where the labels go on the simulation graphs.
                    % It should be a number between 0 and 1. The labels are
                    % "allowed" to "start" anywhere between 0 and SIM_CUT, if
                    % we think of the length of the entire time domain of the
                    % simulation to be normalized by unity.  E.g., for SIM_CUT = 0,
                    % the labels will always be at the left-most end, which is
                    % thus probably not a good choice.
SIM_PLOT_RAW  = 0;  % set = 1, if you want a plot of the raw, unfiltered series, 
                    % even though you have chosen to HP-filter, DO_HP_FILTER = 1.
                    % Note, that if you have chosen to save figures, then
                    % previously saved figures will be overwritten.
                    
%============> options for MOMENTS.M
DO_CALC_XX = 1;         % Set to 0, if XX is given
N_GRIDPOINTS = 64;      % gridpoints for inverse Fourier transform.  Needs to be a power of 2.
MOM_TOL =  0.0001;      % If PP has roots below TOL in absolute value,
                        % then PP may be singular, and a safer but slower route to calculate
                        % everything is taken.
if exist('HP_SELECT')~=1,
   HP_SELECT  = 1:(m_states+n_endog+k_exog); 
      % Selecting the variables for the HP Filter calculations.
else if DO_MOMENTS,
      message = ['OPTIONS.M: I will proceed with your chosen HP_SELECT.  If that leads to'
                 '           a crash in moments.m or mom_out.m or to not a complete      '
                 '           output there, then either set HP_SELECT to something        '
                 '           sensible or type clear before  running your program.        '];
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
   end;  
end;

%============> options for MOM_OUT.M
if exist('DO_HP_GRAPH')~=1, % To plot spectral densities 
    DO_HP_GRAPH=1;
end;
if (DO_HP_GRAPH & DO_MOMENTS),
   message = ['OPTIONS.M: Note, that I plot spectral densities, which you may dislike.'
              '           If you do not want to have them plotted, set DO_HP_GRAPH =0 '];
   warnings = [warnings;message];
   message = ['OPTIONS.M: For moments you can choose between joint-, sub-, and        ' 
              '           single plots. Try, SIM_JOINT=1, SIM_SUBPLOT=1, SIM_SINGLE=1 ']; 
   if DISPLAY_IMMEDIATELY, disp(message); end;
   warnings = [warnings;message];   
else
   message = ['OPTIONS.M: Note, that I will not plot spectral densities.              '
              '           If you want to have them plotted, set DO_HP_GRAPH =1        '];
   warnings = [warnings;message];
   message = ['OPTIONS.M: For moments you can choose between joint-, sub-, and        ' 
              '           single plots. Try, SIM_JOINT=1, SIM_SUBPLOT=1, SIM_SINGLE=1 ']; 
   if DISPLAY_IMMEDIATELY, 
       disp(message); 
   end;
   warnings = [warnings;message];   
end;

DO_DISP1    = 1;    % Set to = 1 to see printout of the autocorrelation matrix. 
DO_DISP2    = 0;    % Set to = 1 to see printout of the variance-covariance matrix.
DO_DISP3    = 1;    % Set to = 1 to see printout of the vector of variances.
   % For the graphics output in MOM_OUT.M:
MOM_TXT_MARKER = round(N_GRIDPOINTS/24); 
                    % Where to put the labels on the graph, 
                    % a number between 1 and N_GRIDPOINTS/2.                                    
if exist('MOM_SUBPLOT')~=1, % if spectral densities should be plotted in a subplot 
    MOM_SUBPLOT=0;
end;
if exist('MOM_JOINT')~=1, % if all spectral densities should be plotted in one figure, =0, if not 
    MOM_JOINT=1;
end; 
if exist('MOM_SINGLE')~=1, % if spectral densities should be plotted in separate figures
    MOM_SINGLE=0;
end;
if exist('MOM_SUB_FONT')~=1,  % Size of the title and labels in the subplots
    MOM_SUB_FONT=7;
end;
if exist('MOM_PLOT_RAW')~=1, % if you like non-HP-filtered spectral density plots
    MOM_PLOT_RAW=0;
end;

                    