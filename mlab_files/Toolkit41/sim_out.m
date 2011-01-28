% VERSION 4.0, November 2002, COPYRIGHT H. UHLIG.
% SIM_OUT.M produces output for SIMUL.M.  It is controlled 
% by several options to be described
% below.  All these options are set in OPTIONS.M,
% see that file for further details.
%
% The series will be plotted
% if SIM_GRAPH is set to 1.  In that case, further
% modifications can be done in particular with:
%   DO_ENLARGE : = 1, if you want large font sizes for the text on your plots.
%                     Good for slides.
%   SIM_JOINT  : = 1, if you want all series on the same graph, else = 0.
%   PRINT_FIG  : = 1, if you want plots to be printed on your printer
%   SAVE_FIG   : = 1, if you want plots to be saved as encaps. postscript. 
%   SAVE_FIG_JPG : = 1, if you want plots to be saved as jpg-files.
%                     Set PRINT_FIG = 0 also. The filenames are sim_ser1.eps, ...
%                     if SIM_JOINT = 0, and sim_data.eps is SIM_JOINT = 1.
%   SIM_PLOT_RAW : = 1, if you want a plot of the raw, unfiltered series, 
%                       even though you have chosen to HP-filter, DO_HP_FILTER = 1.                       
%                       Note, that if you have chosen to save figures, then
%                       the previously saved figures will be overwritten.
%                       This option is useful, if you want to look at the plot
%                       of the raw simulations, after having already seen the filtered
%                       ones: simply type
%                       SIM_PLOT_RAW = 1;
%                       sim_out;
%
% For printing the numbers of the autocorrelation table, 
% the following options are helpful:
%   
%  SIM_DO_DISP1: Set to = 1 to see printout of the autocorrelation matrix. 
%  SIM_DO_DISP2: Set to = 1 to see printout of the variance-covariance matrix.
%  SIM_DO_DISP3: Set to = 1 to see printout of the vector of variances.


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

n_select = max(size(SIM_SELECT));
if SIM_GRAPH,
   time_axis = (0:(SIM_LENGTH-1))/PERIOD + SIM_DATE0;
   
%---------------------------------------------JOINT PLOTS-----------------------------------------------------------------
   if SIM_JOINT
      if SIM_PLOT_RAW,
         hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
            sim_raw(SIM_SELECT,1:min(SIM_LENGTH,SIM_MAX)));
         [mx,pos]=max(abs(sim_raw(SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX)))'));
         for var_index = 1:n_select,
%            text(time_axis(SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
%              VARNAMES(SIM_SELECT(var_index),:));
            text(time_axis(pos(var_index)), sim_raw(var_index,pos(var_index)),...
              VARNAMES(SIM_SELECT(var_index),:));
         end;
      else
         hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
            sim_xyz(SIM_SELECT,1:min(SIM_LENGTH,SIM_MAX)));
         [mx,pos]=max(abs(sim_xyz(SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX)))'));
         for var_index = 1:n_select,
%            text(time_axis(SIM_TXT_MARKER), sim_xyz(var_index,SIM_TXT_MARKER),...
%              VARNAMES(SIM_SELECT(var_index),:));
            text(time_axis(pos(var_index)), sim_xyz(var_index,pos(var_index)),...
              VARNAMES(SIM_SELECT(var_index),:));
         end;
      end;
      set(hndl,'LineWidth',2);
      grid;
      if DO_HP_FILTER & ~SIM_PLOT_RAW,
          title('Simulated data (HP-filtered)');
      else
          title('Simulated data');
      end;
      xlabel('Year');
      ylabel('Percent deviation from steady state');
      if DO_ENLARGE,
          enlarge;
      end;
      if PRINT_FIG,
          disp('SIM_OUT: Printing plot of Simulated data');
          print;
      elseif (SAVE_FIG | SAVE_FIG_JPG),
          if SAVE_FIG,
              sim_filename=['sim_data_joint_',FILENAMEPIECE,'.eps'];
              disp(['SIM_OUT: Saving plot of Simulated data.  Filename is ',sim_filename]);
              if COLOR_FIG 
                  print('-depsc',sim_filename);
              else
                  print('-deps',sim_filename);
              end;      
          end;
          if SAVE_FIG_JPG,
              sim_filename=['sim_data_joint_',FILENAMEPIECE,'.jpg'];
              disp(['SIM_OUT: Saving plot of Simulated data.  Filename is ',sim_filename]);
              print('-djpeg90',sim_filename);    
          end;
      end;                    
          disp('Inspect figure. Hit key when ready...');
          pause;
      end; % if SIM_JOINT; 
%--------------------------------------------END JOINT PLOTS-------------------------------------------------------------

%--------------------------------------------SIM SUBPLOT-----------------------------------------------------------------
if SIM_SUBPLOT, 
        size_plot_x=ceil(sqrt(n_select));
        size_plot_y=round(sqrt(n_select));
      for var_index = 1:n_select,
          subplot(size_plot_y,size_plot_x,var_index);                   
          
         if SIM_PLOT_RAW,
            hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
               sim_raw(SIM_SELECT(var_index),1:min(SIM_LENGTH,SIM_MAX)));
            % text(time_axis(SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
            %    VARNAMES(SIM_SELECT(var_index),:));
         else
            hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
               sim_xyz(SIM_SELECT(var_index),1:min(SIM_LENGTH,SIM_MAX)));
            % text(time_axis(SIM_TXT_MARKER), sim_xyz(var_index,SIM_TXT_MARKER),...
            %    VARNAMES(SIM_SELECT(var_index),:));
         end;
         set(hndl,'LineWidth',1);
         set(gca,'FontSize',SIM_SUB_FONT);   
         %grid;
         if DO_HP_FILTER & ~SIM_PLOT_RAW,
            title(['Sim. data(HP-filt.): ',VARNAMES(SIM_SELECT(var_index),:)]);
         else
            title(['Sim. data: ',VARNAMES(SIM_SELECT(var_index),:)]);
         end;
         xlabel('Year');
         ylabel('Percent dev. from SS');
         if DO_ENLARGE,
            enlarge;
         end;         
          
      end; %for varindex
      if PRINT_FIG,
             disp(['SIM_OUT.M: Printing plot of simulated data']);
             print;
         elseif ( SAVE_FIG | SAVE_FIG_JPG ),
             if SAVE_FIG,
                 sim_filename=['sim_subplot_',FILENAMEPIECE,'.eps'];
                 disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
                 disp(['         as encapsulated postscript file. Filename is ',sim_filename]);
                 if COLOR_FIG 
                     print('-depsc',sim_filename);
                 else
                     print('-deps',sim_filename);
                 end;
             end;
             if SAVE_FIG_JPG,
                 sim_filename=['sim_subplot_',FILENAMEPIECE,'.jpg'];
                 disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
                 disp(['         as jpg-file. Filename is ',sim_filename]);
                 print('-djpeg90',sim_filename);
             end;
         end;  
         disp('Inspect figure. Hit key when ready...');
         pause;
         close all;
      
  end; %SIM_SUBPLOT

%---------------------------------------------END SIM SUBPLOT------------------------------------------------------------

%---------------------------------------------SINGLE PLOTS---------------------------------------------------------------
      
   if SIM_SINGLE,    
      for var_index = 1:n_select,
         if SIM_PLOT_RAW,
            hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
               sim_raw(SIM_SELECT(var_index),1:min(SIM_LENGTH,SIM_MAX)));
            % text(time_axis(SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
            %    VARNAMES(SIM_SELECT(var_index),:));
         else
            hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
               sim_xyz(SIM_SELECT(var_index),1:min(SIM_LENGTH,SIM_MAX)));
            % text(time_axis(SIM_TXT_MARKER), sim_xyz(var_index,SIM_TXT_MARKER),...
            %    VARNAMES(SIM_SELECT(var_index),:));
         end;
         set(hndl,'LineWidth',2);
         grid;
         if DO_HP_FILTER & ~SIM_PLOT_RAW,
            title(['Simulated data (HP-filtered): ',VARNAMES(SIM_SELECT(var_index),:)]);
         else
            title(['Simulated data: ',VARNAMES(SIM_SELECT(var_index),:)]);
         end;
         xlabel('Year');
         ylabel('Percent deviation from steady state');
         if DO_ENLARGE,
            enlarge;
         end;
         if PRINT_FIG,
             disp(['SIM_OUT.M: Printing simulation of ',VARNAMES(var_index,:),'...']);
             print;
         elseif ( SAVE_FIG | SAVE_FIG_JPG ),
             if SAVE_FIG,
                 sim_filename=['sim_ser_',VARNAMES(var_index,:),FILENAMEPIECE,'.eps'];
                 disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
                 disp(['         as encapsulated postscript file. Filename is ',sim_filename]);
                 if COLOR_FIG 
                     print('-depsc',sim_filename);
                 else
                     print('-deps',sim_filename);
                 end;
             end;
             if SAVE_FIG_JPG,
                 sim_filename=['sim_ser_',VARNAMES(var_index,:),FILENAMEPIECE,'.jpg'];
                 disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
                 disp(['         as jpg-file. Filename is ',sim_filename]);
                 print('-djpeg90',sim_filename);
             end;
         end;                    
            disp('Inspect figure. Hit key when ready...');
            pause;
        end;
   end;
   
   %--------------------------------------------------END SIM SINGLE----------------------------------------------------------
end; %If SIM_GRAPH




disp(' ');
if SIM_DO_DISP1 | SIM_DO_DISP2 | SIM_DO_DISP3,
   disp(        'SIM_OUT.M: Simulation-based calculation of moments');
   disp(sprintf('           Simulation length = %10d',SIM_LENGTH));
   if SIM_MODE == 2,
      disp(sprintf('           repeated %10d times.',SIM_N_SERIES));
   end;      
   disp('The variables are:');
   disp(VARNAMES(SIM_SELECT,:));
   disp(' ');
end;
if SIM_DO_DISP1,
   if DO_HP_FILTER,
      disp('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j');
      disp('(Simulation-based calculations)');
   else
      disp('Autocorrelation Table, corr(v(t+j),GNP(t)).  Last row shows j');
      disp('(Simulation-based calculations)');
   end;
   for var_index = 1 : n_select,
      disp(sprintf('  %5.2f',autcor_sim(var_index,:)));
   end; 
   disp(sprintf('  %5.0f',autcor_sim(n_select+1,:)));
   disp(' ');
   if (SIM_N_SERIES > 3) & (SIM_MODE == 2),
      disp('Small sample standard errors for the Autocorrelation Table:');
      disp('(Simulation-based calculations)');
      for var_index = 1 : n_select,
         disp(sprintf('  %5.2f',autcor_std(var_index,:)));
      end; 
      disp(sprintf('  %5.0f',autcor_sim(n_select+1,:)));
      disp(' ');
   end;
end;
if SIM_DO_DISP2,
   if DO_HP_FILTER,
      disp('Variance-Covariance Matrix (HP-filtered series):');
      disp('(Simulation-based calculations)');
   else
      disp('Variance-Covariance Matrix:');
      disp('(Simulation-based calculations)');
   end;
   for var_index = 1 : n_select,
      disp(sprintf(' %6.3f',covmat_sim(var_index,:)));
   end;
   disp(' ');
   if (SIM_N_SERIES > 3) & (SIM_MODE == 2),
      disp('Small sample standard errors for the Variance-Covariance Matrix:');
      disp('(Simulation-based calculations)');
      for var_index = 1 : n_select,
         disp(sprintf(' %6.3f',covmat_std(var_index,:)));
      end;
      disp(' ');
   end;
end;
if SIM_DO_DISP3,
   if DO_HP_FILTER
      disp('Standard deviations (HP-filtered series):');
      disp('(Simulation-based calculations)');
   else
      disp('Standard deviations:');
      disp('(Simulation-based calculations)');
   end;
   disp(stdvec_sim);
   if (SIM_N_SERIES > 3) & (SIM_MODE == 2),
      disp('Small sample standard errors for the Standard deviations:');
      disp('(Simulation-based calculations)');
      disp(stdvec_std);
   end;
end;
