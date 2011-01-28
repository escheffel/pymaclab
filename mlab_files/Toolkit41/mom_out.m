% VERSION 4.1, May 2003, COPYRIGHT H. UHLIG.
% Changes: Some improvements concerning graphics.
% Plus no "HP-filtered" in title of tables, if no HP Filtering was done.
% 
% MOM_OUT produces output from the calculations done with MOMENTS.M,
% which is assumed to have been run just before.
% This program should be modified to suit tastes and needs.  Some options
% are given in the first few lines of this program.
% It is assumed that
% VARNAMES, a matrix with (m+n+k) rows, containing the variable names, has been set.
% The program overwrites freq1, diag_select, hndl, var_index


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

if DO_HP_GRAPH,
    
%-----------------------------------MOM JOINT------------------------------------------------------------------------    
   if MOM_JOINT,
       freq1 = floor(N_GRIDPOINTS/24);
       diag_select = (n_select+1)*(0:(n_select-1)) + 1; % Selects the diagonal
          if MOM_PLOT_RAW,                  
                  hndl = plot(freqs(freq1:N_GRIDPOINTS/2), real(svv_raw(freq1:N_GRIDPOINTS/2,diag_select)));
                    for var_index = 1:n_select,
                    text(freqs(MOM_TXT_MARKER), real(svv_raw(MOM_TXT_MARKER,diag_select(var_index))),...
                    VARNAMES(HP_SELECT(var_index),:));
                    end;
              else %MOM_PLOT_RAW
      
                  hndl = plot(freqs(1:N_GRIDPOINTS/2), real(svv_fil(1:N_GRIDPOINTS/2,diag_select)));
                     for var_index = 1:n_select,
                     text(freqs(MOM_TXT_MARKER), real(svv_fil(MOM_TXT_MARKER,diag_select(var_index))),...
                     VARNAMES(HP_SELECT(var_index),:));
                     end;
           end; %MOM_PLOT_RAW    
           
        set(hndl,'LineWidth',2);
        grid;
      if MOM_PLOT_RAW,
          title('Spectral densities, unfiltered');
      else
          title('Spectral densities, Hodrick-Prescott filtered');
      end;
      xlabel('Frequency');
       
      if DO_ENLARGE,
          enlarge;
      end;
      if PRINT_FIG,
          disp('MOM_OUT: Printing plot of spectral densities');
          if DO_COLOR_PRINT,
            print -dwinc
          else
            print;
          end; 
      elseif (SAVE_FIG | SAVE_FIG_JPG),
          if SAVE_FIG,
              mom_filename=['mom_joint_',FILENAMEPIECE,'.eps'];
              disp(['MOM_OUT: Saving plot of spectral densities.  Filename is ',mom_filename]);
              if COLOR_FIG 
                  print('-depsc',mom_filename);
              else
                  print('-deps',mom_filename);
              end;      
          end;
          if SAVE_FIG_JPG,
              mom_filename=['mom_joint_',FILENAMEPIECE,'.jpg'];
              disp(['MOM_OUT: Saving plot of spectral densities.  Filename is ',mom_filename]);
              print('-djpeg90',mom_filename);    
          end;
      end;   
      disp('Inspect figure. Hit key when ready...');
      pause;
  end; %if MOM_JOINT
%----------------------------------------------END MOM JOINT------------------------------------------------------------------------


%-----------------------------------MOM SUBPLOT-------------------------------------------------------------------------------------   
   if MOM_SUBPLOT,
       
       
       freq1 = floor(N_GRIDPOINTS/24);
       diag_select = (n_select+1)*(0:(n_select-1)) + 1; % Selects the diagonal
       diag_length=length(diag_select);
       
        size_plot_x=ceil(sqrt(n_select));
        size_plot_y=round(sqrt(n_select));
       
        
           
          
       for mom_index=1:diag_length,
          if MOM_PLOT_RAW,
                  subplot(size_plot_y,size_plot_x,mom_index);                  
                  hndl = plot(freqs(freq1:N_GRIDPOINTS/2), real(svv_raw(freq1:N_GRIDPOINTS/2,diag_select(1,mom_index))));
                  text(freqs(MOM_TXT_MARKER), real(svv_raw(MOM_TXT_MARKER,diag_select(1,mom_index))),...
                  VARNAMES(HP_SELECT(mom_index),:));
            else %MOM_PLOT_RAW
                  subplot(size_plot_y,size_plot_x,mom_index);                
                  hndl = plot(freqs(1:N_GRIDPOINTS/2), real(svv_fil(1:N_GRIDPOINTS/2,diag_select(1,mom_index))));
                     text(freqs(MOM_TXT_MARKER), real(svv_fil(MOM_TXT_MARKER,diag_select(1,mom_index))),...
                     VARNAMES(HP_SELECT(mom_index),:));
             end; %MOM_PLOT_RAW    
           
         set(hndl,'LineWidth',1);
         set(gca,'FontSize',MOM_SUB_FONT);   
         %grid;
         
      if MOM_PLOT_RAW,
          title('Spectrum, unfiltered');
      else
          title('Spectrum,HP-filtered');
      end;
      xlabel('Frequency');
       
      if DO_ENLARGE,
          enlarge;
      end;      
   end; % for mom_index
   
   if PRINT_FIG,
          disp(['MOM_OUT: Printing spectral density of ',VARNAMES(HP_SELECT(mom_index),:)]);
          if DO_COLOR_PRINT,
            print -dwinc
          else
            print;
          end; 
      elseif (SAVE_FIG | SAVE_FIG_JPG),
          if SAVE_FIG,
              mom_filename=['mom_subplot',FILENAMEPIECE,'.eps'];
              disp(['MOM_OUT: Saving spectral density of. Filename is ',mom_filename]);
              if COLOR_FIG 
                  print('-depsc',mom_filename);
              else
                  print('-deps',mom_filename);
              end;      
          end;
          if SAVE_FIG_JPG,
              mom_filename=['mom_subplot',FILENAMEPIECE,'.jpg'];
              disp(['MOM_OUT: Saving spectral density of. Filename is ',mom_filename]);
              print('-djpeg90',mom_filename);    
          end;
      end;   
      
      disp('Inspect figure. Hit key when ready...');
      pause;
      close all;
  
  end; %if MOM_SUBPLOT
%----------------------------------------------END MOM SUBPLOT-------------------------------------------------------------------

   

%----------------------------------------------MOM SINGLE------------------------------------------------------------------------    
   if MOM_SINGLE,
       freq1 = floor(N_GRIDPOINTS/24);
       diag_select = (n_select+1)*(0:(n_select-1)) + 1; % Selects the diagonal
       diag_length=length(diag_select);
       for mom_index=1:diag_length,
          if MOM_PLOT_RAW,                  
                  hndl = plot(freqs(freq1:N_GRIDPOINTS/2), real(svv_raw(freq1:N_GRIDPOINTS/2,diag_select(1,mom_index))));
                    text(freqs(MOM_TXT_MARKER), real(svv_raw(MOM_TXT_MARKER,diag_select(1,mom_index))),...
                    VARNAMES(HP_SELECT(mom_index),:));
            else %MOM_PLOT_RAW      
                  hndl = plot(freqs(1:N_GRIDPOINTS/2), real(svv_fil(1:N_GRIDPOINTS/2,diag_select(1,mom_index))));
                     text(freqs(MOM_TXT_MARKER), real(svv_fil(MOM_TXT_MARKER,diag_select(1,mom_index))),...
                     VARNAMES(HP_SELECT(mom_index),:));
             end; %MOM_PLOT_RAW    
           
        set(hndl,'LineWidth',2);
        grid;
      if MOM_PLOT_RAW,
          title('Spectral densities, unfiltered');
      else
          title('Spectral densities, Hodrick-Prescott filtered');
      end;
      xlabel('Frequency');
       
      if DO_ENLARGE,
          enlarge;
      end;
      if PRINT_FIG,
          disp(['MOM_OUT: Printing spectral density of ',VARNAMES(HP_SELECT(mom_index),:)]);
          if DO_COLOR_PRINT,
            print -dwinc
          else
            print;
          end; 
      elseif (SAVE_FIG | SAVE_FIG_JPG),
          if SAVE_FIG,
              mom_filename=['mom_',VARNAMES(HP_SELECT(mom_index),:),'_',FILENAMEPIECE,'.eps'];
              disp(['MOM_OUT: Saving spectral density of ',VARNAMES(HP_SELECT(mom_index),:),'.  Filename is ',mom_filename]);
              if COLOR_FIG 
                  print('-depsc',mom_filename);
              else
                  print('-deps',mom_filename);
              end;      
          end;
          if SAVE_FIG_JPG,
              mom_filename=['mom_',VARNAMES(HP_SELECT(mom_index),:),'_',FILENAMEPIECE,'.jpg'];
              disp(['MOM_OUT: Saving spectral density of ',VARNAMES(HP_SELECT(mom_index),:),'.  Filename is ',mom_filename]);
              print('-djpeg90',mom_filename);    
          end;
      end;   
      disp('Inspect figure. Hit key when ready...');
      pause;
  end; % for mom_index
  end; %if MOM_SINGLE
%----------------------------------------------END MOM SINGLE------------------------------------------------------------------------


end; %DO_HP_GRAPH


if DO_DISP1 | DO_DISP2 | DO_DISP3,
   disp(' ');
   disp('MOM_OUT.M: Frequency-domain method based calculation of moments');
   disp('The variables are:');
   disp(VARNAMES(HP_SELECT,:));
   disp(' ');
end;
if DO_HP_FILTER,
   if DO_DISP1,
      disp('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j');
      for var_index = 1 : n_select,
         disp(sprintf('  %5.2f',autcor_fil(var_index,:)));
      end; 
      disp(sprintf('  %5.0f',autcor_fil(n_select+1,:)));
      disp(' ');
   end;
   if DO_DISP2,
      disp('Variance-Covariance Matrix, HP-filtered series:');
      for var_index = 1 : n_select,
         disp(sprintf(' %6.3f',covmat_fil(var_index,:)));
      end;
      disp(' ');
   end;
   if DO_DISP3,
      disp('Standard deviations, HP-filtered series:');
     disp(sqrt(varvec_fil));
   end;
else
   if DO_DISP1,
       disp('Autocorrelation Table, corr(v(t+j),GNP(t)).  Last row shows j');
      for var_index = 1 : n_select,
         disp(sprintf('  %5.2f',autcor_raw(var_index,:)));
      end; 
      disp(sprintf('  %5.0f',autcor_raw(n_select+1,:)));
      disp(' ');
   end;
  if DO_DISP2,
      disp('Variance-Covariance Matrix:');
      for var_index = 1 : n_select,
         disp(sprintf(' %6.3f',covmat_raw(var_index,:)));
      end;
      disp(' ');
   end;
   if DO_DISP3,
      disp('Standard deviations:');
     disp(sqrt(varvec_raw));
   end;
end;  
