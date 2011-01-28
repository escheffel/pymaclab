% VERSION 4.1, MAY 2003, COPYRIGHT H. UHLIG.
% README.M contains hints on how to use the programs in this directory.
% Start MATLAB and type
% readme
% on a single line.  Alternatively, inspect this file with any
% text editor.


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

disp('The files here are VERSION 4.1 of my "toolkit" programs."');
disp('TO SEE WHAT IS NEW HERE, TYPE IN');
disp('whatsnew');
disp('The files in this directory perform the calculations described');
disp('in Harald Uhlig (1995), "A Toolkit for Analyzing Nonlinear Dynamic');
disp('Stochastic Models Easily", Discussion Paper, Institute for');
disp('Empirical Macroeconomics, Federal Reserve Bank of Minneapolis.');
disp('To see quickly, how they work, start MATLAB and type');
disp('example1');
disp('to calculate through example 1, i.e. Hansen RBC model.');
disp('There are other examples: type');
disp('help exampl1.m');
disp('to get help on any of these examples. Use the example files');
disp('as templates for your own work.  Alternatively, declare all');
disp('needed matrices and type in');
disp('do_it');
disp('to do all calculations.');
disp(' ');
disp('Aside from the example files, the files in this directory are:');
disp('do_it.m      : does it all, once all needed matrices are defined.');
disp('options.m    : sets the options for all programs.  It is called');
disp('               by do_it and needs to be called, if any of the');
disp('               following routines is used in isolation.');
disp('impresp.m    : calculates and shows impulse responses to shocks.');
disp('enlarge.m    : allows you to manipulate letter sizes on plots');
disp('mom_out.m    : produces output.  To be called after moments.m');
disp('moments.m    : calculates second moment properties.');
disp('readme.m     : This file.  It tells you what to do.');
disp('sol_out.m    : produces output.  To be called after solve.m');
disp('solve.m      : solves for the coefficent matrix PP.');
disp('qzswitch.m   : necessary for QZ decomposition.');
disp('qzdiv.m      : necessary for QZ decomposition.');
disp('calc_qrs.m   : calculates the other coefficient matricies for the recursive law of motion');
disp(' ');
disp('All files are extensively documented.  Type');
disp('help filename');
disp('in MATLAB to get more information. Note that these files');
disp('set some additional variables, which you may have used before:');
disp('thus, be careful not to use names appearing in the programs.');
disp('If you discover a serious mistake, please contact me at uhlig@wiwi.hu-berlin.de.');
disp('If you have a question, please do NOT contact me.  Please read the paper.');
disp('Feel free to copy, modify and use these files at your own risk.');
disp('There is absolutely no guarantee that this stuff works the way it is');
disp('supposed to.  Have fun.');
disp('Again: The files here are VERSION 4.1 of my "toolkit" programs."');
disp('TO SEE WHAT IS NEW HERE, TYPE IN');
disp('whatsnew');
disp('                Harald Uhlig, Berlin, May 2003');
