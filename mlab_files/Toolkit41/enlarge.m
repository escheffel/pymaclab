% VERSION 2.0, MAY 2001, COPYRIGHT H. UHLIG.
% ENLARGE.M modifies an existing plot with a graphics handle
% called hndl by enlarging the character sizes.  This is useful
% for projecting the graphs via an overhead projector, for example.
% It is particularly useful, when creating encapsulated postscript graphics
% with 
% print -deps figure.eps
% for example.
%
% ENLARGE must be preceded by a command of the type 
% hndl = plot( ... );
% and, if desired, the title, etc. commands.
% 
% The size of the letters are controlled by LET_SIZE for titles
% etc. and LAB_SIZE for labels.  If you want to change them BEFORE you have
% run enlarge.m for a particular graph, just set LET_SIZE and LAB_SIZE
% as you please.  Otherwise, the program chooses some sensible default values.
% There are various, obvious ways to further fine-tune this by changing the program,
% if you so desire.

if exist('LET_SIZE')~=1,
    LET_SIZE    = 18;  
end; 
if exist('LAB_SIZE')~=1,
    LAB_SIZE    = 24;  
end; 

if exist('DO_THICK_LINES')~=1,
    DO_THICK_LINES = 1;
end;

par_hndl = get(hndl(1),'Parent');
old_size = get(par_hndl,'FontSize');
txt_hndl = findobj(par_hndl,'FontSize',old_size);
for hndl_j = 1 : max(size(txt_hndl)),
   set(txt_hndl(hndl_j),'FontSize',LET_SIZE);
end;

if DO_THICK_LINES,
    for hndl_j = 1 : max(size(hndl)),
        hndl_struc = get(hndl(hndl_j));
        if isfield(hndl_struc,'LineWidth'),
            set(hndl(hndl_j),'LineWidth',3);
        end;
    end;
end;

par_struc = get(par_hndl);

if isfield(par_struc,'Title'),
    tit_hndl = get(par_hndl,'Title');
    old_size = get(tit_hndl,'FontSize');
    txt_hndl = findobj(tit_hndl,'FontSize',old_size);
    for hndl_j = 1 : max(size(txt_hndl)),
        set(txt_hndl(hndl_j),'FontSize',LET_SIZE);
    end;
end;

if isfield(par_struc,'XLabel'),
    x_hndl = get(par_hndl,'XLabel');
    old_size = get(x_hndl,'FontSize');
    txt_hndl = findobj(x_hndl,'FontSize',old_size);
    for hndl_j= 1 : max(size(txt_hndl)),
        set(txt_hndl(hndl_j),'FontSize',LET_SIZE);
    end;
end;

if isfield(par_struc,'YLabel'),
    y_hndl = get(par_hndl,'YLabel');
    old_size = get(y_hndl,'FontSize');
    txt_hndl = findobj(y_hndl,'FontSize',old_size);
    for hndl_j = 1 : max(size(txt_hndl)),
        set(txt_hndl(hndl_j),'FontSize',LET_SIZE);
    end;
end;


if isfield(par_struc,'ZLabel'),
    z_hndl = get(par_hndl,'ZLabel');
    old_size = get(z_hndl,'FontSize');
    txt_hndl = findobj(z_hndl,'FontSize',old_size);
    for hndl_j = 1 : max(size(txt_hndl)),
        set(txt_hndl(hndl_j),'FontSize',LET_SIZE);
    end;
end;

if isfield(par_struc,'Children'),
    chil = get(par_hndl,'Children');
    for hndl_j = 1 : max(size(chil)),
        chil_struc = get(chil(hndl_j));
        if isfield(chil_struc,'FontSize'),
            old_size = get(chil(hndl_j),'FontSize');
            label_chil = findobj(chil(hndl_j),'FontSize',old_size);
            set(label_chil,'FontSize',LAB_SIZE);
        end;
    end;
end;
