function [t,list] = read_out_to_table(outdir)
%% function t = read_out_to_table(outdir)
% Looks for the inputed root file name, finds matching csv files, loads
% them into a datatable for manipulation and plotting.
%

if isunix
    sep_slash = '/';
else
    sep_slash = '\';
end

list = dir(strcat(outdir,'*.csv'));

if isempty(list)
   error(strcat('read_out_to_table: ',outdir,', no files found')); 
end

tcount = 0;
for i = 1:numel(list)
   if ~list(i).isdir && ~strcmp(list(i).name,'dummy.txt')
       fprintf('Adding: %s\n',list(i).name);
       
       filetmp = strcat(list(i).folder,sep_slash,list(i).name);
       
       if tcount == 0
           t = readtable(filetmp,'ReadVariableNames',true);
           tcount = tcount + 1;
           fprintf('%s\t%d\t%d\n',t.method{1},t.iter(1),t.num_pps(1));
       else
           ttmp = readtable(filetmp,'ReadVariableNames',true);
           fprintf('%s\t%d\t%d\n',ttmp.method{1},ttmp.iter(1),ttmp.num_pps(1));
           t = [t;ttmp];
           tcount = tcount + 1;
       end
   end 
end

% Filter raw table into something more useful:
t = condense_table(t,'mu');
t = condense_table(t,'sigma');
t = condense_table(t,'p_lower');
t = condense_table(t,'p_upper');
t = condense_table(t,'Pfinal');
t = condense_table(t,'VPChar');

% Scale the parameters:
t.PScaled = (t.Pfinal - t.p_lower)./(t.p_upper - t.p_lower);

if (max(max(t.PScaled)) > 1 || min(min(t.PScaled)) < 0)
    error('problem with parameter search, parameters out of bounds');
end

t.meth_cat = categorical(t.method); % categorical version comes in handy

end

function tout = condense_table(tin,str_to_find)
% Some inputs are vectors that have been split across multiple columns in the table. This function condenses them.
tout = tin;
r = contains(tin.Properties.VariableNames,str_to_find);
tout.(str_to_find) = [tin{:,r}];
tout(:,r) = [];

end