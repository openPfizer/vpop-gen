function a_make_paper_figs(load_file_name)
%% function a_make_paper_figs
%   Creates figures from the manuscript. Only input is the full path and
%   root filename to the .csv files created by a_run_vpop_fit.m. For
%   example, if a_run_vpop_fit wrote myrun_xxx.csv to /Documents/me/files/
%   the input would be '/Documents/me/files/myrun_' (the numbers and .csv
%   are added later).

close all;rng('shuffle');
method_order = {'SA';'NSA';'MH';'GA'}; % order of the methods for plotting

%% LOAD FILE NAME
if nargin ~= 1
    error('a_make_paper_figs: requires a filestring input.');
end

%% Read data files to Matlab table:
t = read_out_to_table(load_file_name); 

%% Plots for each individual method:
plot_histograms(t,method_order);
plot_violin(t,method_order); % parameters are scaled within function

%% Summary plots across methods:
plot_gof_vp_vs_pp(t,method_order); % Comparison of plausible population --> virtual population
plot_orthog(t,method_order); % plot the orthogonality of the Vpops
plot_time_summary_bar(t,method_order); % plot time to pp/vps
plot_gof_vs_npps(t,method_order);

end
