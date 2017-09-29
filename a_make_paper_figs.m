%% Startup
clear;clc;close all;rng('shuffle');

%% LOAD FILE NAME
load_file_name = 'txtout/github_run1_';
method_order = {'SA';'NSA';'MH';'GA'}; % order of the methods for plotting

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

%% EoS
