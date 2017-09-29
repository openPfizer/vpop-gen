%% a_run_vpop_fit -- main simulation script
% This script generates the data from the manuscript by looping through a
% user-specfied number of iterations for any number of target PPs for the
% specified methods. The default options are the ones used in generation of
% the manuscript figures (results will differ slightly due to
% stochasticity). See the Key Inputs section for the user-specified inputs
% including:
%   resume_run - boolean (0 or 1) for starting a new run (0) or resuming a
%       previous run (1).
%   txt_name - string for the root of the file names, this will be used for
%       writing and reading back results as well as determinig a resume
%       point.
%   num_pps - vector of number of plausible patients to target
%   num_iters - number of iterations for each method and num_pp entry. Use
%       this option to gather better statistics.
%   num_regions - NSA method only option for setting the number of
%       ellipses.
%   methods - cell array specifying which methods to run. Accepted options
%       are 'SA', 'NSA', 'MH', and 'GA'.
%   mdl_mat - location of the exported SimBiology model file
%

%% Startup
clear;clc;close all;rng('shuffle');
addpath('./NHANES');
addpath('./model');
addpath('./allcomb');

% Check if output directory exists (do not change this name without also 
% correcting ga_interim.m and ga_generate_pps.m):
if exist('txtout','dir') ~= 7
    mkdir('txtout');
end

%% Locate save files and if this is a new run:
resume_run = 0; % <--- SETS NEW RUN OR RESUME (0 = NEW RUN, 1 = RESUME BASED ON TXTFILES)
txt_name = 'txtout/github_run2_'; % Running text files for extremely long runs or basis for resume

%% Key inputs, ONLY relevant if resume_run == 0 (otherwise overwritten):
num_pps     = [100;500;1000;5000;10000]; % Number of plausible patients to attempt to create
num_iters   = 1;                        % <--- Number of iterations
num_regions = 5;                        % NSA-method only input, un-needed otherwise
methods      = {'SA';'NSA';'MH';'GA'};  % Method to use, must be a precise input (e.g., 'NSA','SA','GA')
mdl_mat     = 'My_Model.mat';           % Location of the exported SimBiology model

%% Find VPs de novo, or load last results:
if ~resume_run
    i_start = 1;
    
    % NHANES Log-Normal Fit - see correlate_NHANES_chol.m
    m       = correlate_nhanes_chol(0); % fetch the NHANES data
    mu      = m.mu;                     % log-normal distribution parameters
    sigma   = m.Sigma;                  % log-normal distribution parameters
    
    % Load the model
    load(mdl_mat);
    num_pps = sort(num_pps); % just in case
    atmp = allcomb(1:numel(methods),1:num_iters,1:numel(num_pps)); % define all the cases to be run
    i_end = size(atmp,1);
    save(strcat(txt_name,'_run_setup.mat'));
else
    load(strcat(txt_name,'_run_setup.mat'));
    list = dir(strcat(txt_name,'*.csv'));
    i_start = numel(list)+1;
end

%% Main loop, go through for each case (method/#pps), generate results
for i_all_iter = i_start:i_end
    fprintf('%d of %d: %s, %d\n',i_all_iter,size(atmp,1),...
        methods{atmp(i_all_iter,1)},num_pps(atmp(i_all_iter,3))); % Update on progress
    % Call test_speed the master function:
    j(i_all_iter).sm = test_speed(num_pps(atmp(i_all_iter,3)), ...
        methods{atmp(i_all_iter,1)}, van_de_pas_mod1, mu, sigma, num_regions);
    j(i_all_iter).num_pps = num_pps(atmp(i_all_iter,3)); % nominal number of pps
    j(i_all_iter).method = methods{atmp(i_all_iter,1)};
    write_j_to_file(txt_name,j(i_all_iter),i_all_iter); % write the results to a textfile
end

%% EoS