function [mvmodel, r] = correlate_nhanes_chol(do_plots)
%clear;clc;close all;rng('shuffle');
if nargin == 0
    do_plots = 0;
end

% Load the NHANES data:
thdl = xptread('HDL_G.XPT');
ttc = xptread('TCHOL_G.XPT');
ttrig = xptread('TRIGLY_G.XPT');

% Merge based on the key SEQN:
c1 = join(thdl,ttc,'Keys','SEQN');
c1 = join(ttrig,c1,'Keys','SEQN');

% Log-transform the relevant columns:
lnHDL = log(c1.LBDHDD);
lnLDL = log(c1.LBDLDL);
lnTC = log(c1.LBXTC);

if do_plots
    % Output some correlations and univariate distributions (justifies
    % log-transformation):
    figure('Name','TC vs. HDL');
    scatterhist(lnHDL,lnTC);
    xlabel('log(HDL)');
    ylabel('log(TC)');
    figure('Name','TC vs. LDL');
    scatterhist(lnLDL,lnTC);
    xlabel('log(LDL)');
    ylabel('log(TC)');
    figure('Name','HDL vs. LDL');
    scatterhist(lnLDL,lnHDL);
    xlabel('log(LDL)');
    ylabel('log(HDL)');
end

%% Find the joint probability distribtion:
r = [lnHDL lnLDL lnTC]; % Gather the variables of interest
r = r(~any(isnan(r'))',:); % Clear any rows with NaN, not absolutely necessary, but avoids a warning message

% Fit a multivariate Gaussian model to the log-transformed data:
mvmodel = fitgmdist(r,1); % model parameters should match mean(r) and cov(r).

if do_plots
    fprintf('Joint Gaussian Fit\nlog(HDL)\tlog(LDL)\tlog(TC)\nmu:\n');
    disp(mvmodel.mu);
    fprintf('covariance matrix:\n');
    disp(mvmodel.Sigma);
end
% Clear the raw file data:
%clear thdl ttc ttrig
