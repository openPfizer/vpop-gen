function [select,hist_score,hist_mu,hist_std,p_incl,sf] = ...
    get_prevalence(mu,sigma,pps_obs,p_pp,sf_in,just_score)
%% function get_prevalence(mu,sigma,pps_obs)
%   Main selection function for converting plausible patients into virtual
%   patients. Requires inputs:
%       mu          - m x 1 vector, log-normal parameter of actual distribution
%       sigma       - m x 1 vector, log-normal parameter of actual distribution
%       pps_obs     - n x m matrix, plausible patient values where
%               n is the number of plausible patients (1 row per patient)
%               m is the number of observations to be matched (same size as
%               mu).
%
%% Check inputs:
find_sf = 1;

if nargin >= 5
    find_sf = 0; % no need to find the scaling factor
end

do_select = 1;
if nargin == 6
    do_select = ~just_score;
end

%% 1. Setup:
num_pps = numel(pps_obs(:,1)); % number of plausible patients
num_obs = numel(mu); % number of distributions to match (1d histograms)
mu1     = repmat(mu,num_pps,1); % for use in scaling VPs

%% 2. Generate Test Data:
% Based on mu and sigma, generate multi-variate PDF (assume Gaussian).
% For real data will need individual data, or infer distribution using known
% distr. and cross-correlations. Assumining multi-variate normal for now.

% Data input check:
sig_pps    = repmat(diag(sigma)',num_pps,1);
pps_obs_scal = (pps_obs-mu1)./sig_pps.^0.5;                % Normalize VPs
data_pdf   = mvnpdf(pps_obs_scal.*sig_pps.^0.5+mu1,mu,sigma); % Data density

%% 3. Calculate the prob. of inclusion for each plausible patient:
% Probability density (PDF) at the VP's location divided by the density of
% plausible patients in the same region. The plausible patient density is
% created by taking the volume of the sphere to the num_pts nearest
% neighbors (default = 5).
num_pts     = 5; % Nearest neighbor points for density calculation

try
    vp_tree     = KDTreeSearcher(pps_obs_scal,'BucketSize',1);
catch
    save('results/fail.mat','pps_obs','pps_obs_scal','p_pp');
    error('failed at KDTreeSearcher');
end

[~,dpps]    = knnsearch(vp_tree,pps_obs_scal,'K',num_pts);  % Assumes non-identical VPs
ball_vol    = get_nsphere_vol(num_obs,max(dpps'));
vp_density  = num_pts./ball_vol;        % PP density based on num_pts nearest neighbors
p_incl      = (data_pdf)./vp_density';  % Uniform probability of inclusion

%% 4. Optimize the scaling factor for selection of VPs:
runs    = 10; % Number of times to try the fitting
fopt    = @(p)get_vp_select(p,p_incl,pps_obs_scal,runs,do_select);
sf_max  = 1/max(p_incl); % Maximum scaling factor
sf_lower = 0; % lower bound for the scaling factor on probability of inclusion
sf_upper = log10(1000*sf_max); % upper bound for the scaling factor on the probability of inclusion

options_sa = saoptimset('TolFun',1e-15,'ReannealInterval',50000,...
    'InitialTemperature',0.5,'MaxIter',1000,'Display','off',...
    'TemperatureFcn',@temperatureboltz,'AnnealingFcn',@annealingboltz,...
    'AcceptanceFcn',@acceptancesa);
%   options=saoptimset('TolFun',1e-15,'ReannealInterval',50000,'InitialTemperature',1,'MaxIter',1000,'TemperatureFcn',@temperatureboltz,'AnnealingFcn', @annealingboltz,'AcceptanceFcn',@acceptancesa);

if find_sf
    k   = simulannealbnd(fopt,log10(sf_max),sf_lower,sf_upper,options_sa); % optimize the log10(scaling factor) value
    sf  = 10.^k(1); % scaling factor transformed back from log10
else
    k = log10(sf_in);
    sf = sf_in;
end

%% 5. With the Optimal Scaling Factor, [Re]select VPs:
% Iterate several times (default = 100) and select the best fit and mean
% fits.
num_tests = 100;
hist_score_temp = zeros(num_tests,1);
selected_temp   = zeros(num_pps,num_tests);
for i=1:num_tests
    [hist_score_temp(i),selected_temp(:,i)] = get_vp_select(k,p_incl,pps_obs_scal,1,do_select);
end
[hist_score,I]  = min(hist_score_temp);
select          = logical(selected_temp(:,I));
hist_mu         = nanmean(hist_score_temp);
hist_std        = nanstd(hist_score_temp);

end % function get_prevalence

% ************************************************************************
% Helper functions:

function vol = get_nsphere_vol(n,radius)
%% function get_nsphere_vol
% Calculates the volume of a n-dimensional sphere
% formula can be found on Wikipedia.

vol = pi^(n/2).*(1/gamma(n/2+1)).*radius.^n;

end

% ************************************************************************
function [hist_score,select,p_out] = ...
    get_vp_select(sf_log10, p_incld, pps, runs,do_select)
%% function get_vp_select(sf_log10, p_incld, pps, runs)
%   Selects virtual patients from the cohort of plausible patients
%   previously generated. And scores them for normality using a series of
%   K-S tests.
%
%   Inputs:
%       sf_log10 - 1x1 scalar, log10(scaling factor) for prob. of inclusion
%       p_incld - nx1 vector, probability of inclusion for each plausible
%           patient (pre-scaling-factor).
%       pps - nxm matrix, plausible patients for inclusion
%       runs - 1x1 scalar, number of times to attempt the fitting
%
%   Where n is the number of plausible patients amd m is the number of
%   observations being matched to data.
%
%   Outputs:
%       hist_score  - 1x1 scalar, sum of all univariate K-S tests across
%           all runs.
%       select      - nx1 boolean,
%       p_out       - 1x1 scalar, p-value for K-S tests
%

%% 1. Setup
sf = 10.^sf_log10(1); % scale factor input is log10
data_dim = numel(pps(1,:));

% Pre-allocate output vectors:
p_val        = zeros(data_dim,1);
ksstat      = zeros(data_dim,1);
cv          = zeros(data_dim,1);
hist_score  = zeros(runs,1);
num_pps     = size(pps,1);

if ~do_select
    p_incld = ones(num_pps,1)/sf;
end

%% 2. Perform Selection and Univariate K-S Tests:
for j = 1:runs
    r = rand(num_pps,1);
    select = r < p_incld*sf;
    num_vp = sum(select);
    if num_vp > 1
        for i = 1:data_dim
            %Test if dist. is standard normal
            [~,p_val(i),ksstat(i),cv(i)] = kstest(pps(select,i));
        end
        hist_score(j) = sum(ksstat);
        p_out = p_val;
    else
        hist_score(j) = 1*data_dim;
    end
end
hist_score = sum(hist_score)/runs;

end % function get_vp_select
