function sm = test_speed(num_pps, method, vdp_model, mu, sigma, num_regions)
%% function test_speed.m
%   Evaluates the performance of 'GenerateVPs', and prevelance_select_func and
%   generates a few key metrics as a output struct including:
%
%   Inputs:
%       num_pps - desired number of plausible patients (scalar)
%       method - selector for solution method
%           valid inputs for method: 'simanneal','ellipse','ga'
%   Outputs:
%       sm structure:
%       sm.gof - goodness-of-fit to the histograms
%       sm.selection_efficiency - how efficient was plausible patients --> VPs
%       sm.time_per_pp - how much time did it take to generate all of the PPs
%       sm.time_per_vp - time_per_pp * selection_efficiency
%       sm.cputime_pp - total cputime for generating pp
%       sm.cputime_total - total cputime for calling TestSpeed()
%       sm.n_vps - number of vps generated
%       sm.n_pps - number of pps generated, just for convenience
%

%% 1. Choose the method for finding plausible patients:
switch method
    case 'SA'
        gen_pp_f = @(x)simanneal_generate_pps(x,vdp_model);
    case 'GA'
        gen_pp_f = @(x)ga_generate_pps(x,vdp_model,mu,sigma);
    case 'MH'
        gen_pp_f = @(x)mh_generate_pps(x,vdp_model,mu,sigma);
    case 'NSA'
        gen_pp_f = @(x)ellipse_generate_pps(x,vdp_model,mu,sigma, num_regions);
    otherwise
        fprintf('Valid methods are: ''SA'',''GA'', ''MH'', and ''NSA''.\n');
        fprintf('Method requested: %s.\n',method);
        error('Requested method not implemented, yet');
end

%% 2. Generate plausible patient cohort:
% Start the timing clock
cpu_time_start = cputime;
pp_start_time1 = cputime;
[p_pp,p_names,mdl_config,pp_yield,p_bnds] = gen_pp_f(num_pps);
pp_time = cputime - pp_start_time1;

%% 3. Simulate the new plausible patients to steady state:
for ipp = 1:size(p_pp,2)
    logp = log10(p_pp(:,ipp));
    [~,xpp(:,ipp)] = run_variant_and_score(logp,p_names,vdp_model,mdl_config,0);
end

% Create NHANES outputs:
pp_obs(:,1) = xpp(5,:)' + xpp(6,:)';                % HDL
pp_obs(:,2) = xpp(9,:)';                            % LDL
pp_obs(:,3) = xpp(7,:)' + xpp(5,:)' + xpp(6,:)';    % TC

% Unit conversion to match NHANES
pp_obs = pp_obs*38.66/2.79; % [mg/dl], unit conversion from [mM]

%% 4. Select the virtual population based on the NHANES data distribution
[select,gof,~,~,p_incld,sf] = get_prevalence(mu,sigma,log(pp_obs),p_pp);

%% 5. Get key metrics on performance
total_cpu_time = cputime - cpu_time_start;
%How well does your histograms match the data? (lower the better)
sm.gof = gof;
%How many plausible patients do you need to make 1 VP? (lower the better)
sm.selection_efficiency = numel(select)/sum(select);
%How long to make 1 plausible patient (shorter the better)
sm.time_per_pp = (pp_time)/num_pps;
%How long to make 1 VP? (shorter the better)
sm.time_per_vp = sm.time_per_pp*sm.selection_efficiency;
%How long is the cpu time for generating pps?
sm.cputime_pp = pp_time;
%How long is the cpu time for all this procedure?
sm.cputime_total = total_cpu_time;
%How many virtual patients you generated in the end?
sm.n_vps = sum(select);
%Record the number of PPs generated.
sm.n_pps = num_pps;

%Store the generated characteristics for selection
sm.VPChar      = pp_obs ;
%Store the selected population
sm.selection   = select;
%Store the probability of inclusion
sm.ProbInclude = p_incld;
%Store the final parameters for the plausible population
sm.Pfinal      = p_pp;
%Store the scale-factore to multiply p_include
sm.betaS       = sf;
%Mean of MVN distribution for fitting
sm.mu          = mu;
%Covariance matrix of MVN distribution for fitting
sm.sigma       = sigma;
%Simbiology model
sm.VPmodel     = vdp_model;
%Method Name
sm.method_name = method;
% Limits on parameter values:
sm.p_bnds = p_bnds;

%fprintf('\nPP yield = %3.2f\n',pp_yield);

end