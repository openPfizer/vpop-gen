function [p_pp,p_names,mdl_config,pp_yield,p_bnds] = ellipse_generate_pps(num_pps, vdp_model, mu, s, num_regions)
%% function ellipse_generate_pps
%   Nested ellipse specific version of generate pps
%
%   Inputs:
%       num_pps: scalar for the number of PPs to approximately generate
%       vdp_model: simbiology object for the VDP model
%       mu: multivariate-lognormal expectation of the NHANES data
%       sigma: multivariate-lognormal covariance of the NHANES data
%
%   Outputs:
%       p_pp: n x num_pps matrix of plausible patient parameters found
%
%   n is the number of parameters in the VDP model.
%
% Note: this function is method specific.
%

do_score = 1; % Flag for turning on scoring of PPs.

%% 1. Get VDP model ready for simulation

% Accelerate the model - need to setup mex compliler to do this:
sbioaccelerate(vdp_model);

% Simulate to generate baseline VP, grab names:
[~, ~, state_names] = sbiosimulate(vdp_model);

%% 2. Load free parameters
% Do not include 'parameters' redefined by SimBiology in repeated
% assignments, or parameters which are known with high accuracy.
[p_names,p_values] = get_mdl_p_to_vary(vdp_model);
num_p = numel(p_values);
% Calculate/fetch the bounds for the parameters and states:
% All state values are steady states.
search_mag = 5;
[p_bnds,state_bnds] = ...
    get_p_and_state_bnds(vdp_model,numel(state_names),p_values,search_mag);

%% 3. Get/set integration configuration
mdl_config = set_mdl_config(vdp_model);

%% 4. Prepare scoring function
[ck,p_ck] = ellipse_find_cks(num_regions,0,5); % hard-coded stop at 95% of observations.

opt_sa = saoptimset('ObjectiveLimit',0.005,'TolFun',1e-5,'Display','off',...
    'ReannealInterval',500,'InitialTemperature',0.5,'MaxIter',400,'TemperatureFcn',...
    @temperatureboltz,'AnnealingFcn', @annealingboltz,'AcceptanceFcn',@acceptancesa);

% Predefine number of PPs you want and allocate storage arrays:
p_pp = zeros(numel(p_names),num_pps);
score_v = zeros(num_pps,1);

%% 5. Vary parameters to form initial guess and perform optimization:
logp_lower  = log10(p_bnds(:,1));
logp_upper  = log10(p_bnds(:,2));
logp0       = log10(p_values);
k = 1; % counter for plausible patients
iter = 1; % iteration counter
num_pps_r = num_pps*p_ck./sum(p_ck);%num_pps/num_regions; % Expected number of PPs/region
max_iter = num_pps * 100; % maximum number of iterations (hard-coded for now)
%fprintf('Region\tPPs found\tAttempts\n');
for j = 1:num_regions
    %fprintf('Region: %d\n',j);
    % Update functions:
    score_fcn = @(a,b)ellipse_score_model(a,b,mu,s,ck(j));
    f = @(p)run_variant_and_score(p,p_names,vdp_model,mdl_config,...
        do_score,state_bnds,score_fcn);
    kr = 0; % Reset region counter
    while kr < num_pps_r(j) && iter <= max_iter
        p = p_bnds(:,1) + (p_bnds(:,2) - p_bnds(:,1)).*rand(num_p,1);
        logp = log10(p);
        logp_sa = ...
            simulannealbnd(f,logp,logp_lower,logp_upper,opt_sa);
        score_pp = f(logp_sa);
        if score_pp < 2
            p_pp(:,k) = 10.^logp_sa;
            score_v(k) = score_pp;
%             if ~mod(k,50)
%                 fprintf('%6d\t%6d\t%9d\t%8d\n',j,kr,k,iter);
%             end
            k = k + 1;
            kr = kr + 1;
        end
        iter = iter + 1; % total loops
    end
end

if iter >= max_iter
    error('pp-selection failed due to too many iterations');
end

pp_yield = (k-1)./(iter-1);

end % function ellipse_generate_pps
