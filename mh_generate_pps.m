function [p_pp,p_names,mdl_config,pp_yield,p_bnds] = mh_generate_pps(num_pps,vdp_model,mu,sigma)
%% function mh_generate_pps
%   Metropolis-Hastings algorithm-specific version of generate pps
%
%   Inputs:
%       num_pps: scalar for the number of PPs to approximately generate
%       vdp_model: simbiology object for the VDP model
%       mu: observable expectation (vector)
%       sigma: observable covariance (matrix)
%
%   Outputs:
%       p_pp: n x num_pps matrix of plausible patient parameters found
%
%   n is the number of parameters in the VDP model.
%
% Note: this function is method-specific.
%

do_score = 1; % Flag for turning on scoring of PPs.
score_fcn = @(x,y)mh_score_model(x,y,mu,sigma); % Set the scoring function for the model, this could also be a passed parameter

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
logp0 = log10(p_values);

% Calculate/fetch the bounds for the parameters and states:
% All state values are steady states.
search_mag = 5;
[p_bnds,state_bnds] = ...
    get_p_and_state_bnds(vdp_model,numel(state_names),p_values,search_mag);
logpl = log10(p_bnds(:,1));
logpu = log10(p_bnds(:,2));

%% 3. Get/set integration configuration
mdl_config = set_mdl_config(vdp_model);

%% 4. Prepare scoring function
f = @(p)run_variant_and_score(p,p_names,vdp_model,mdl_config,...
    do_score,state_bnds,score_fcn);

k = 1;
iter = 1;
max_iter = num_pps*1000;

%q(:,1)=p_lower+(p_upper-p_lower).*rand(numel(p_upper),1);
p = p_values;%p_bnds(:,1) + (p_bnds(:,2)-p_bnds(:,1)).*rand(num_p,1);
%logp = log10(p);
%logp = min(logpu,max(logpl,logp0+0.25*randn(num_p,1).*abs(logp0)));
s1 = f(log10(p));
while k <= num_pps && iter <= max_iter
    % Propose new parameter set:
    %logq = min(logpu,max(logpl,logp+0.25*randn(num_p,1).*abs(logp)));
    %q(:,i)=q(:,i-1)+(pSS(:,2)-pSS(:,1))/10.*(2*rand(m,1)-1);
    q = max(p_bnds(:,1),min(p_bnds(:,2),p + (p_bnds(:,2)-p_bnds(:,1))/10.*(2*rand(num_p,1)-1)));
    %logq = log10(q);
    s2  = f(log10(q));
    if s2 == 1e256
        s2 = eps; % correct the error checking from run_variant
    end
    r = rand();
    if r < s2/s1 && s2 > eps
%         if mod(k,100) == 0
%             fprintf('Found: %d PPs\n',k);
%         end
        %fprintf('Accepted: %3.1d over %3.1d\n',s2,s1);;
        p_pp(:,k) = q(:);
        p = q;
        %logp = log10(p);
        s1 = s2;
        k = k + 1;
    end
    iter = iter + 1;
end

if iter > max_iter
   error('Exceeded maximum iterations'); 
end

pp_yield = (k-1)/(iter-1);

end % function mh   _generate_pps
