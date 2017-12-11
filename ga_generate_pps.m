function [p_pp,p_names,mdl_config,pp_yield,p_bnds] = ga_generate_pps(num_pps,vdp_model,mu,sigma)
%% function ga_generate_pps
%   Genetic algorithm-specific version of generate pps
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

num_gen = 80; % Number of generations, can be problem-tuned
num_pop = max(num_pps,1e4); % Population size/gen, can be problem-tuned
do_score = 1; % Flag for turning on scoring of PPs.
score_fcn = @(x,y)ga_score_model(x,y,mu,sigma); % Set the scoring function for the model, this could also be a passed parameter
out_fcn = @(a,b,c)ga_output(a,b,c,num_pps);

if exist('txtout/ga_interim.mat')==2
    delete('txtout/ga_interim.mat')
end

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

%% 3. Get/set integration configuration
mdl_config = set_mdl_config(vdp_model);

%% 4. Prepare scoring function

f = @(p)run_variant_and_score(p,p_names,vdp_model,mdl_config,...
    do_score,state_bnds,score_fcn);

% Set initial population using a latin-hypercube design for the ICs:
% hn = lhsdesign(num_pop,size(p_bnds,1))
% logp_init = log10(bsxfun(@plus,p_bnds(:,1),bsxfun(@times,hn,(p_bnds(:,2)-p_bnds(:,1)))));
pl = log10(p_bnds(:,1));
pu = log10(p_bnds(:,2));
logp_init = zeros(num_p,num_pop);
for i = 1:num_pop
   %logp_init(:,i) = min(pu,max(pl,logp0 + 0.25*randn(num_p,1).*abs(logp0)));
    logp_init(:,i) = log10(p_bnds(:,1) + (p_bnds(:,2)-p_bnds(:,1)).*rand(num_p,1));
end

% Options for optimization - note, problem specific:
opt_ga = gaoptimset('Display', 'off', 'Generations', num_gen, ...
        'PopulationSize', num_pop,...
        'OutputFcn', out_fcn);

[~,~,~,output] = ga(f,num_p,[],[],[],[],pl,pu,[],[],opt_ga);
%disp(output.message)
load('txtout/ga_interim.mat'); % bring back the plausible patients (p_pp)
p_pp = p_pp'; % transpose the plausible patients
pp_yield = size(p_pp,2)./(output.generations*num_pop);

end % function ga_generate_pps
