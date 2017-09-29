function [score_final,x_mdl] = run_variant_and_score(logp,p_names,...
    vdp_model,mdl_config,do_score,state_range,score_f)
%% function run_variant_and_score
%   Takes in a plausible patient values and simulates to steady state,
%   variant is then score using a passed scoring function.
%
%   Inputs:
%       logp - log10 transformed parameter vector for PP/VP
%       p_names - names of parameters in logp vector
%       vdp_model - SimBiology exported model
%       mdl_config -
%       state_range - Valid ranges for state variables (from observables)
%       state_flag - Logical selector for states to be scored against model
%       score_f - Scoring function to turn model values and observable
%           ranges into a fitness value. Required inputs are the model
%           steady state values, observable ranges, and state_flag logical.
%   Outputs:
%       score_final - score of the plausible patient assessed in score_f
%       x_mdl - steady state (final value) point for model
%

%% 1. Setup
score_final = 1e256;
x_mdl = NaN(9,1);
fail_flag = 0;
p = 10.^logp;       % Parameters are passed as log10(p)

%% 2. Create New Variant
% This step adds a "variant" to the simbiology model, that is an
% alternative parameterization of the model, which is simulated below.
try
    new_var = addvariant(vdp_model,'NewP2');
catch
    delete(vdp_model.variant(numel(vdp_model.variant)));
    new_var = addvariant(vdp_model,'NewP2');
end

% Writes the proposed parameters to the new variant:
for i=1:numel(p)
    addcontent(new_var,{'parameter',p_names{i},'Value',p(i)})
end

%% 2. Simulate and Score Results
% Use try/catch - may fail due to integration issues.
% If integration fails, assume something bad has happened and we don't
% want that PP (add a very high score)

try
    [~,x1] = sbiosimulate(vdp_model,mdl_config,new_var); % Simulated the new variant, usually to steady state
catch
    fail_flag = 1;
    score_final = 1e256;
    x_mdl = NaN(9,1);
end

if ~fail_flag && ~any(imag(x1(end,:))) && ~any(any(x1<=0))
    x_mdl = x1(end,:)';
    if do_score
        score_final = feval(score_f,x1(end,:),state_range);
    end
end

%% 3. Clean-up (delete variants)
% Delete the last variant that was added to the model - need to do this
% to keep things organized
delete(vdp_model.variant(numel(vdp_model.variant)));

end

