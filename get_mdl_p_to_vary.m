function [p_name, p_values] = get_mdl_p_to_vary(mdl)
%% function get_mdl_p_to_vary
%   Filters the SimBiology model for parameters that are actual constants
%   of the model and not rules or zero @ baseline.
%   
%       Inputs:
%           mdl - SimBiology object exported from SimBiology and loaded
%               previously.
%       Outputs:
%           p_name - Vector of parameter names in the model that can be
%               varied by the optimization.
%           p_values - Vector of parameter values (bsln.) in the model that
%               can be varied by the optimization.
%
%   This function is a cleaned-up version of the original ParametersToVary function.
%

%% Get a List of the repeated assignments to exclude.
rules = mdl.Rules;
all_p = sbioselect(mdl.parameters); % List of all parameters in the model

% Loop through all parameters in the model. Compare the name of each
% parameter against the names of the Rules List fetched from the model.
% Those parameters should be ignored in the fitting, the other parameter
% values ("constants") are fair game for fitting:
k = 1; % constant parameter counter
for i=1:numel(all_p)
    check(i) = 0; % reset the "found" flag
    tmp_name = all_p(i).Name; % save the current parameter name
    
    % Check if rule exists, or is on the no fly list
    for j = 1:numel(rules)
        name_chk = strcmp(rules(j).parserule, tmp_name);
        if name_chk == 1
            check(i) = check(i) + 1;
        end
    end
    
    if check(i) == 0
        % We did not find the parameter amongst the Rules, save the
        % parameter name and value:
        p_name{k} = tmp_name;
        p_values(k) = all_p(i).Value;
        k = k + 1;
    end
end

%% Tidy up outputs
p_name = p_name';
p_values = p_values';
% Remove parameters that are 0 @ baseline:
zero_p = p_values==0; % Find parameter values that are zero
p_name = p_name(~zero_p); % Filter out zero parameter names
p_values = p_values(~zero_p); % Filter out zero parameters

end % function get_mdl_p_to_vary
