function mdl_config = set_mdl_config(mdl)
%% function set_mdl_config
%
%

%% Fetch current configset from mdl and set solver options:
mdl_config = getconfigset(mdl);
set(mdl_config, 'SolverType', 'ode15s');
set(mdl_config.SolverOptions, 'RelativeTolerance', 1.0e-4);
set(mdl_config.SolverOptions, 'AbsoluteTolerance', 1.0e-6);

end
