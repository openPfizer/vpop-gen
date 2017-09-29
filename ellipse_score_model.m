function [score,score_ellip,score_sa] = ellipse_score_model(mdl,state_bnds, mu, s, ck)
%% function ellipse_score_model
%   Scores plausible patients from the ellipse method. Built on top of the
%   simulated annealing method.
%
%   Inputs:
%       mdl - Vector of model states.
%       state_bnds - Matrix (nx2) of lower and upper bounds on model states
%       mu  - maximum likelihood for observables from NHANES data
%       s   - covariance matrix for NHANES observables
%       ck  - ellipse boundary
%
%   Outputs:
%       score       - overall score (simulated annealing + ellipse)
%       score_ellip - contribution of ellipse method to scoring
%
%   Notes:
%       1)  lambda is a weighting parameter for reweighing the simulated
%           annealing vs. the ellipse score. Presently, this parameter is
%           hard-coded here, but it could easily be passed as a tunable
%           parameter.
%       2)  PP is only considered invalid if any of the states is < 0.
%           Otherwise they are scored as is.
%

%% Score calculation
lambda = 0.1; % hard-coded weight on ellipse method score.
mdls = mdl(:); % make sure the inputs are a column-vector

if ~any(mdls<0)
    
    % Simulated annealing score, uses same function as simulated annealing
    % method:
    score_sa = score_model(mdls,state_bnds);
    
    % Create the log(observables) with a unit conversion:
    lnhdl = log((mdls(5) + mdls(6))*(38.66/2.79));
    lntc =  log((mdls(7) + mdls(5) + mdls(6))*(38.66/2.79));
    lnldl = log(mdls(9)*(38.66/2.79));
    lnobs = [lnhdl lnldl lntc];
    
    % Ellipse-specific score:
    score_ellip = max(0,(lnobs - mu)/s*(lnobs(:) - mu(:)) - ck^2);
    
    % Total score, weighting by lambda:
    score = score_sa + lambda*score_ellip;
else
    % A model value is out-of-bounds, flag as bad and return:
    score_sa = 1e256;
    score_ellip = 1e256;
    score = 1e256;
end

end % function ellipse_score_model
