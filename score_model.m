function score = score_model(mdl_ss,state_range)
%   function score_model
%       Generic scoring function for model assuming a PP/VP's score is 0 if
%       all model states are within observation limits and > 0 otherwise.
%       This function does not assess the likelihood of a PP/VP, only if
%       they are plausible in that they are "in-range" of the observables.
%
%   Inputs:
%       mdl_ss - steady state states from the model of interest.
%       state_range - steady state obserable ranges (col. 1 - lower bnds,
%           col. 2 - upper bnds).
%
%   Outputs:
%       score - scalar, 0 = all model values within observable ranges
%       (plausible patient). > 0 otherwise.
%
%   Note: this is a slightly modified version of the original
%   "ScoreModelSS" function.
%

%% Filter for the date we are interested in for fitting:
state_midpt = (state_range(:,1)+state_range(:,2))/2;
%mdl_ss = mdl_ss(state_flg);
%state_range = state_range(state_flg);
%state_midpt = state_midpt(state_flg);

%% Calculate score:
score = (state_midpt - mdl_ss(:)).^2 - (state_midpt - state_range(:,1)).^2;
score = score./state_midpt.^2;
score = sum(score(score>=0));

end % function score_model
