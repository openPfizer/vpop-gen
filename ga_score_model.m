function score = ga_score_model(mdl_ss,state_range,mu,sigma)
%   function ga_score_model
%       Builds on generic score_model function by adding in some knowledge
%       of the multivariate PDF
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
%   "ScoreModelSSGA" function.
%

%% Filter for the date we are interested in for fitting:
state_midpt = (state_range(:,1)'+state_range(:,2)')/2;

%% Calculate score:
score = (state_midpt - mdl_ss).^2 - (state_midpt - state_range(:,1)').^2;
score = score./state_midpt.^2;
score = sum(score(score>=0));

% Create the NHANES observables:
mdl_ss = mdl_ss*(38.66/2.79); %convert units to match NHANES
lnhDL = log((mdl_ss(5) + mdl_ss(6)));
lntc =  log((mdl_ss(7) + mdl_ss(5) + mdl_ss(6)));
lnldl = log((mdl_ss(9)));

score = score - mvnpdf([lnhDL lnldl lntc],mu,sigma);

end % function score_model
