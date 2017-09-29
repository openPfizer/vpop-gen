function score = mh_score_model(mdl,state_bnds,mu,sigma)

mdls = mdl(:);
if all(mdls(:)>=state_bnds(:,1) & mdls<=state_bnds(:,2))
    % Create the model observables:
    mdls = mdls*(38.66/2.79); %convert units to match NHANES
    lnhDL = log((mdls(5) + mdls(6)));
    lntc =  log((mdls(7) + mdls(5) + mdls(6)));
    lnldl = log((mdls(9)));
    % In the M-H algorithm, higher is better for score, unlike other
    % algorithms:
    score = mvnpdf([lnhDL lnldl lntc],mu,sigma);
else
    %[state_bnds(:,1) mdls state_bnds(:,2)]
    score = eps; % Flag as a really, really small value
    %error('full stop.');
end

end % function mh_score_model
