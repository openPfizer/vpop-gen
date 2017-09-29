function [ck_p,p,dsq] = ellipse_find_cks(num_r, do_plot, thresh)
%% function ellipse_find_cks(num_r, do_plot)
%   Calculate the Ck distance for ellipse contours using equal percentiles
%   of the observables from NHANES.
%
%   Inputs:
%       num_r - number of regions to seek, the CDF will be equally divided
%       into this number of regions (e.g., num_r = 4, p = 25, 50, 75, 97.5).
%       do_plot - boolean for optional plotting of CDF vs Ck values.
%       Default is 0.
%   Outputs:
%       ck_p - Cks for percentile p
%       p - Percentiles for the Cks (e.g., 25, 50, 75, 97.5...)
%       d - Distances for all NHANES observations
%
%   Note: There is a hard-coded upper cutoff of 2.5%, this can be changed
%   in the line where p is calculated, suggest not using 0 though due to
%   the logarithmic nature of the data.
%

%% Check inputs:
if nargin == 1
    do_plot = 0;
    thresh = 0;
elseif nargin == 2
    thresh = 0;
end

%% Retrieve NHANES data:
[s,obs] = correlate_nhanes_chol(0);
nobs = size(obs,1);

p = 100 - linspace(100,thresh,num_r+1)';
p(1) = []; % we don't need 0th percentile

%% Calculate the distanace for all observations:
dsq = zeros(nobs,1);
for i = 1:nobs
    dsq(i) = (obs(i,:) - s.mu)/s.Sigma*(obs(i,:)-s.mu)'; % squared distance
end

%% Find Cks from percentiles p:
ck_p = sqrt(prctile(dsq,p)'); % taking the sqrt

%% Optional plotting:
if do_plot
    figure;
    %plot the CDF of d...
    [f,x] = ecdf(dsq);
    plot(x,f);
    hold on;
    plot(ck_p.^2,p/100,'or','MarkerSize',12,'MarkerFaceColor','r');
    %xlim([0 20]);
end

end % function ellipse_find_cks
