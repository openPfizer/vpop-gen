function [mdl, r] = correlate_nhanes_chol(do_plots)
%clear;clc;close all;rng('shuffle');
if nargin == 0
    do_plots = 0;
end

% Load the NHANES data:
thdl = xptread('HDL_G.XPT');
ttc = xptread('TCHOL_G.XPT');
ttrig = xptread('TRIGLY_G.XPT');

% Merge based on the key SEQN:
c1 = join(thdl,ttc,'Keys','SEQN');
c1 = join(ttrig,c1,'Keys','SEQN');

% Log-transform the relevant columns:
lnHDL = log(c1.LBDHDD);
lnLDL = log(c1.LBDLDL);
lnTC = log(c1.LBXTC);

%% Find the joint probability distribtion:
r = [lnHDL lnLDL lnTC]; % Gather the variables of interest
r = r(~any(isnan(r'))',:); % Clear any rows with NaN, not absolutely necessary, but avoids a warning message
c = corrcoef(r);
% Fit a multivariate Gaussian model to the log-transformed data:
mdl = fitgmdist(r,1); % model parameters should match mean(r) and cov(r).
labels = {'log(HDL)';'log(LDL)';'log(TC)'};

fs1 = 28;
if do_plots
    figure('Name','Correlations','Units','Inches','Position',[1 1 20 10]);
    
    subplot(3,3,1);
    histogram(r(:,1));
    ylabel(labels{1});
    title(labels{1});
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    
    subplot(3,3,2);
    title(labels{2});
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    
    subplot(3,3,3);
    title(labels{3});
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    
    % Row 2:
    subplot(3,3,4);
    plot(r(:,2),r(:,1),'o');
    ylabel(labels{2});
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    
    subplot(3,3,5);
    histogram(r(:,2));
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    
    subplot(3,3,6);
    
    % Row 3:
    subplot(3,3,7);
    plot(r(:,3),r(:,1),'o');
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    ylabel(labels{3});
    xlabel(labels{1});
    
    subplot(3,3,8);
    plot(r(:,2),r(:,3),'o');
    xlabel(labels{2});
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    
    subplot(3,3,9);
    histogram(r(:,3));
    xlabel(labels{3});
    
    % Upper-diagonal correlations:
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    subplot(3,3,2);
    text(0.5,0.5,...
        strjoin({'\rho = ' num2str(round(c(1,2),1))}),...
        'HorizontalAlignment','center','FontSize',28);
    set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'Box','on');
    
    subplot(3,3,3);
    text(0.5,0.5,...
        strjoin({'\rho = ' num2str(round(c(1,3),1))}),...
        'HorizontalAlignment','center','FontSize',28);
    set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'Box','on');
    
    subplot(3,3,6);
    text(0.5,0.5,...
        strjoin({'\rho = ' num2str(round(c(2,3),1))}),...
        'HorizontalAlignment','center','FontSize',28);
    set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'Box','on');
    
    % Diagnostics on fitting:
    fprintf('Joint Gaussian Fit\nlnHDL\tlnLDL\tlnTC\nmu:\n');
    disp(mdl.mu);
    fprintf('covariance matrix:\n');
    disp(mdl.Sigma);
    fprintf('correlation matrix:\n');
    disp(corrcoef(r));
end

end
