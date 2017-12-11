function [mvmodel, r] = correlate_nhanes_chol(do_plots)
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
r = [lnLDL lnHDL lnTC]; % Gather the variables of interest
r = r(~any(isnan(r'))',:); % Clear any rows with NaN, not absolutely necessary, but avoids a warning message
c = corrcoef(r);
% Fit a multivariate Gaussian model to the log-transformed data:
mvmodel = fitgmdist(r,1); % model parameters should match mean(r) and cov(r).

fs1 = 28;

if do_plots
    figure('Name','Correlations','Units','Inches','Position',[1 1 20 10]);
    subplot(3,3,1);
    histogram(lnLDL);
    ylabel('log(LDL)');
    title('log(LDL)');
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    subplot(3,3,2);
    title('log(HDL)');
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    %plot(lnHDL,lnLDL,'o');
    subplot(3,3,3);
    title('log(TC)');
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    %plot(lnTC,lnLDL,'o');
    subplot(3,3,4);
    plot(lnLDL,lnHDL,'o');
    ylabel('log(HDL)');
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    subplot(3,3,5);
    histogram(lnHDL);
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    subplot(3,3,6);
    %plot(lnTC,lnHDL,'o');
    subplot(3,3,7);
    plot(lnLDL,lnTC,'o');
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    ylabel('log(TC)');
    xlabel('log(LDL)');
    subplot(3,3,8);
    plot(lnHDL,lnTC,'o');
    xlabel('log(HDL)');
    set(gca,'Box','on','XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'FontSize',fs1);
    subplot(3,3,9);
    histogram(lnTC);
    xlabel('log(TC)');
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
    
    fprintf('Joint Gaussian Fit\nlnHDL\tlnLDL\tlnTC\nmu:\n');
    disp(mvmodel.mu);
    fprintf('covariance matrix:\n');
    disp(mvmodel.Sigma);
    fprintf('correlation matrix:\n');
    disp(corrcoef(r));
end
% Clear the raw file data:
%clear thdl ttc ttrig

end
