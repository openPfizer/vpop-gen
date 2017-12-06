function plot_histograms(t,method_order)
%% function plot_histograms
% This script loads the data generated for the paper, or that generated from
% RunScript.m and makes figures similar to the main results of the paper.

% For calculating 95% error surfaces from 3D data, this script calls
% error_ellipse.m by AJ Johnson (2004). Available from the MATLAB File
% Exchange http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse
% see error_license.txt for the license for this file.
addpath('./error_ellipse');

npps_lvls = unique(t.num_pps);

r = t.num_pps == max(npps_lvls); % just plotting one iteration here
ts1 = t(r,:); % thin down

for i = 1:numel(method_order)
    ts2 = ts1(ts1.meth_cat == method_order(i),:); % thin for method only
    iter_lvls = unique(ts2.iter);
    ts2 = ts2(ts2.iter == iter_lvls(1),:); % just choose the first iteration
    % Unpack the function arguments:
    pp_obs          = ts2.VPChar;
    select          = logical(ts2.selection);
    p_include       = ts2.ProbInclude;
    betaS           = ts2.betaS;
    mu              = ts2.mu(1,:);
    sigma           = reshape(ts2.sigma(1,:),numel(mu),numel(mu));
    method          = ts2.method{1};
    
    fidelity = 1;
    num_bin = 10;
    num_pps = numel(pp_obs(:,1));
    fig_name = strcat(method,', PP = ',num2str(num_pps),', ');
    
    %% Plotting options:
    fs1 = 24; % FontSize
    fs2 = 14;
    fs3 = 28; % subplot labels
    lw1 = 2; % LineWidth
    xtgap = 0.05;
    
    %% Plot an invisible error ellipse:
    figure(100)
    hellipse = error_ellipse(sigma,mu,'conf',0.95);
    
    xdatas1 = get(hellipse(1), 'XData');
    ydatas1 = get(hellipse(1), 'YData');
    zdatas1 = get(hellipse(1), 'ZData');
    
    xdatas2 = get(hellipse(2), 'XData');
    ydatas2 = get(hellipse(2), 'YData');
    zdatas2 = get(hellipse(2), 'ZData');
    
    xdatas3 = get(hellipse(3), 'XData');
    ydatas3 = get(hellipse(3), 'YData');
    zdatas3 = get(hellipse(3), 'ZData');
    close 100 % clears the invisible error ellipse figure
    
    %% Probability of Inclusion Plot
    a = figure('units','inches','Position',[1 1 20 10],'Name',[fig_name,' Probability of PP Inclusion']);
    
    subplot(1,2,2);
    try
    p_include2 = p_include.*betaS;
    catch
        size(p_include)
        size(betaS)
    end
    bins = -10:0.1:0.5;
    %bins2 = -10:0.1:0.5;
    temp = p_include2>1;
    p_include2(temp) = 1;
    
    histogram(log10(p_include2),bins,'FaceAlpha',0.2,'EdgeAlpha',0.2);
    %bar1 = bar(x,n,'b','LineWidth',1,'BarWidth',1);
    %set(bar1,'FaceColor','none')
    xlim([-10  2]);
    %xlabel('log(p(x=VP|x=PP))');
    xlabel('log(Prob. PP in Vpop)');
    ylabel('Freq.');
    hold on;
    histogram(log10(p_include2(select)),bins,'FaceColor','r','EdgeColor','r');
    %bar2 = bar(x,n,'r','LineWidth',1,'BarWidth',1);
    %set(bar2,'FaceColor','r','EdgeColor','k');
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'M','FontSize',fs3,'VerticalAlignment','Top');
    
    set(gca,'FontSize',fs1,'LineWidth',lw1);
    
    axes('parent',gcf,'position',[0.65 0.6 0.15 0.15]);
    histogram(log10(p_include2(select)),bins,'FaceColor','r','EdgeColor','r');
    xlim([-3  1])
    
    vp_select = pp_obs(select,:);
    
    %xlabel('log(Prob. PP in Vpop)');
    %ylabel('Freq.');
    
    set(gca,'FontSize',fs1,'LineWidth',lw1);
    
    %% Plausible Population Fit Plot
    figure(a);
    subplot(4,6,1)
    plot_pdf_over_bar(log(pp_obs(:,1)),mu(1),sigma(1,1)^0.5,num_bin,'log(HDL_c)')
    fig_a_lim = xlim;
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'A','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    %LDL_c
    subplot(4,6,2);
    plot_pdf_over_bar(log(pp_obs(:,2)),mu(2),sigma(2,2)^0.5,num_bin,'log(LDL_c)')
    fig_b_lim = xlim;
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'B','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    %TC
    subplot(4,6,3);
    plot_pdf_over_bar(log(pp_obs(:,3)),mu(3),sigma(3,3)^0.5,num_bin,'log(TC)')
    fig_c_lim = xlim;
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'C','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    %TC vs HDL_c
    subplot(4,6,7);
    hold on;
    scatter(log(pp_obs(1:fidelity:num_pps,3)),log(pp_obs(1:fidelity:num_pps,1)),2,'r')
    xlabel('log(TC)')
    ylabel('log(HDL_c)');
    
    plot(zdatas3,xdatas3,'k--','LineWidth',2)
    xlim([fig_c_lim])
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'D','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    %TC vs LDL_c
    subplot(4,6,8);
    hold on;
    scatter(log(pp_obs(1:fidelity:num_pps,2)),log(pp_obs(1:fidelity:num_pps,3)),2,'r')
    xlabel('log(LDL_c)')
    ylabel('log(TC)')
    plot(ydatas2,zdatas2,'k--','LineWidth',2)
    xlim([fig_b_lim])
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'E','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    
    %HDL_c vs LDL_c
    subplot(4,6,9)
    hold on;
    scatter(log(pp_obs(1:fidelity:num_pps,1)),log(pp_obs(1:fidelity:num_pps,2)),2,'r')
    xlabel('log(HDL_c)')
    ylabel('log(LDL_c)')
    plot(xdatas1,ydatas1,'k--','LineWidth',2)
    xlim([fig_a_lim]);
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'F','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    %% VPop Fit Plot
    % A repeat of the plausible patient plot
    %figure('units','inches','position',[1 1 20 10],'Name',[fig_name,' Virtual Population Fit']);
    
    %HDL_c
    subplot(4,6,13)
    plot_pdf_over_bar(log(vp_select(:,1)),mu(1),sigma(1,1)^0.5,num_bin,'log(HDL_c)')
    fig_a_lim = xlim;
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'G','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    %LDL_c
    subplot(4,6,14)
    plot_pdf_over_bar(log(vp_select(:,2)),mu(2),sigma(2,2)^0.5,num_bin,'log(LDL_c)')
    fig_b_lim = xlim;
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'H','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    %TC
    subplot(4,6,15);
    plot_pdf_over_bar(log(vp_select(:,3)),mu(3),sigma(3,3)^0.5,num_bin,'log(TC)')
    fig_c_lim = xlim;
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'I','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    %
    % xdatas = get(hellipse(3), 'XData');
    % ydatas = get(hellipse(3), 'YData');
    % zdatas = get(hellipse(3), 'ZData');
    
    
    %TC vs HDL_c
    subplot(4,6,19)
    hold on;
    
    scatter(log(vp_select(:,3)),log(vp_select(:,1)),6,'r');
    plot(zdatas3,xdatas3,'k--','LineWidth',2);
    xlim([fig_c_lim]);
    xlabel('log(TC)');
    ylabel('log(HDL_c)');
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'J','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    %
    % xdatas = get(hellipse(2), 'XData');
    % ydatas = get(hellipse(2), 'YData');
    % zdatas = get(hellipse(2), 'ZData');
    
    %TC vs LDL_c
    subplot(4,6,20)
    hold on;
    scatter(log(vp_select(:,2)),log(vp_select(:,3)),6,'r')
    plot(ydatas2,zdatas2,'k--','LineWidth',2)
    xlim([fig_b_lim])
    xlabel('log(LDL_c)')
    ylabel('log(TC)')
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'K','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    %
    % xdatas = get(hellipse(1), 'XData');
    % ydatas = get(hellipse(1), 'YData');
    % zdatas = get(hellipse(1), 'ZData');
    
    %HDL_c vs LDL_c
    subplot(4,6,21)
    hold on;
    scatter(log(vp_select(:,1)),log(vp_select(:,2)),6,'r')
    plot(xdatas1,ydatas1,'k--','LineWidth',2)
    xlim([fig_a_lim])
    xlabel('log(HDL_c)');
    ylabel('log(LDL_c)');
    
    h = get(gca);
    text(h.XLim(1)+xtgap*(h.XLim(2)-h.XLim(1)),h.YLim(2),'L','FontSize',fs3,'VerticalAlignment','Top');
    ylim(h.YAxis.Limits);
    set(gca,'FontSize',fs2,'LineWidth',lw1);
    
    figure(a);
    hold on;
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    h = text(ax1,0.067,0.74,'PPs','FontSize',fs3,'VerticalAlignment','Middle','HorizontalAlignment','Center');
    h.Rotation = 90;
    h = text(ax1,0.067,0.30,'VPs','FontSize',fs3,'VerticalAlignment','Middle','HorizontalAlignment','Center');
    h.Rotation = 90;
    
end

end % plot_histograms

% *************************************************************************

function plot_pdf_over_bar(data,mu,sigma,num_bins,x_label)

hold on;

[n,x] = hist(data,num_bins);

area2 = sum(n.*(x(2)-x(1)));

stairs((x-(x(2)-x(1))/2),n/area2,'r','LineWidth',3)
xlim([min(x) max(x)])

%% Plot the PDF
xp = min(x):0.01:max(x);
y = normpdf(xp,mu,sigma);
plot(xp,y,'--k','LineWidth',5)

xlabel(x_label);
ylabel('PDF');
ylim([0,1.2*max(y)])

end %


