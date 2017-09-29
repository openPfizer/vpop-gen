function plot_violin(t,method_order)
%% function plot_violin
%   Creates violin plots as used in the Supplementary Materials
addpath('./distributionPlot');

fs1 = 28;
fs2 = 16;

npps_lvls = unique(t.num_pps);
r = t.num_pps == max(npps_lvls); % just plotting one iteration here
ts1 = t(r,:); % thin down to just maximum n_pps

for i = 1:numel(method_order)
    ts2 = ts1(ts1.meth_cat == method_order(i),:); % thin for method only
    iter_lvls = unique(ts2.iter);
    ts2 = ts2(ts2.iter == iter_lvls(1),:); % just choose the first iteration
    
    %% Unpack sm structure:
    %sm = t.sm;
    p_final     = ts2.Pfinal';
    select      = logical(ts2.selection);
    p_bnds      = [ts2.p_lower(1,:)' ts2.p_upper(1,:)'];
    method      = ts2.method{1};
    p_scaled_select = ts2.PScaled(select,:);
    p_scaled = ts2.PScaled;
    
    num_pps = numel(p_final(1,:));
    plot_head = strcat(method,' PP = ',num2str(num_pps));
    
    a = figure('units','inches','position',[1 1 20 10],'Name',[plot_head ' Violin Plot']);
    
    %% Correlation Heatmaps:
    % p_upper_mat = repmat(p_bnds(:,2),1,num_pps);
    % p_lower_mat = repmat(p_bnds(:,1),1,num_pps);
    % p_scaled = (p_final-p_lower_mat)./(p_upper_mat-p_lower_mat);
    % p_scaled = p_scaled';
    % p_scaled_select = p_scaled(select,:);
    % % Prepare an output:
    % sm_out = sm;
    % sm_out.P_scaled = p_scaled_select; % Note: selected version is passed back
    
    C1 = corrcoef(p_scaled);
    C2 = corrcoef(p_scaled_select);
    Cmax = max(max(max(C1-eye(size(C1)))),max(max(C2-eye(size(C2)))));
    Cmin = min(min(min(C1-eye(size(C1)))),min(min(C2-eye(size(C2)))));
    
    %% Plausible Patients Heatmap:
    subplot(2,2,2);
    imagesc(C1,[Cmin Cmax])
    set(gca,'xtick',[],'ytick',[])
    colorbar
    axis square;
    xlabel('Scaled Parameter');
    ylabel('Scaled Parameter');
    title('Parameter Correlation');
    set(gca,'FontSize',fs2,'LineWidth',2);
    
    %% Virtual Patients Heatmap:
    subplot(2,2,4)
    imagesc(C2,[Cmin Cmax])
    set(gca,'xtick',[],'ytick',[])
    colorbar
    axis square
    xlabel('Scaled Parameter');
    ylabel('Scaled Parameter');
    set(gca,'FontSize',fs2,'LineWidth',2);
    
    %% Plausible Patient Violin:
    subplot(2,2,1);
    
    for i = 1:size(p_bnds,1)
        ind = num2str(i);
        p_name_tmp = strcat('k_{',ind,'}');
        p_names_plot{i} = p_name_tmp;
    end
    
    distributionPlot(p_scaled,'showMM',6,'xNames',p_names_plot,'histOpt',1)
    ylim([0 1]);
    alpha(0.3);
    ylabel('Scaled Parameter');
    title('Parameter Distribution');
    set(gca,'FontSize',fs2,'LineWidth',2);
    set(gca,'Position',[0.13 0.584 0.425 0.341]);
    
    %% Virtual Patient Violin Plot:
    subplot(2,2,3);
    
    distributionPlot(p_scaled_select,'showMM',6,'xNames',p_names_plot,'histOpt',1)
    ylim([0 1]);
    alpha(0.3);
    
    ylabel('Scaled Parameter');
    set(gca,'FontSize',fs2,'LineWidth',2);
    hold on;
    set(gca,'Position',[0.13 0.11 0.425 0.341]);
    
    %% Add Row Labels and Subpanel Labels
    figure(a);
    hold on;
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    h = text(ax1,0.033,0.74,'PPs','FontSize',36,'VerticalAlignment','Middle','HorizontalAlignment','Center');
    h.Rotation = 90;
    h = text(ax1,0.033,0.30,'VPs','FontSize',36,'VerticalAlignment','Middle','HorizontalAlignment','Center');
    h.Rotation = 90;
    text(ax1,0.085,0.925,'A','FontSize',fs1);
    text(ax1,0.6,0.925,'B','FontSize',fs1);
    text(ax1,0.085,0.425,'C','FontSize',fs1);
    text(ax1,0.6,0.425,'D','FontSize',fs1);
    
end % function plot_violin

end