function plot_time_summary_bar(t,method_order)
%% Generates Figure 3 in the manuscript
fs1 = 18;
lw1 = 1;
lw2 = 2;

plot_color = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

ymean_pp = zeros(numel(method_order),1);
ymean_vp = zeros(numel(method_order),1);
yerr_pp = zeros(numel(method_order),1);
yerr_vp = zeros(numel(method_order),1);

numpps_lvls = unique(t.num_pps);

figure('units','Inches','position',[1 1 20 10],'Name','Computational Time Summary');
for i = 1:numel(method_order)
    r = t.meth_cat == method_order(i) & t.num_pps == max(numpps_lvls) & t.selection;
    js = t(r,:);
    
    iter_lvls = unique(js.iter);
    iter_num = numel(iter_lvls);
    
    time_per_pp = zeros(iter_num,1);
    time_per_vp = zeros(iter_num,1);
    for k = 1:iter_num
        r2 = js.iter == iter_lvls(k);  
        time_per_pp(k) = js(r2,:).time_per_pp(1);
        time_per_vp(k) = js(r2,:).time_per_vp(1);
    end
    
    ymean_pp(i) = mean(time_per_pp);
    ymean_vp(i) = mean(time_per_vp);

    bar(i-0.2,ymean_pp(i),0.3,'FaceColor','none','EdgeColor',plot_color(i,:),'LineWidth',lw2);
    hold on;
    bar(i+0.2,ymean_vp(i),0.3,'FaceColor',plot_color(i,:),'EdgeColor',plot_color(i,:),'LineWidth',lw2);
    
    if numel(js) > 1
        % Adding error bars to bar plots in Matlab remains incredibly
        % frustrating...
        yerr_pp(i)  = std(time_per_pp);
        yerr_vp(i) = std(time_per_vp);
        plot([i-0.2 i-0.2],[ymean_pp(i) yerr_pp(i)+ymean_pp(i)],'Color',plot_color(i,:)','LineWidth',lw2);
        plot([i+0.2 i+0.2],[ymean_vp(i) yerr_pp(i)+ymean_vp(i)],'Color',plot_color(i,:)','LineWidth',lw2);
        plot([i-0.3 i-0.1],[yerr_pp(i)+ymean_pp(i) yerr_pp(i)+ymean_pp(i)],'Color',plot_color(i,:)','LineWidth',lw2);
        plot([i+0.1 i+0.3],[yerr_pp(i)+ymean_vp(i) yerr_pp(i)+ymean_vp(i)],'Color',plot_color(i,:)','LineWidth',lw2);
    end
    labels{i} = method_order{i};
end

h(1) = bar(1000,1,'FaceColor','none','EdgeColor','k','LineWidth',lw2);
h(2) = bar(1002,1,'FaceColor','k','EdgeColor','k','LineWidth',lw2);

set(gca, 'XTick', 1:numel(method_order), 'XTickLabel', labels);
set(gca,'FontSize',fs1,'LineWidth',lw1);
set(gca, 'YScale','log');
%set(gca,'YScale','log');
ylabel('Time (seconds)')
legend(h,{'Time Per Plausible Patient','Time Per Virtual Patient'});
legend boxoff;
xlim([0 numel(method_order)+1]);

print('figures/rieger-fig03-time-bar.pdf','-dpdf','-bestfit');

end % function plot_time_summary_bar
