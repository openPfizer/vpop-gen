function plot_gof_vp_vs_pp(t,method_order)
%% Generates Figure 4 in the manuscript
meth_num = numel(method_order);
numpps_lvls = unique(t.num_pps);

plot_color = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

m = figure('units','inches','position',[1 1 20 10],'name','GoF VP vs. PP');
for i = 1:meth_num
    labels{i} = method_order{i};
    r = t.meth_cat == method_order(i) & t.num_pps == max(numpps_lvls);
    js = t(r,:);
    iter_lvl = unique(js.iter);
    iter_num = numel(iter_lvl);
    
    yvp_iter = zeros(iter_num,1);
    ypp_iter = zeros(iter_num,1);
    
    for k = 1:iter_num
        r2 = js.iter == iter_lvl(k);
        js2 = js(r2,:); % for each iter
        yvp_iter(k) = mean(js2.gof); % should all be the same entry again and again
        
        % Need to calculate the GoF for PPs, not pre-computed:
        mu = js2.mu(1,:);
        sigma = reshape(js2.sigma(1,:),numel(mu),numel(mu));
        [~,ypp_iter(k)] = get_prevalence(mu,sigma,log(js2.VPChar),0,js2.betaS,1);
    end
    
    yvp = mean(yvp_iter);
    ypp = mean(ypp_iter);
    
    bar(i-0.2,ypp,0.3,'FaceColor','none','EdgeColor',plot_color(i,:),'LineWidth',3);
    hold on;
    bar(i+0.2,yvp,0.3,'FaceColor',plot_color(i,:),'EdgeColor',plot_color(i,:),'LineWidth',3);
    
    if iter_num > 1
        yvp_err = std(yvp_iter);
        ypp_err = std(ypp_iter);
        plot([i+0.2 i+0.2],[yvp yvp_err+yvp],'Color',plot_color(i,:)','LineWidth',3);
        plot([i+0.1 i+0.3],[yvp_err+yvp yvp_err+yvp],'Color',plot_color(i,:)','LineWidth',3);
        plot([i-0.2 i-0.2],[ypp ypp_err+ypp],'Color',plot_color(i,:)','LineWidth',3);
        plot([i-0.1 i-0.3],[ypp_err+ypp ypp_err+ypp],'Color',plot_color(i,:)','LineWidth',3);
    end
end

h(1) = bar(1000,1,'FaceColor','none','EdgeColor','k','LineWidth',3);
h(2) = bar(1002,1,'FaceColor','k','EdgeColor','k','LineWidth',3);

set(gca, 'XTick', 1:meth_num, 'XTickLabel', labels);
set(gca,'FontSize',28,'LineWidth',2);
%set(gca,'YScale','log');
ylabel('Goodness-of-Fit')
legend(h,{'Plausible Patients','Virtual Patients'},'Location','NorthEast');
legend boxoff;
xlim([0 meth_num+1]);

end