function plot_gof_vp_vs_pp(t,method_order)
%% Generates Figure 4 in the manuscript
meth_num = numel(method_order);
numpps_lvls = unique(t.num_pps);

p = plot_style;

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
        r3 = r2 & js.selection == 1;
        js2 = js(r2,:); % for each iter
        js3 = js(r3,:);
        yvp_iter(k) = mean(js3.gof); % should all be the same entry again and again
        
        % Need to calculate the GoF for PPs, not pre-computed:
        mu = js2.mu(1,:);
        sigma = reshape(js2.sigma(1,:),numel(mu),numel(mu));
        [~,ypp_iter(k)] = get_prevalence(mu,sigma,log(js2.VPChar),0,js2.betaS,1);
    end
    
    yvp = mean(yvp_iter);
    ypp = mean(ypp_iter);
    
    bar(i-0.2,ypp,0.3,'FaceColor','none',...
        'EdgeColor',p.plot_color(i,:),'LineWidth',p.lw1);
    hold on;
    bar(i+0.2,yvp,0.3,...
        'FaceColor',p.plot_color(i,:),...
        'EdgeColor',p.plot_color(i,:),'LineWidth',p.lw1);
    
    if iter_num > 1
        yvp_err = std(yvp_iter);
        ypp_err = std(ypp_iter);
        plot([i+0.2 i+0.2],[yvp yvp_err+yvp],...
            'Color',p.plot_color(i,:)','LineWidth',p.lw1);
        plot([i+0.1 i+0.3],[yvp_err+yvp yvp_err+yvp],...
            'Color',p.plot_color(i,:)','LineWidth',p.lw1);
        plot([i-0.2 i-0.2],[ypp ypp_err+ypp],...
            'Color',p.plot_color(i,:)','LineWidth',p.lw1);
        plot([i-0.1 i-0.3],[ypp_err+ypp ypp_err+ypp],...
            'Color',p.plot_color(i,:)','LineWidth',p.lw1);
    end
end

h(1) = bar(1000,1,'FaceColor','none','EdgeColor','k','LineWidth',p.lw1);
h(2) = bar(1002,1,'FaceColor','k','EdgeColor','k','LineWidth',p.lw1);

set(gca, 'XTick', 1:meth_num, 'XTickLabel', labels);
set(gca,'FontSize',p.fs1,'LineWidth',p.lw2);
%set(gca,'YScale','log');
ylabel('Goodness-of-Fit')
legend(h,{'Plausible Patients','Virtual Patients'},'Location','NorthEast');
legend boxoff;
xlim([0.5 meth_num+0.5]);

print('figures/rieger-fig05-gof-pp-vs-vp.pdf','-dpdf','-bestfit');

end