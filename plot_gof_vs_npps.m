function plot_gof_vs_npps(t,method_order)
%% Generates Figure 2 in the manuscript
line_style = {'-';'--';':';'-.'};
marker_style = {'o';'s';'d';'x'};

plot_color = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

a = figure('units','inches','position',[1 1 20 10],'Name','GoF vs. NPPs');

numpps_lvls = unique(t.num_pps);
numpps_num = numel(numpps_lvls);

for i = 1:numel(method_order)
    ym = [];
    ymerr = [];
    kcount = 0;
    for k = 1:numpps_num
        r = t.meth_cat == method_order(i) & t.selection & t.num_pps == numpps_lvls(k);
        if sum(r) >= 1
            kcount = kcount + 1;
            x(kcount) = numpps_lvls(k);
            ym(kcount) = mean(t.gof(r));
            if sum(r) > 1
                ymerr(kcount) = std(t.gof(r));
            else
                ymerr(kcount) = 0;
            end
        end
    end
    
    if kcount > 0
        figure(a);
        ls = strcat(line_style{i},marker_style{i});
        h(i) = plot(x,ym,ls,'LineWidth',3,'Color',plot_color(i,:),...
            'MarkerFaceColor',plot_color(i,:),'MarkerSize',20);
        hold on;
        for k = 1:numel(x)
            plot([x(k) x(k)],[ym(k)-ymerr(k) ym(k)+ymerr(k)],line_style{i},'LineWidth',3,'Color',plot_color(i,:));
        end
    end
end % Method Loop

%% Finish the plot:
figure(a);
xlabel('Number of Plausible Patients');
ylabel('Goodness-of-Fit');
legend(h,method_order,'Location','NorthWest');
legend boxoff
set(gca,'FontSize',28,'LineWidth',2);
set(gca,'XScale','log');

end % EoF
