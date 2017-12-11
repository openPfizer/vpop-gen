function plot_gof_vs_npps(t,method_order)
%% Generates Figure 2 in the manuscript
p = plot_style;

a = figure('units','inches','position',[1 1 20 10],'Name','GoF vs. NPPs');

numpps_lvls = unique(t.num_pps);
numpps_num = numel(numpps_lvls);

for i = 1:numel(method_order)
    ym = [];
    ymerr = [];
    x = [];
    kcount = 0;
    for k = 1:numpps_num
        r = t.meth_cat == method_order(i) & t.selection & t.num_pps == numpps_lvls(k);
        atmp = [];
        if sum(r) > 0
            kcount = kcount + 1;
            s = t(r,:);
            iterlvls = unique(s.iter);
            for j = 1:numel(iterlvls)
                r2 = s.iter == iterlvls(j);
                atmp(j) = mean(s.gof(r2));
            end
            x(kcount) = numpps_lvls(k);
            ym(kcount) = mean(atmp);
            if numel(iterlvls) > 1
               ymerr(kcount) = std(atmp); 
            else
               ymerr(kcount) = NaN;
            end
        end
    end
    
    if kcount > 0
        figure(a);
        ls = strcat(p.line_style{i},p.marker_style{i});
        h(i) = plot(x,ym,ls,'LineWidth',p.lw1,'Color',p.plot_color(i,:),...
            'MarkerFaceColor','none','MarkerSize',p.marker_size);
        hold on;
        for k = 1:numel(x)
            plot([x(k) x(k)],[ym(k)-ymerr(k) ym(k)+ymerr(k)],...
                p.line_style{i},'LineWidth',p.lw1,'Color',p.plot_color(i,:));
        end
    end
end % Method Loop

%% Finish the plot:
figure(a);
xlabel('Number of Plausible Patients');
ylabel('Goodness-of-Fit');

if kcount > 1
    legend(h,method_order,'Location','NorthEast');
    legend boxoff
end

set(gca,'FontSize',p.fs1,'LineWidth',p.lw2);
set(gca,'XScale','log');

print('figures/rieger-fig02-gof-vs-pps.pdf','-dpdf','-bestfit');

end % EoF
