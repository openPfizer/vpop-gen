function plot_orthog(t,method_order)
%% function plot_orthog
% Function to check the pairwise projections of vectors.
% Worst results are skewed to the right, better to the left.
% Used to think about how independent large parameter sets are. In
% particular, for virtual populations we would prefer orthogonal parameter
% sets if possible for hypothesis exploration

% Too expensive to do every pair-wise comparison.
% assume that if this number is large enough we get a good sample of the true distribution
num_case = 5e4;
ecdf_save = linspace(0.001,0.999,50);

p = plot_style;

std_thin = 5; % factor to thin out the points plotted

num_pps_lvls = unique(t.num_pps);
r = t.selection & t.num_pps == max(num_pps_lvls);
js = t(r,:);

%% Create the ecdf for the random vector:

for i_uni = 1:10
    p_rand = rand(size(js.PScaled,2),1000) - 0.5; % zero-centered
    [~,~,rand_x(:,i_uni)] = get_ecdf_dot(p_rand,num_case,ecdf_save);
end
h = figure('units','inches','position',[1 1 20 10],'Name','Orthogonality of Parameters, CDF');
g(1) = plot(mean(rand_x,2),ecdf_save,'-k','LineWidth',p.lw1);
hold on;

meth_x_mean = zeros(numel(ecdf_save),numel(method_order));
meth_x_std = zeros(numel(ecdf_save),numel(method_order));

for i_met = 1:numel(method_order)
    r2 = js.meth_cat == method_order(i_met);
    js2 = js(r2,:);
    iter_lvls = unique(js2.iter);
    num_iter = numel(iter_lvls);    
    meth_x = zeros(numel(ecdf_save),num_iter);
    
    for i_iter = 1:num_iter
        r3 = js2.iter == iter_lvls(i_iter);
        p_meth = (js2{r3,'PScaled'} - 0.5)'; % zero-centered
        %size(p_meth)
        [~,~,meth_x(:,i_iter)] = get_ecdf_dot(p_meth,num_case,ecdf_save);
    end
    
    meth_x_mean(:,i_met) = mean(meth_x,2);
    meth_x_std(:,i_met) = std(meth_x,[],2);
    ls = strcat(p.line_style{i_met},p.marker_style{i_met});
    
    plot(meth_x_mean(:,i_met),...
        ecdf_save,p.line_style{i_met},...
        'Color',p.plot_color(i_met,:),'LineWidth',p.lw1);
    
    g(i_met+1) = plot(meth_x_mean(std_thin:std_thin:numel(ecdf_save),i_met),...
        ecdf_save(std_thin:std_thin:numel(ecdf_save)),p.marker_style{i_met},...
        'Color',p.plot_color(i_met,:),...
        'LineWidth',p.lw1,'MarkerSize',p.marker_size);
    
    for k = std_thin:std_thin:numel(ecdf_save)
        ytmp = [ecdf_save(k) ecdf_save(k)];
        xtmp = [meth_x_mean(k,i_met)-meth_x_std(k,i_met) meth_x_mean(k,i_met)+meth_x_std(k,i_met)];
        ytmp_hat = [ytmp(1)-0.01 ytmp(2)+0.01];
        plot(xtmp,ytmp,'-','Color',p.plot_color(i_met,:),'LineWidth',p.lw1);
        plot([xtmp(1) xtmp(1)],ytmp_hat,'-',...
            'Color',p.plot_color(i_met,:),'LineWidth',p.lw1);
        plot([xtmp(2) xtmp(2)],ytmp_hat,'-',...
            'Color',p.plot_color(i_met,:),'LineWidth',p.lw1);
    end  
end % for i_met

figure(h);
xlabel('Dot Product of Scaled VP Parameters');
ylabel('Cumulative Distribution');
ylim([0 1.05]);
set(gca,'Ytick',0:0.2:1);
set(gca,'FontSize',p.fs1,'LineWidth',p.lw2);
legend(g,['Uniform';method_order],'Location','NorthWest');
legend boxoff;

print('figures/rieger-fig04-orthogonality-of-vps.pdf','-dpdf','-bestfit');

end % function plot_orthog
%% ***********************************************************************

function [f,x,xout] = get_ecdf_dot(p,n,fin)

y = zeros(n,1);
for i = 1:n
    % Choose two indices:
    idx(1) = 1; idx(2) = 1;
    while idx(1) == idx(2)
        idx = randi([1 size(p,2)],2,1);
    end
    y(i) = dot(p(:,idx(1)),p(:,idx(2)))/(norm(p(:,idx(1))*norm(p(:,idx(2)))));
end

[f,x] = ecdf(y);

xout = interp1(f,x,fin); % interpolate to exact points

end