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

plot_color = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840]; % hard-coding in the plot color order

line_style = {'-';'--';':';'-.'};
marker_style = {'o';'s';'d';'x'};

std_thin = 5; % factor to thin out the 

num_pps_lvls = unique(t.num_pps);
r = t.selection & t.num_pps == max(num_pps_lvls);
js = t(r,:);

%% Create the ecdf for the random vector:
p_rand = rand(size(js.PScaled,2)) - 0.5; % zero-centered
[~,~,rand_x] = get_ecdf_dot(p_rand,num_case,ecdf_save);
h = figure('units','inches','position',[1 1 20 10],'Name','Orthogonality of Parameters, CDF');
plot(rand_x,ecdf_save,'-k','LineWidth',3);
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
        [~,~,meth_x(:,i_iter)] = get_ecdf_dot(p_meth,num_case,ecdf_save);
    end
    
    meth_x_mean(:,i_met) = mean(meth_x,2);
    meth_x_std(:,i_met) = std(meth_x,[],2);
    ls = strcat(line_style{i_met},marker_style{i_met});
    
    plot(meth_x_mean(:,i_met),...
        ecdf_save,line_style{i_met},...
        'Color',plot_color(i_met,:),'LineWidth',3);
    
    g(i_met) = plot(meth_x_mean(std_thin:std_thin:numel(ecdf_save),i_met),...
        ecdf_save(std_thin:std_thin:numel(ecdf_save)),marker_style{i_met},...
        'Color',plot_color(i_met,:),'LineWidth',3,'MarkerSize',20);
    
    for k = std_thin:std_thin:numel(ecdf_save)
        ytmp = [ecdf_save(k) ecdf_save(k)];
        xtmp = [meth_x_mean(k,i_met)-meth_x_std(k,i_met) meth_x_mean(k,i_met)+meth_x_std(k,i_met)];
        ytmp_hat = [ytmp(1)-0.01 ytmp(2)+0.01];
        plot(xtmp,ytmp,'-','Color',plot_color(i_met,:),'LineWidth',3);
        plot([xtmp(1) xtmp(1)],ytmp_hat,'-','Color',plot_color(i_met,:),'LineWidth',3);
        plot([xtmp(2) xtmp(2)],ytmp_hat,'-','Color',plot_color(i_met,:),'LineWidth',3);
    end  
end % for i_met

figure(h);
% legend(g,g_names,'Location','NorthWest');
% legend boxoff;
xlabel('Dot Product of Scaled VP Parameters');
ylabel('Cumulative Distribution');
ylim([0 1.05]);
set(gca,'Ytick',0:0.2:1);
set(gca,'FontSize',28,'LineWidth',2);
legend(g,method_order,'Location','NorthWest');
legend boxoff;

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