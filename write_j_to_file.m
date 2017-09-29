function jtab = write_j_to_file(fname,j,iter)
%% Writes the main output structure to csv file
%

ntab = numel(j.sm.selection);

iter_col = iter*ones(ntab,1);
method_col = categorical(repmat({j.method},ntab,1));
numpps_col = j.num_pps*ones(ntab,1); % writing down the NOMINAL value
gof_col = j.sm.gof*ones(ntab,1);
betaS_col = j.sm.betaS*ones(ntab,1);

mu_col = repmat(j.sm.mu,ntab,1);
sigma_col = repmat(j.sm.sigma(:)',ntab,1);
p_bnds_lower_col = repmat(j.sm.p_bnds(:,1)',ntab,1);
p_bnds_upper_col = repmat(j.sm.p_bnds(:,2)',ntab,1);
time_pp_col = j.sm.time_per_pp*ones(ntab,1);
time_vp_col = j.sm.time_per_vp*ones(ntab,1);

jtab = table(iter_col,method_col,numpps_col,j.sm.selection, ...
    gof_col,mu_col,sigma_col,betaS_col,j.sm.ProbInclude, ...
    p_bnds_lower_col,p_bnds_upper_col, j.sm.Pfinal', ...
    j.sm.VPChar,time_pp_col,time_vp_col,...
    'VariableNames',{'iter','method','num_pps','selection',...
    'gof','mu','sigma','betaS','ProbInclude','p_lower',...
    'p_upper','Pfinal','VPChar','time_per_pp','time_per_vp'});

%% Find file name:
if iter < 10
    itername = strcat('000',num2str(iter));
elseif iter >= 10 && iter < 100
    itername = strcat('00',num2str(iter));
elseif iter >= 100 && iter < 1000
    itername = strcat('0',num2str(iter));
else
    itername = num2str(iter);
end  
filename = strcat(fname,itername,'.csv');

writetable(jtab,filename);

end