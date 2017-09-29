function [state, options, optchanged] = ga_output(options, state, ~, num_pps)
%% ga_output.m : This is function that can be passed as an option to 
% Matlab's genetic algorthim routine, inputs and outputs are defined at
% http://www.mathworks.com/help/gads/genetic-algorithm-options.html#f17837

first_flag = 0;
optchanged = false;

if state.Generation == 0
    if exist('results/ga_interim.mat','file') > 0
        delete('results/ga_interim.mat');
    end
end

%% Find plausible patients, i.e. those with a score < 0
pp_idx = find(state.Score < 0);

if ~isempty(pp_idx)
    pps_new = state.Population(pp_idx,:);
    pps_scores_new = state.Score(pp_idx);
else
    %disp(state.Score);
    return;
end

%% Append to existing list of plausible patients if possible
if exist('results/ga_interim.mat','file') == 2
    % Generation 2+, load and append:
    
    try
        load('results/ga_interim.mat');
    catch
        disp(first_flag);
        disp(state.Generation);
        error('no load');
    end
    
    p_pp = [p_pp; 10.^pps_new];
    score_pp = [score_pp; pps_scores_new];
    
    % Filter out close-ish duplicates:
    z = p_pp./mean(p_pp,1);
    [~,ia,~] = uniquetol(z,1e-2,'ByRows',true);
    p_pp = p_pp(ia,:);
    score_pp = score_pp(ia,:);
    
else   %create new list of plausible patients
    p_pp = 10.^pps_new;
    score_pp = pps_scores_new;  
end

n_pp = size(p_pp,1);

%% Save results back to temporary file
save('txtout/ga_interim.mat', 'p_pp','score_pp','n_pp');
%disp(['Number of plausible patients found: ', num2str(n_pp)]);

%% Flag exit if we have enough PPs
if (n_pp >= num_pps)
    state.StopFlag = 'y';   
end

end % function ga_output
