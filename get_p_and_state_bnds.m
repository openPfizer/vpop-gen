function [p_bnds,state_bnds] = ...
    get_p_and_state_bnds(mdl,num_state,p_values,search_mag)
%% function get_p_and_state_bnds
%   Calculate the lower and upper bounds on the fitting parameters and
%   states of the model. This function is entirely model dependent and must
%   be modified for each individual model.
%
%       Inputs:
%           mdl - SimBiology model object
%           names - Names of the model states, extracted from the model
%               object
%           p_values - Baseline values of the parameters to be changed
%           search_mag - scalar for parameter with poor literature bounds.
%
%       Outputs:
%           p_lower - Parameter lower bounds for optimization.
%           p_upper - Parameter upper bounds for optimization.
%           state_lower - Observable state lower bounds.
%           state_upper - Observable state upper bounds.
%
% Note on calculations: followed through analysis and calculations from
% Van de Pas 2012 with mean +/- 3*sigma (from the experimental data).
%
% NOTE(1): THIS FUNCTION IS HIGHLY MODEL-SPECIFIC AND MUST BE MODIFIED FOR OTHER
% MODELS.
%
% NOTE(2): This function is a modified version of the original function
% "Input_Ranges".
%

%% Initalize
p_bestfit = p_values;
%num_state = numel(state_names);
num_p = numel(p_values);

% Initialize arrays:
state_bestfit = -1*ones(num_state,1);
state_lower = state_bestfit;
state_upper = state_bestfit;
p_lower = -1*ones(num_p,1);
p_upper = -1*ones(num_p,1);

%% Set bounds on model states:
% Volumes of distribution:
v_intestine = mdl.compartments(1).capacity;
v_liver = mdl.compartments(2).capacity;
v_plasma = mdl.compartments(3).capacity;
v_periphery = mdl.compartments(4).capacity;

% intestine free cholesterol
state_bestfit(1) = 1.2736;
state_lower(1) = 0;
state_upper(1) = 4.8*v_intestine;

% intestesine cholesterol ester
state_bestfit(2) = 0.16;
state_lower(2) = 0;
state_upper(2) = 0.88*v_intestine;

% liver free cholesterol
state_bestfit(3) = 14.4;
state_lower(3) = 6.8*v_liver;
state_upper(3) = 9.2*v_liver;

% liver cholesterol ester
state_bestfit(4) = 9.54;
state_lower(4) = 3.5*v_liver;
state_upper(4) = 7.1*v_liver;

% Choose:
%%%  HDL_fc + HDL_ce < upperHDL = max(NHANES) = 4.0415{ to ensure
%%%  coverage of range}
% van de Pas claim 1:3 fc:ce ratio from Groener, et al. (1998) Atherosclerosis 137, 311-319
% assume this is accurate... then 4*HDL_fc <upperHDL & 1.33*HDL_ce
%%%

% plasma HDL free cholesterol
state_bestfit(5) = 0.837;
state_lower(5) = 0.3627/4;
state_upper(5) = v_plasma*4.0415/4;
% plasma HDL cholesterol ester
state_bestfit(6) = 2.4831;
state_lower(6) = 0.3627/1.333;
state_upper(6) = v_plasma*4.0415/1.333;

% plasma non-HDL cholesterol
state_bestfit(7) = 11.2437;
state_lower(7) = 0.5959*v_plasma;
state_upper(7) = 10.18*v_plasma;  %% NHAHNES

% peripherary cholesterol
state_bestfit(8) = 77.76;
state_lower(8) = 0.34*v_periphery;
state_upper(8) = 1.74*v_periphery;

% LDL
state_bestfit(9) = state_bestfit(7)*0.85;
state_lower(9) = 0.232*v_plasma; 
state_upper(9) = 8.571*v_plasma; %% NHAHNES

%% Set bounds on model parameters:
% paper/simbiology reaction name: "Hepatic Cholesterol Synthesis"
p_lower(1) = 0.0362; 
p_upper(1) = 1.44; %p_bestfit(1)=0.44;

% paper/simbiology reaction name: "Peripheral Cholesterol Synthesis"
p_lower(2) = 0; 
p_upper(2) = 12.1; %p_bestfit(2)=3.79;

% paper/simbiology reaction name: "Intestinal Cholesterol Synthesis"
p_lower(3) = 0.0692; 
p_upper(3) = 0.4995; %p_bestfit(3)=0.18;

% paper/simbiology reaction name: "Dietary Cholesterol Intake"
p_lower(4) = 0.01; 
p_upper(4) = 2.19; %p_bestfit(4)=1.09;

% paper/simbiology reaction name: "HDL-associated cholesterol esterification"
p_lower(9) = 0; 
p_upper(9) = 17.58/state_bestfit(5);    %p_bestfit(9)=9.39068;

% paper/simbiology reaction name: "biliary cholesterol excretion"
p_lower(14) = 1.28/state_bestfit(3); 
p_upper(14) = 6.03/state_bestfit(3); %p_bestfit(14)=0.2542;

% paper/simbiology reaction name: "fecal cholesterol excretion"
p_lower(15) = 0/state_bestfit(1); 
p_upper(15) = 3.71/state_bestfit(1);%%% range corrected for best fit steady state %p_bestfit(15)=1.4526;

% paper/simbiology reaction name: "intestinal cholesterol transport to HDL"
%p_lower(16)=0.0314; p_upper(16)=1.83/SS_bestfit(1);%p_bestfit(16)=0.0314;
%% Parameters with bounds that were more difficult to estimate:

% paper/simbiology reaction name: "hepatic cholesterol catabolism"
p_lower(18) = 0/state_bestfit(3); 
p_upper(18) = 2.23/state_bestfit(3);  %p_bestfit(18)=0.0715278;

% paper/simbiology reaction name: "CE transfer from HDL to LDL"
p_lower(21) = 3.46*v_plasma/(state_bestfit(6)*state_bestfit(7)); 
p_upper(21) = 6.64*v_plasma/(state_bestfit(6)*state_bestfit(7)); %p_bestfit(21)=0.49266;

% ratio of LDL to total non-HDL cholesterol (by definitiion between 0 and 1).
p_bestfit(22) = 0.85; 
p_lower(22) = 0;
p_upper(22) = 1;

%% Parameters with poorly defined literature ranges:
% Apply a uniform scalar to constrain.
p_2_scale = p_lower == -1 & p_upper == -1;
p_lower(p_2_scale) = p_bestfit(p_2_scale)/search_mag;
p_upper(p_2_scale) = p_bestfit(p_2_scale)*search_mag;

%% Package outputs:
p_bnds = [p_lower(:) p_upper(:)];
state_bnds = [state_lower(:) state_upper(:)];

end

