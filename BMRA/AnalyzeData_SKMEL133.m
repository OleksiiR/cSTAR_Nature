%% set general parameters

% specify names of modules in network

proteins = {'ERK','AKT','mTOR','SRC','CDK','PKC','IRS','MYC','DPD'};

% hyperparameters
nsamples = 200000; % number of MCMC sampling iterations
OK = 5000; % Occams window
repeat = 1; % Number of repeat runs

row=8; 
col=8;


%% SKMEL 133 network inference (melanoma cell line), with Myc

% load data 
G = (readmatrix('data_SKMEL_133/SKMEL_133_prior_network_withMyc.txt'));
SKMEL_133_pert = readmatrix('data_SKMEL_133/SKMEL_133_perturbation_matrix_withMyc.txt');
SKMEL_133_R = readmatrix('data_SKMEL_133/SKMEL_133_R_global_withMyc.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(SKMEL_133_R,SKMEL_133_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'results_SKMEL_133/SKMEL133_Am_log_200_5K_withMyc.csv');
writetable(As,'results_SKMEL_133/SKMEL133_As_log_200_5K_withMyc.csv');
writetable(rm,'results_SKMEL_133/SKMEL133_rm_log_200_5K_withMyc.csv');
writetable(rs,'results_SKMEL_133/SKMEL133_rs_log_200_5K_withMyc.csv');

%% SKMEL 133 network inference (melanoma cell line), no Myc

proteins = {'ERK','AKT','mTOR','SRC','CDK','PKC','IRS','DPD'};

% load data 
G = (readmatrix('data_SKMEL_133/SKMEL_133_prior_network_noMyc.txt'));
SKMEL_133_pert = readmatrix('data_SKMEL_133/SKMEL_133_perturbation_matrix_noMyc.txt');
SKMEL_133_R = readmatrix('data_SKMEL_133/SKMEL_133_R_global_noMyc.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(SKMEL_133_R,SKMEL_133_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'results_SKMEL_133/SKMEL133_Am_log_200_5K_noMyc.csv');
writetable(As,'results_SKMEL_133/SKMEL133_As_log_200_5K_noMyc.csv');
writetable(rm,'results_SKMEL_133/SKMEL133_rm_log_200_5K_noMyc.csv');
writetable(rs,'results_SKMEL_133/SKMEL133_rs_log_200_5K_noMyc.csv');
