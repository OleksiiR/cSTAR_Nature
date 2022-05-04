%% set general parameters

% specify names of modules in network
proteins={'TRK','ERK','AKT','JNK','S6K','RSK','RTK','STV'};
% proteins ={'ActiveCaspase3','AKT','EGFR','ERK','JNK','p21','p38','RSK','S6K'};
% proteins ={'ActiveCaspase3','AKT','EGFR','ERK','JNK','p21','p38','RAS','RSK','S6K'};

% hyperparameters
nsamples = 200000; % number of MCMC sampling iterations
OK = 5000; % Occams window (comment if configuring hyperparameters, else leave uncommented)
repeat = 1; % Number of repeat runs

% if plotting posterior distribution of connection coefficients, select
% corresponding coordinates
row=7; 
col=7;

%% TrkA network inference (with example data)

% 45 min load data
G=(readmatrix('data/prior_network_new.txt')); % prior network
trkA_pert = readmatrix('data/perturbation_matrix.txt'); % pertubation matrix
trkA_R = readmatrix('data/R_log_TrkA_45.txt'); % global responses

% execute network inference task
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(trkA_R,trkA_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% plot posterior distribution of connection coefficients
% histogram(rtot); 

% store & export results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'results/TrkA_Am_45_log_200_5K.csv');
writetable(As,'results/TrkA_As_45_log_200_5K.csv');
writetable(rm,'results/TrkA_rm_45_log_200_5K.csv');
writetable(rs,'results/TrkA_rs_45_log_200_5K.csv');


% 10 min load data
G=(readmatrix('data/prior_network_new.txt')); % prior network
trkA_pert = readmatrix('data/perturbation_matrix.txt'); % pertubation matrix
trkA_R = readmatrix('data/R_log_TrkA_10.txt'); % global responses

% execute network inference task
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(trkA_R,trkA_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% plot posterior distribution of connection coefficients
% histogram(rtot); 

% store & export results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'results/TrkA_Am_10_log_200_5K.csv');
writetable(As,'results/TrkA_As_10_log_200_5K.csv');
writetable(rm,'results/TrkA_rm_10_log_200_5K.csv');
writetable(rs,'results/TrkA_rs_10_log_200_5K.csv');



%% TrkB network inference (with example data)

% 45 min load data
G=(readmatrix('data/prior_network_new.txt')); % prior network
trkA_pert = readmatrix('data/perturbation_matrix.txt'); % pertubation matrix
trkA_R = readmatrix('data/R_log_TrkB_45.txt'); % global responses

% execute network inference task
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(trkA_R,trkA_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% plot posterior distribution of connection coefficients
% histogram(rtot); 

% store & export results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'results/TrkB_Am_45_log_200_5K.csv');
writetable(As,'results/TrkB_As_45_log_200_5K.csv');
writetable(rm,'results/TrkB_rm_45_log_200_5K.csv');
writetable(rs,'results/TrkB_rs_45_log_200_5K.csv');


% 10 min load data
G=(readmatrix('data/prior_network_new.txt')); % prior network
trkA_pert = readmatrix('data/perturbation_matrix.txt'); % pertubation matrix
trkA_R = readmatrix('data/R_log_TrkB_10.txt'); % global responses

% execute network inference task
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(trkA_R,trkA_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% plot posterior distribution of connection coefficients
% histogram(rtot); 

% store & export results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'results/TrkB_Am_10_log_200_5K.csv');
writetable(As,'results/TrkB_As_10_log_200_5K.csv');
writetable(rm,'results/TrkB_rm_10_log_200_5K.csv');
writetable(rs,'results/TrkB_rs_10_log_200_5K.csv');



%% TrkA network inference (with configuration of hyperparameters)

% % load data
% G=(readmatrix('data/prior_network_new.txt')); % prior network
% trkA_pert = readmatrix('data/perturbation_matrix.txt'); % pertubation matrix
% trkA_R = readmatrix('data/R_log_TrkA_45.txt'); % global responses
% 
% % begin hyperparameter testing (OK parameter)
% occ = [5000 10000 15000 20000 25000 30000 40000 50000]; % specify vecor of values to test
% j = 1;
% 
% AhyperPAm = zeros(length(proteins),length(proteins),length(occ));
% AhyperPAs = zeros(length(proteins),length(proteins),length(occ));
% AhyperPrm = zeros(length(proteins),length(proteins),length(occ));
% AhyperPrs = zeros(length(proteins),length(proteins),length(occ));
% 
% for occ = occ
% 
% % execute network inference task
% cd main_code;
% [Am,As,rm,rs,rtot]=ExecuteMultipleTimes(trkA_R,trkA_pert,G,nsamples,occ,repeat,row,col); 
% cd ..;
% 
% % store & export results
% % Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
% % As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
% % rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
% % rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);
% % 
% % writetable(Am,'results/TrkA_Am_45_log_200_5K.csv');
% % writetable(As,'results/TrkA_As_45_log_200_5K.csv');
% % writetable(rm,'results/TrkA_rm_45_log_200_5K.csv');
% % writetable(rs,'results/TrkA_rs_45_log_200_5K.csv');
% 
% AhyperPAm(:,:,j) = Am;
% AhyperPAs(:,:,j) = As;
% AhyperPrm(:,:,j) = rm;
% AhyperPrs(:,:,j) = rs;
% 
% j = j + 1;
% end
% 
% Aprop = zeros(size(AhyperPrm));
% 
% for i = 1:size(AhyperPrm,3)
%     
% Aprop(:,:,i) = AhyperPrs(:,:,i)./abs(AhyperPrm(:,:,i));
% 
% end

%% TrkB network inference

% % load data
% G=(readmatrix('data/prior_network_new.txt')); % prior network
% trkB_pert = readmatrix('data/perturbation_matrix.txt'); % pertubation matrix
% trkB_R = readmatrix('data/R_log_TrkB_45.txt'); % global responses
% 
% %prohibited = readmatrix('data/prohibited_interactions.txt');
% 
% % % hyperparameter testing
% % occ = [5000 10000 15000 20000 25000 30000 40000 50000]; % number of selected samples
% % j = 1;
% % 
% % BhyperPAm = zeros(length(proteins),length(proteins),length(occ));
% % BhyperPAs = zeros(length(proteins),length(proteins),length(occ));
% % BhyperPrm = zeros(length(proteins),length(proteins),length(occ));
% % BhyperPrs = zeros(length(proteins),length(proteins),length(occ));
% % 
% % for occ = occ
% 
% % execute network inference
% cd main_code;
% [Am,As,rm,rs,rtot]=ExecuteMultipleTimes(trkB_R,trkB_pert,G,nsamples,OK,repeat,row,col); % OV: nsamples/2 instead of OK
% cd ..;
% 
% %histogram(rtot);
% 
% %writematrix(rtot,'results/p_dist_Rtk_Erk_45min.csv');
% 
% % store results
% Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
% As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
% rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
% rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);
% 
% writetable(Am,'results/TrkB_Am_45_log_200_5K.csv');
% writetable(As,'results/TrkB_As_45_log_200_5K.csv');
% writetable(rm,'results/TrkB_rm_45_log_200_5K.csv');
% writetable(rs,'results/TrkB_rs_45_log_200_5K.csv');
% 
% % BhyperPAm(:,:,j) = Am;
% % BhyperPAs(:,:,j) = As;
% % BhyperPrm(:,:,j) = rm;
% % BhyperPrs(:,:,j) = rs;
% % 
% % j = j + 1;
% % end
% % 
% % 
% % Bprop = zeros(size(BhyperPrm));
% % 
% % for i = 1:size(BhyperPrm,3)
% %     
% % Bprop(:,:,i) = BhyperPrs(:,:,i)./abs(BhyperPrm(:,:,i));
% % 
% % end

