%% set general parameters

% specify names of modules in network

proteins = {'TGFBR1','EGFR','VGFR/PDGFR/FGFR','MEK/ERK','PI3K/mTOR','PKC','p38','JNK','Aurora-A','RIPK1','IKK2','JAK','GSK3B','DPD'};

% hyperparameters
nsamples = 200000; % number of MCMC sampling iterations
OK = 5000; % Occams window
repeat = 1; % Number of repeat runs

row=13; 
col=13;

%% OVCA420 network inference: TNF perturbation (complete)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TNF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
OVCA420_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_OVCA420_TNF.txt');
OVCA420_R = readmatrix('EMT_data/R_OVCA420_TNF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(OVCA420_R,OVCA420_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/OVCA420_Am_log_200_5K_TNF.csv');
writetable(As,'EMT_results/OVCA420_As_log_200_5K_TNF.csv');
writetable(rm,'EMT_results/OVCA420_rm_log_200_5K_TNF.csv');
writetable(rs,'EMT_results/OVCA420_rs_log_200_5K_TNF.csv');

%% OVCA420 network inference: TGFb perturbation (complete)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TGFb.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
OVCA420_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_OVCA420_TGFb.txt');
OVCA420_R = readmatrix('EMT_data/R_OVCA420_TGFb.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(OVCA420_R,OVCA420_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/OVCA420_Am_log_200_5K_TGFb.csv');
writetable(As,'EMT_results/OVCA420_As_log_200_5K_TGFb.csv');
writetable(rm,'EMT_results/OVCA420_rm_log_200_5K_TGFb.csv');
writetable(rs,'EMT_results/OVCA420_rs_log_200_5K_TGFb.csv');

%% OVCA420 network inference: EGF perturbation (PKCi missing)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_EGF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
OVCA420_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_OVCA420_EGF.txt');
OVCA420_R = readmatrix('EMT_data/R_OVCA420_EGF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(OVCA420_R,OVCA420_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/OVCA420_Am_log_200_5K_EGF.csv');
writetable(As,'EMT_results/OVCA420_As_log_200_5K_EGF.csv');
writetable(rm,'EMT_results/OVCA420_rm_log_200_5K_EGF.csv');
writetable(rs,'EMT_results/OVCA420_rs_log_200_5K_EGF.csv');

%% MCF7 network inference: TNF perturbation (PI3Ki, PKCi missing)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TNF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
MCF7_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_MCF7_TNF.txt');
MCF7_R = readmatrix('EMT_data/R_MCF7_TNF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(MCF7_R,MCF7_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/MCF7_Am_log_200_5K_TNF.csv');
writetable(As,'EMT_results/MCF7_As_log_200_5K_TNF.csv');
writetable(rm,'EMT_results/MCF7_rm_log_200_5K_TNF.csv');
writetable(rs,'EMT_results/MCF7_rs_log_200_5K_TNF.csv');

%% MCF7 network inference: TGFb perturbation (PI3Ki, Aurora-Ai, IKK2i missing)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TGFb.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
MCF7_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_MCF7_TGFb.txt');
MCF7_R = readmatrix('EMT_data/R_MCF7_TGFb.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(MCF7_R,MCF7_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/MCF7_Am_log_200_5K_TGFb.csv');
writetable(As,'EMT_results/MCF7_As_log_200_5K_TGFb.csv');
writetable(rm,'EMT_results/MCF7_rm_log_200_5K_TGFb.csv');
writetable(rs,'EMT_results/MCF7_rs_log_200_5K_TGFb.csv');

%% MCF7 network inference: EGF perturbation (complete)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_EGF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
MCF7_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_MCF7_EGF.txt');
MCF7_R = readmatrix('EMT_data/R_MCF7_EGF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(MCF7_R,MCF7_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/MCF7_Am_log_200_5K_EGF.csv');
writetable(As,'EMT_results/MCF7_As_log_200_5K_EGF.csv');
writetable(rm,'EMT_results/MCF7_rm_log_200_5K_EGF.csv');
writetable(rs,'EMT_results/MCF7_rs_log_200_5K_EGF.csv');

%% DU145 network inference: TNF perturbation (complete)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TNF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
DU145_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_DU145_TNF.txt');
DU145_R = readmatrix('EMT_data/R_DU145_TNF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(DU145_R,DU145_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/DU145_Am_log_200_5K_TNF.csv');
writetable(As,'EMT_results/DU145_As_log_200_5K_TNF.csv');
writetable(rm,'EMT_results/DU145_rm_log_200_5K_TNF.csv');
writetable(rs,'EMT_results/DU145_rs_log_200_5K_TNF.csv');

%% DU145 network inference: TGFb perturbation (complete)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TGFb.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
DU145_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_DU145_TGFb.txt');
DU145_R = readmatrix('EMT_data/R_DU145_TGFb.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(DU145_R,DU145_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/DU145_Am_log_200_5K_TGFb.csv');
writetable(As,'EMT_results/DU145_As_log_200_5K_TGFb.csv');
writetable(rm,'EMT_results/DU145_rm_log_200_5K_TGFb.csv');
writetable(rs,'EMT_results/DU145_rs_log_200_5K_TGFb.csv');

%% DU145 network inference: EGF perturbation (IKK2i missing)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_EGF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
DU145_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_DU145_EGF.txt');
DU145_R = readmatrix('EMT_data/R_DU145_EGF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(DU145_R,DU145_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/DU145_Am_log_200_5K_EGF.csv');
writetable(As,'EMT_results/DU145_As_log_200_5K_EGF.csv');
writetable(rm,'EMT_results/DU145_rm_log_200_5K_EGF.csv');
writetable(rs,'EMT_results/DU145_rs_log_200_5K_EGF.csv');

%% A549 network inference (lung cancer cell line): TNF perturbation (EGFRi (partially), VGFRi, PI3Ki, PKCi, Aurora-Ai, IKK2i missing)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TNF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
A549_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_A549_TNF.txt');
A549_R = readmatrix('EMT_data/R_A549_TNF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(A549_R,A549_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/A549_Am_log_200_5K_TNF.csv');
writetable(As,'EMT_results/A549_As_log_200_5K_TNF.csv');
writetable(rm,'EMT_results/A549_rm_log_200_5K_TNF.csv');
writetable(rs,'EMT_results/A549_rs_log_200_5K_TNF.csv');

%% A549 network inference (lung cancer cell line): TGFb perturbation (PI3K inhibitor missing (partially))

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_TGFb.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
A549_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_A549_TGFb.txt');
A549_R = readmatrix('EMT_data/R_A549_TGFb.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(A549_R,A549_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/A549_Am_log_200_5K_TGFb.csv');
writetable(As,'EMT_results/A549_As_log_200_5K_TGFb.csv');
writetable(rm,'EMT_results/A549_rm_log_200_5K_TGFb.csv');
writetable(rs,'EMT_results/A549_rs_log_200_5K_TGFb.csv');

%% A549 network inference (lung cancer cell line): EGF perturbation (complete)

% load data 
G = (readmatrix('EMT_data/EMT_scRNA_prior_network_EGF.txt'));
% G = (readmatrix('EMT_data/EMT_scRNA_prior_network.txt'));
A549_pert = readmatrix('EMT_data/EMT_scRNA_perturbation_matrix_A549_EGF.txt');
A549_R = readmatrix('EMT_data/R_A549_EGF.txt');

% network inference
cd main_code;
[Am,As,rm,rs,rtot]=ExecuteMultipleTimes(A549_R,A549_pert,G,nsamples,OK,repeat,row,col);
cd ..;

% results
Am = array2table(Am,'RowNames',proteins,'VariableNames',proteins);
As = array2table(As,'RowNames',proteins,'VariableNames',proteins);
rm = array2table(rm,'RowNames',proteins,'VariableNames',proteins);
rs = array2table(rs,'RowNames',proteins,'VariableNames',proteins);

writetable(Am,'EMT_results/A549_Am_log_200_5K_EGF.csv');
writetable(As,'EMT_results/A549_As_log_200_5K_EGF.csv');
writetable(rm,'EMT_results/A549_rm_log_200_5K_EGF.csv');
writetable(rs,'EMT_results/A549_rs_log_200_5K_EGF.csv');
