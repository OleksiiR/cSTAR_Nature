function [A1,r1,LMP1] = Occam(A,r,LPA,LPR,K)
% OCCAMS WINDOW
% Selects the K models with the highest posterior probabilities
% LPA = matrix of posterior sub-network probabilities
% LPR = matrix of posterior connection coefficient probabilities
% A=model;
% r=connection coefficient; K number of models, returns top K models

LMP=LPA+LPR; % log model posterior probability
LMP1=zeros(size(LMP,1),K); % initialize matrix of log-posterior probabilities of model
A1=zeros(size(A,1),size(A,2),K); % initializes new matrix for networks of selected models
r1=zeros(size(r,1),size(r,2),K); % initializes new matrix for connection coefficients of selected models
for i=1:size(A,1) % iterate over number of nodes
   MPi=LMP(i,:); % select model probabilites for node i
   % sort model probabilities in descending order
   % SMPi : probability of model
   % Is : index of model in original vector (MPi)
   [SMPi,Is]=sort(MPi,'descend'); 
   Is1=Is(1:K); % select first K indices (= K models with highest probability)
   %size(A,3)
   Ai=A(i,:,Is1); % select sub-networks of K highest probability models
   ri=r(i,:,Is1); % select connection coefficients of K highest probability models
   A1(i,:,:)=Ai; % update matrix of sub-networks
   r1(i,:,:)=ri; % update matrix of connection coefficients
   LMP1(i,:)=SMPi(1:K); % update matrix of posterior model probabilities with highest probability models
end

end