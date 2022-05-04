function [Am,As,rm,rs,rtot] = ExecuteMultipleTimes(data,pert,G,noit,OK,times,row,col)
%% Inputs:  
%%(a)data= (N x Np) Global response matrix; N= Number of nodes; Np = Number of perturbation experiments; 
%%(b)pert = N x Np binary matrix where pert(i,j)=1 indicates that data(:,j) is the global response of the network when node i is perturbed.
%%(c)noit = number of Gibbs scans; 
%%(d)times =number of times the program is intended to be execute
%% Outputs:
%%(a) PS= posterior edge probability matrix .

% prior parameters
a=1; % a & b: parameters of inverse Gamma distribution 
b=1; 
lambda=0.2; % parameter for the variance of gaussian prior of interaction strengths;
alpha1=1;
beta1=1;
%OK=10; % number of top models for Occams window

sd=size(data); % returns dimensions of global response matrix (number of nodes & number of perturbations)
AS=zeros(sd(1),sd(1),OK*times); % 3D matrix of 0s (nb of nodes x nb of nodes x OK*times)
RS=zeros(sd(1),sd(1),OK*times);
LMPs=zeros(sd(1),OK*times); % 2D matrix of 0s (nb of nodes x OK*times)

for i=1:times
    [As,Rs,LMP,r]=Net_inf_corrected(data,pert,G,a,b,lambda,alpha1,beta1,noit,OK);
    AS(:,:,(i-1)*(OK)+1:i*(OK))=As;
    RS(:,:,(i-1)*(OK)+1:i*(OK))=Rs;
    LMPs(:,(i-1)*(OK)+1:i*(OK))=LMP;
    rtot(:,:,(i-1)*noit+1:i*(noit))=r; % for plotting posterior distribution
end

Am=zeros(sd(1));
As=zeros(sd(1));
rm=zeros(sd(1));
rs=zeros(sd(1));

for i=1:sd(1) % iterate over number of nodes
    %size(As)
    A1=reshape(AS(i,:,:),sd(1),times*OK)';
    r1=reshape(RS(i,:,:),sd(1),times*OK)';
    LMP1=LMPs(i,:);
    [Ami, Asi]=mean_standard_deviation(A1,LMP1);
    [rmi, rsi]=mean_standard_deviation(r1,LMP1);
    Am(i,:)=Ami; % posterior mean of model
    As(i,:)=Asi; % posterior SD of model
    rm(i,:)=rmi; % posterior mean of connection coefficients
    rs(i,:)=rsi; % posterior SD of connection coefficients
end

rtot=rtot(row,col,:); % for plotting posterior distribution
rtot=(reshape(rtot,1,[]))';
end
