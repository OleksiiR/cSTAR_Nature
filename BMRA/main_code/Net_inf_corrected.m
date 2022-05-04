function [As1,Rs1,LMP1,r] = Net_inf_corrected(data,pert,G,a,b,lambda,alpha1,beta1,noit,OK)
% NETWORK INFERENCE
%data = number of measured nodes X number of perturbations (global response matrix)
%pert = perturbation matrix (indicates which node is perturbed)
Nodes=size(data,1);
As=zeros(Nodes,Nodes,noit); % initializes network topography cuboid
Rs=zeros(Nodes,Nodes,noit); % initializes local response cuboid
LPAs=zeros(Nodes,noit);
LPRs=zeros(Nodes,noit);

%mrb=zeros(Nodes);
%srb=zeros(Nodes);
%G=createInitialGraph1(size(data,1),size(data,2)-1);

rr=randn(size(data,1)).*G; % prior * gaussian noise

for i=1:Nodes % !!! original version: parfor !!!
    Ai=G(i,:); % vector of prior interactions leading to node i
    if(sum(Ai)>0) % test whether prior contains interactions leading to i (something to change if want to propose new interations, sampling skipped if prior network empty)
        
    ri=rr(i,:)';
    XY=data(:,logical(1-pert(i,:))); % extract global responses not associated with a perturbation of node i
    Y=XY(i,:)'; % global responses of node i
    X=XY';
    
    %initialize Ai and ri
    %make sure the number of nonzero connection coefficients does not
    %exceed the number of perturbation experiments
    if(sum(Ai)>length(Y)) % number of interacions > number of perturbations ? Maybe replace "if" with "while"?
        indx=find(Ai==1); % find index of all existing connections
        %indxe=ceil(rand(1,sum(Ai)-length(Y))*length(indx));
        
        %if you have too many interactions randomly delete the extra edges
        I1=randperm(length(indx)); % randomly permutate integer values up to the number of connections
        d1=sum(Ai)-length(Y); % determine number of excess connections to delete
        
        ri(I1(d1))=0; % replace connections with 0
        Ai(I1(d1))=0;
    end
   % Ai
    %start MCMC sampling
    [r, A, LPA, LPR]=sample_corrected_student_t(X,Y,Ai,ri,a,b,lambda,alpha1,beta1,noit);
    %store data
    %mr=mean(r(burnin+1:noit,:));
    %sr=std(r(burnin+1:noit,:));
    %pi=mean(A(burnin+1:noit,:));
    %P(i,:)=pi;
    %mrb(i,:)=mr;
    %srb(i,:)=sr;
    As(i,:,:)=A; % updates network topology cuboid for node/layer i using MCMC samples
    Rs(i,:,:)=r; % updates local response cuboid for node/layer i using MCMC samples
    LPAs(i,:)=LPA;
    LPRs(i,:)=LPR;
    end
end
%size(As)

r=Rs; % retain all iterations for r to plot posterior distributions

[As1,Rs1,LMP1]=Occam(As,Rs,LPAs,LPRs,OK); % select OK models with highest posterior probability
end
