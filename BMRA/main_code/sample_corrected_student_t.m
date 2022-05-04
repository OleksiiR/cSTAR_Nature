function [r,A,LPA,LPR] = sample_corrected_student_t( X,Y,A_init,r_init,a,b,lambda,alpha1,beta1,noit)
% Samples network connections and local responses directed towards a given node using
% random walk Metropolis sampling
r=zeros(size(X,2),noit); % initializes "iteration-local response" matrix (nodes x noit)
A=zeros(size(X,2),noit); % initializes "iteration-network topography" matrix (nodes x noit)
LPA=zeros(noit,1);
LPR=zeros(noit,1);
%index=find(A_init==1);
delta=1; % penalization constant of Hamming distance
%%initialize log posterior
NKMAX=min(sum(A_init),size(X,1)-1); % min(nb of prior connections, nb of perturbations)== max nb of interactions => also change if new interactions wanted
A0=A_init; % keep prior (Ai)
A_init=InitialA(A0,NKMAX); % propose initial network topology (random sub-network of Ai)
r_init=r_init'.*A_init; % propose initial state of local responses (based on structure of sub-network)
Xe=X(:,logical(A_init)); % select global responses of nodes that interact with node i based on initial sub-network
[c, p, n, V1, Vb_1, Vs_1, Vs, mus, as, bs, na] = Prepare(Xe,Y,A_init,a,b,lambda); % prepare expression of posterior
%Lp_A_Y=Lp_AGivenY(p, n,  Vb_1, Vs, as, bs,alpha1,beta1, na);
Lp_A_Y=Lp_AGivenY_new(A_init, A0, Vb_1, Vs, as, bs,delta); % calculate posterior log-probability of initial sub-network given the global responses
for i=1:noit % initiate MCMC sampling
    temp_A=proposeA1(A_init,A0,NKMAX); % generate new sub-network
    Xe_temp=X(:,logical(temp_A)); % update global responses (Rs of nodes that interact with node i according to proposed sub-network)
    [c_t, p_t, n_t, V1_t, Vb_1_t, Vs_1_t, Vs_t, mus_t, as_t, bs_t, na_t] = Prepare(Xe_temp,Y,temp_A,a,b,lambda); % calculate new expression of posterior
    Lp_A_Y_t=Lp_AGivenY_new(temp_A,A0,Vb_1_t, Vs_t, as_t, bs_t,delta); % calculate posterior log-probability of proposed network given the global responses
        % sub-network sampling:
        % compares posterior log-probabilities of proposed and current A (new
        % candidate is always accepted if its posterior probability is
        % higher than that of the previously accepted candidate)
    if(Lp_A_Y_t>Lp_A_Y) 
        A_init=temp_A;
        Lp_A_Y=Lp_A_Y_t;
        c=c_t;
        p=p_t;
        n=n_t;
        V1=V1_t;
        Vb_1=Vb_1_t;
        Vs_1=Vs_1_t;
        Vs=Vs_t;
        mus=mus_t;
        as=as_t;
        bs=bs_t;
        na=na_t;
    else % if posterior probability of new candidate is lower or equal to that of the previously accepted candidate => acceptance/rejection step
        if sample_mini([Lp_A_Y_t Lp_A_Y],1)==1
            A_init=temp_A;
            Lp_A_Y=Lp_A_Y_t;
            c=c_t;
            p=p_t;
            n=n_t;
            V1=V1_t;
            Vb_1=Vb_1_t;
            Vs_1=Vs_1_t;
            Vs=Vs_t;
            mus=mus_t;
            as=as_t;
            bs=bs_t;
            na=na_t;
        end
    end
	% connection coefficient sampling:
    A(:,i)=A_init;
    r1=zeros(size(r_init)); % initialize vector of local responses
    %sum(A_init)
    %sample_r_st(a,bs,p,n,Vs,mus)
    [rt,rp]=sample_r_st(a,bs,p,n,Vs,mus); % random walk Metropolis sampling (NO acceptance/rejection)
    r1(logical(A_init))=rt; % update vector of local responses (while correcting with current sub-network structure)
    LPA(i)=Lp_A_Y; % add posterior log-probability of current sub-network
    LPR(i)=rp; % add posterior probability of current connection coefficients
    r(:,i)=r1; % add current connection coefficients
end

end
