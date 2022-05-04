function [c, p, n, V1, Vb_1, Vs_1, Vs, mus, as, bs, na] = Prepare(X,Y,A,a,b,lambda)
% MULTIVARIATE STUDENT DISTRIBUTION (proportional to target
% distribution/posterior)

c=length(Y); % number of perturbations
p=size(X,2); % number of interactions leading to i
n=length(Y); % number of perturbations
V1=X'*X; % R*R' (here R = global responses of nodes that interact with node i)

% calculate expression of mu1 (location parameter) :
Vb_1=(V1+lambda*eye(size(X,2)))/(c); 
Vs_1=V1+Vb_1; % sigma1 (scale matrix )
Vs=Vs_1\eye(size(Vs_1));
%size(Vs)
%size(X)
%size(Y)
%lambda
mus=Vs*X'*Y; % mu (location parameter)

as=a+n/2; % a1 (degrees of freedom)
bs=b+0.5*(Y'*Y-mus'*Vs_1*mus);

na=length(A)-1; % number of nodes - 1

end
