function [A] = InitialA(A_init,NKMAX)
% Initial sub-network 
% Proposes a random sub-network of the prior network
% A_init = prior network
% NKMAX = maximum number of connections allowed (limited by prior & number
% of perturbations)
index=find(A_init==1); % returns indices of prior interactions
%rand(1,NKMAX)
%length(index)
l=ceil(rand(1,NKMAX)*length(index)); % rounds up to next higher integer value of (vector of uniformly distributed random numbers of length NKMAX * number of prior interactions)
l=unique(l); % extract unique values of l (no repetitions)
A=zeros(size(A_init)); % initializes 0 matrix A of same dimensions as A_init
A(index(l))=1;  

end
