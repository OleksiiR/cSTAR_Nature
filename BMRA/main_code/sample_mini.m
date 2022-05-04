function [i] = sample_mini(l,T)
% ACCEPTANCE-REJECTION FUNCTION FOR A (Random Walk Metropolis)
% Calculates the acceptance probability of the new sub-network candidate
% and samples a uniformly distributed random number in order to accept or
% reject the candidate network

% p0=exp(l(1));
% p1=exp(l(2));
% sum=p0+p1;
% P0=[p0/sum p1/sum];
l=l/T;
% f=max(l);
% logsum=f+log(sum(exp(l-f)));
% P=cumsum(exp(l-logsum));
% r=rand();
% i=min(find(P>=r));

% Calculate the acceptance probability a:
% exp(log of posterior probability of Anew - log of posterior probability
% of A current) == ratio of the posterior probabilities of A new and A current
a=exp(l(1)-l(2)); 
    % acceptance-rejection step:
if rand<a % if random uniformly distributed number < a we accept the new network, else we reject the new network
    i=1;
else
    i=2;
end
end