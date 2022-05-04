function [r,p] = sample_r_st(a,bs,p,n,Vs,mus)
% LOCAL RESPONSE VECTOR SAMPLING
% Samples vector of local responses from a multivariate Student distribution using random walk Metropolis
as=a+(n+p)/2; % degrees of freedom (1 + (nb perturbations + nb of interactions leading to i)/2)
%as=a+n/2;
sig=(bs/as)*Vs; % sigma (variance matrix)
%mus
%mvtrnd(sig,as,1);
r1=mvtrnd(sig,as,1); % random walk: (1) samples from MVSt according to updated parameters
r=mus'+r1; % (2) updates location
p=mvtpdf(r1',sig,as); % returns probability density for this set of r
end
