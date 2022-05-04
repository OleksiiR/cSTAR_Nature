function [lp] = Lp_AGivenY_new(A,A0,Vb_1,Vs,as,bs,delta)
% POSTERIOR LOG-PROBABILITY OF A GIVEN R => P(Ai|R)
% Calculates the posterior log-probability of the proposed sub-network
% given the global responses
% A0 : prior network
% A : proposed sub-network
% delta : penalization constant of Hamming distance
%[log(det(Vs))+log(det(Vb_1)) log(n/(n+1))]
lp=0.5*log(det(Vs))+0.5*log(det(Vb_1))-as*log(bs)-delta*sum(abs((A0-A)));%+log(nchoosek(na,p))+log(beta(alpha1+p,beta1+na-p));
end
