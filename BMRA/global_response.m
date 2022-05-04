function [R] = global_response(xt,xs)
% GLOBAL RESPONSE MATRIX
% Calculates the standardized global response matrix from the concentrations at t0 and t1
R=2*(xt-xs)./(xt+xs);
%R=R*50;
%R=log2(xt)-log2(xs);
%R=2*(log2(xt)-log2(xs))./(log2(xt)+log2(xs));
%R=log2(xt./xs);
%R=(R-mean(R(:)))./std(R(:));
R=Standarize(R')';
%R=R./repmat(std(R),size(R,1),1);
%R=R./repmat(std(R,0,2),1,size(R,2));
end
