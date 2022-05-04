function [Am,As] = mean_standard_deviation(A,LMP)
% Each row of A contains a model, each element of LMP contains log
% probability
MP=exp(LMP-max(LMP)); % ratio probability of each model/max probability
MP=MP/sum(MP); % calculate weight for each model
Am=sum(A.*repmat(MP',1,size(A,2))); % weighted mean
As=sqrt(sum(((A-repmat(Am,size(A,1),1)).^2).*repmat(MP',1,size(A,2)))); % SD

end