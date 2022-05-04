function [M1] = Standarize(M)
%% Standarize each Column
% M1=(M-repmat(mean(M),size(M,1),1))./repmat(std(M),size(M,1),1);
M1=(M-repmat(mean(M),size(M,1),1))./repmat(std(M),size(M,1),1);
end
