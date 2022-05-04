function [A1] = proposeA1(A,A0,NKMAX)
% PROPOSAL FUNCTION FOR A 
% Proposes a new candidate sub-network based on the prior and
% uniform sampling
%A1=A;
if(sum(A)==NKMAX) % nb of interactions == max nb of interactions
   A1=remove(A); % remove an interaction
else
    if(sum(A)==0)
        A1=add(A,A0); % add an interaction (only adds interactions that exist in prior)
    else
        toss=rand();
        if(toss<0.5)
          A1=add(A,A0);
        else
           A1=remove(A);
        end
    end
end
    
    function [B] =add(A,A0) % randomly adds an interaction
        B=A; % keep current sub-network
        i=find((A0-A)==1); % find indices of interactions that are in prior but not in current sub-network
        j=ceil(rand()*length(i)); % determine the index of the interaction to be added
        B(i(j))=1; % add interaction to sub-network
    end

    function [B] =remove(A) % randomly removes an interaction
        B=A; % keep current sub-network  
        i=find(A==1); % find indices of interactions in current sub-network
        j=ceil(rand()*length(i)); % determine the index of connection to be removed
        B(i(j))=0; % remove connection
    end

end

%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% A1=A;
% if sum(A)==min(length(A)-1,n-1)
%     index=find(A==1);
%     j=ceil(rand()*length(index));
%     A1(index(j))=0;
% end
% 
% if sum(A)==0
%     index=setdiff(1:length(A),i);
%     j=ceil(rand()*length(index));
%     A1(index(j))=1;
% end
% 
% if sum(A)>=1 && sum(A)<min(length(A)-1,n-1)
%     toss=rand();
%     if(toss<0.5)
%         index=find(A==1);
%         j=ceil(rand()*length(index));
%         A1(index(j))=0;
%     else
%         index=setdiff(find(A==0),i);
%         j=ceil(rand()*length(index));
%         A1(index(j))=1;
%     end
% end
%     
% end
