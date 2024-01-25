function [p_val] = custom_shuffle_test(x_all,x_bi,n,alllow)

%   custom_shuffle_test
%
%   [p_val] = custom_shuffle_test(x_all,x_bi,n,biishigh) performs a 
%             permutation test to compare the means of two groups, all 
%             and bilateral, by randomly sampling from all and comparing it 
%             to bilateral in a one-way manner. The parameter alllow determines 
%             whether the mean of all values is considered lower than the
%             mean of the bilateral.
%
%   Inputs:  x_all,     all data
%            x_bi,      bilateral data
%            n,         iteration amount
%            alllow,    mean of all is lower than mean of bi
%
%   Outputs:    p_val, 
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
% 
%   ____________________________________________________________________
%

temp=zeros(1,n);
parfor i=1:n
    [~,~,~,stats] = ttest2(x_bi,x_all);
    t_obs = stats.tstat;
    
    perm_ind=randperm(length(x_all),length(x_bi));
    if alllow
        [~,~,~,stats] = ttest2(x_all(perm_ind),x_all);
        temp(i)= stats.tstat >= t_obs;
    else
        [~,~,~,stats] = ttest2(x_all(perm_ind),x_all);
        temp(i)=stats.tstat <= t_obs;
    end
end
p_val=sum(temp)/n;

end