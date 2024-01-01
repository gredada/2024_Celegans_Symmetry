function [p_val] = custom_shuffle_test(x_all,x_bi,n,highlow)
%UNTITLED2 이 함수의 요약 설명 위치
%   자세한 설명 위치

temp=zeros(1,n);
parfor i=1:n
    perm_ind=randperm(length(x_all),length(x_bi));
    if highlow
        temp(i)=mean(x_all(perm_ind))>mean(x_bi);
    else
        temp(i)=mean(x_all(perm_ind))<mean(x_bi);
    end
end
p_val=sum(temp)/n;

end

