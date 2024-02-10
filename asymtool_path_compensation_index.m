function [compensation_index_list_bilateral,compensation_index_list_all] = asymtool_path_compensation_index(A,LRUL,is_bilateral,is_all)
% Calculate path-compensation index
%   A: (N by N) adjacency matrix (either directed or undirected)

if nargin < 2
    error('Not enough input arguments!!');
elseif nargin < 4
    is_bilateral = true;
    is_all = true;
end


compensation_index_list_bilateral = [];
compensation_index_list_all = [];
N = size(A,1);
    
left_idx_list = find(LRUL(:,1));
right_idx_list = find(LRUL(:,2));

if is_bilateral
    compensation_index_list_bilateral = zeros(length(left_idx_list),1);
    for n_bi = 1:length(left_idx_list)
        left_idx = left_idx_list(n_bi);
        right_idx = right_idx_list(n_bi);
        
        % calculate original global efficiency
        D = distance_bin(A);
        D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
        D_inv = 1./D;
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE = mean(D_inv_omitInf,'all');

        % % % remove left node % % %
        A_temp = A;
        A_temp(left_idx,:) = 0; A_temp(:,left_idx) = 0;

        D = distance_bin(A_temp);
        D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
        D_inv = 1./D; 
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE_removeLeft = mean(D_inv_omitInf,'all');

        % % % remove right node % % %
        A_temp = A;
        A_temp(right_idx,:) = 0; A_temp(:,right_idx) = 0;

        D = distance_bin(A_temp);
        D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
        D_inv = 1./D;
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE_removeRight = mean(D_inv_omitInf,'all');

        % % % remove both nodes % % %
        A_temp = A;
        A_temp(left_idx,:) = 0; A_temp(:,left_idx) = 0;
        A_temp(right_idx,:) = 0; A_temp(:,right_idx) = 0;

        D = distance_bin(A_temp);
        D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
        D_inv = 1./D;
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE_removeBoth = mean(D_inv_omitInf,'all');
        
        if (GE-GE_removeBoth) > 0
            compensation_index_list_bilateral(n_bi) = ...
                (max((GE-GE_removeBoth),0) - max(GE-GE_removeLeft,0) - max(GE-GE_removeRight,0))/max(GE-GE_removeBoth,0);
        else
            compensation_index_list_bilateral(n_bi) = 0;
        end
    end
end

if is_all
    vulnerability_list_all = zeros(N,1);
    for target_idx = 1:N
        % % % remove target node % % %
        A_temp = A;
        A_temp(target_idx,:) = 0; A_temp(:,target_idx) = 0;

        D = distance_bin(A_temp);
        D_inv = 1./D;
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE_removeTarget = mean(D_inv_omitInf,'all');
        vulnerability_list_all(target_idx) = GE - GE_removeTarget;
    end
    
    compensation_index_list_all = zeros(N*(N-1)/2,1);
    n_temp = 0;
    for left_idx = 1:N
%         disp(left_idx);
        for right_idx = left_idx+1:N
            % calculate original global efficiency
            D = distance_bin(A);
            D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
            D_inv = 1./D;
            D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
            GE = mean(D_inv_omitInf,'all');
            
            % % % remove left node % % %
            A_temp = A;
            A_temp(left_idx,:) = 0; A_temp(:,left_idx) = 0;

            D = distance_bin(A_temp);
            D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
            D_inv = 1./D; 
            D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
            GE_removeLeft = mean(D_inv_omitInf,'all');

            % % % remove right node % % %
            A_temp = A;
            A_temp(right_idx,:) = 0; A_temp(:,right_idx) = 0;

            D = distance_bin(A_temp);
            D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
            D_inv = 1./D;
            D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
            GE_removeRight = mean(D_inv_omitInf,'all');
            
            % % % remove both nodes % % %
            A_temp = A;
            A_temp(left_idx,:) = 0; A_temp(:,left_idx) = 0;
            A_temp(right_idx,:) = 0; A_temp(:,right_idx) = 0;

            D = distance_bin(A_temp);
            D([left_idx,right_idx],:) = []; D(:,[left_idx,right_idx]) = [];
            D_inv = 1./D;
            D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
            GE_removeBoth = mean(D_inv_omitInf,'all');

            n_temp = n_temp + 1;
            if (GE-GE_removeBoth) > 0
                compensation_index_list_all(n_temp) = ...
                    (max((GE-GE_removeBoth),0) - max(GE-GE_removeLeft,0) - max(GE-GE_removeRight,0))/max(GE-GE_removeBoth,0);
            else
                compensation_index_list_all(n_temp) = 0;
            end
        end
    end
end

compensation_index_list_bilateral(isnan(compensation_index_list_bilateral)) = 0;
compensation_index_list_all(isnan(compensation_index_list_all)) = 0;

end

