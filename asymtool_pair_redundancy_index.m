function [redundancy_index_list_bilateral,redundancy_index_list_all] = asymtool_pair_redundancy_index(A,LRU,is_bilateral,is_all)
%   Pair-Redundancy Index
%
%   [redundancy_index_list_bilateral,redundancy_index_list_all] = 
%   asymtool_pair_redundancy_index(A,LRU,is_bilateral,is_all) calculates
%   the pair-redundancy index of A. The pair redundancy of each pair 
%   corresponds to the difference between change of global efficiency
%   when both nodes are removed and sum of each nodes are removed, 
%   normalized by the sum of each nodes are removed.
%
%   Inputs:     A,      nxn Adjacent matrix where n is the number of nodes.
%                       adjacency matrix for a graph where nodes are
%                       organized in a bilateral symmetric manner, with
%                       each left nodes in the front or right nodes.
%               LRU,    nx3 matrix indicating states of (L) left and (R) 
%                       right for bilaterally symmetric neurons and (U)
%                       unilateral.
%               is_bilateral,    false : no / true : yes (default)
%               is_all,          false : no / true : yes (default)
%
%   Outputs:    redundancy_index_list_bilateral, 
%               redundancy_index_list_all,
%
%   Other m-files required: distance_bin from BCT
%   Subfunctions: none
%   MAT-files required: none
%
%   ____________________________________________________________________
%

if nargin < 2
    error('Not enough input arguments!!');
elseif nargin < 4
    is_bilateral = true;
    is_all = true;
end

redundancy_index_list_bilateral = [];
redundancy_index_list_all = [];
N = size(A,1);
    
left_idx_list = find(LRU(:,1));
right_idx_list = find(LRU(:,2));

% calculate original global efficiency
D = distance_bin(A);
D_inv = 1./D;
D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
GE = mean(D_inv_omitInf,'all');

if is_bilateral
    redundancy_index_list_bilateral = zeros(length(left_idx_list),1);
    for n_bi = 1:length(left_idx_list)
        left_idx = left_idx_list(n_bi);
        right_idx = right_idx_list(n_bi);

        % % % remove left node % % %
        A_temp = A;
        A_temp(left_idx,:) = 0; A_temp(:,left_idx) = 0;

        D = distance_bin(A_temp);
        D_inv = 1./D;
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE_removeLeft = mean(D_inv_omitInf,'all');

        % % % remove right node % % %
        A_temp = A;
        A_temp(right_idx,:) = 0; A_temp(:,right_idx) = 0;

        D = distance_bin(A_temp);
        D_inv = 1./D;
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE_removeRight = mean(D_inv_omitInf,'all');

        % % % remove both nodes % % %
        A_temp = A;
        A_temp(left_idx,:) = 0; A_temp(:,left_idx) = 0;
        A_temp(right_idx,:) = 0; A_temp(:,right_idx) = 0;

        D = distance_bin(A_temp);
        D_inv = 1./D;
        D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
        GE_removeBoth = mean(D_inv_omitInf,'all');

        redundancy_index_list_bilateral(n_bi) = ...
            ((GE-GE_removeBoth) - (GE-GE_removeLeft) - (GE-GE_removeRight))/abs(2*GE-GE_removeLeft-GE_removeRight);
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
    
    redundancy_index_list_all = zeros(N*(N-1)/2,1);
    n_temp = 0;
    for left_idx = 1:N
        for right_idx = left_idx+1:N
            % % % remove both nodes % % %
            A_temp = A;
            A_temp(left_idx,:) = 0; A_temp(:,left_idx) = 0;
            A_temp(right_idx,:) = 0; A_temp(:,right_idx) = 0;

            D = distance_bin(A_temp);
            D_inv = 1./D;
            D_inv_omitInf = D_inv(:); D_inv_omitInf(isinf(D_inv_omitInf))=[];
            GE_removeBoth = mean(D_inv_omitInf,'all');

            n_temp = n_temp + 1;
            redundancy_index_list_all(n_temp) = ...
            ((GE-GE_removeBoth) - vulnerability_list_all(left_idx) - vulnerability_list_all(right_idx))/abs(vulnerability_list_all(left_idx)+vulnerability_list_all(right_idx));
        end
    end
end

redundancy_index_list_bilateral(isnan(redundancy_index_list_bilateral)) = 1;
redundancy_index_list_all(isnan(redundancy_index_list_all)) = 1;

end

