function [MFdiff_struc_list_bilateral,MFdiff_struc_list_all,MFdiff_func_list_bilateral,MFdiff_func_list_all] = asymtool_motif_fingerprint_difference(A,LRU,motif_node,is_bilateral,is_all,is_struc,is_func)

%   Motif-fingerprint Difference Index
%
%   asymtool_motif_fingerprint_difference calculates the structural and 
%   functional motif-fingerprint difference between two nodes (bilateral
%   pair or all pair). The motif fingerprints are conceptualized as
%   probabilities and difference were measured by Jenson-Shannon divergence
%    
%   Inputs:     A,      nxn Adjacent matrix where n is the number of nodes.
%                       adjacency matrix for a graph where nodes are
%                       organized in a bilateral symmetric manner, with
%                       each left nodes in the front or right nodes.
%               LRU,    nx3 matrix indicating states of (L) left and (R) 
%                       right for bilaterally symmetric neurons and (U)
%                       unilateral.
%               motif_node,      3 / 4 
%               is_bilateral,    false / true (default)
%               is_all,          false / true (default)
%               is_struc,        false / true (default)
%               is_func,         false / true (default)
%
%   Outputs:    MFdiff_struc_list_bilateral, 
%               MFdiff_struc_list_all,
%               MFdiff_func_list_bilateral,
%               MFdiff_func_list_all,
%
%   Other m-files required: motif3struct_bin, motif3funct_bin,
%                           motif4struct_bin, motif4funct_bin from BCT
%                           JSDiv
%   Subfunctions: none
%   MAT-files required: none
%
%   ____________________________________________________________________
%

if nargin < 3
    error('Not enough input arguments!!');
elseif nargin < 5
    is_bilateral = true;
    is_all = true;
    is_struc = true;
    is_func = true;
elseif nargin < 7
    is_struc = true;
    is_func = true;
end

if motif_node == 3
    motif_struc_cal_fun = @motif3struct_bin;
    motif_func_cal_fun = @motif3funct_bin;
elseif motif_node == 4
    motif_struc_cal_fun = @motif4struct_bin;
    motif_func_cal_fun = @motif4funct_bin;
else
    error('Only 3 or 4 node-motif can be entered!!');
end
    

MFdiff_struc_list_bilateral = [];
MFdiff_func_list_bilateral = [];
MFdiff_struc_list_all = [];
MFdiff_func_list_all = [];

N = size(A,1);
    
left_idx_list = find(LRU(:,1));
right_idx_list = find(LRU(:,2));

if is_struc
    [~,F]=motif_struc_cal_fun(A);
    F_prob = F + 1;
    F_prob = F_prob ./ sum(F_prob);
    
    if is_bilateral
        % % %
        MFdiff_struc_list_bilateral = zeros(length(left_idx_list),1);
        for n_bi = 1:length(left_idx_list)
            left_idx = left_idx_list(n_bi);
            right_idx = right_idx_list(n_bi);

            F_prob_left = F_prob(:,left_idx);
            F_prob_right = F_prob(:,right_idx);

            MFdiff_struc_list_bilateral(n_bi) = JS_div(F_prob_left,F_prob_right);
        end
    end
    
    if is_all
        % % %
        MFdiff_struc_list_all = zeros(N*(N-1)/2,1);
        n_temp = 0;
        for left_idx = 1:N
            for right_idx = left_idx+1:N

                F_prob_left = F_prob(:,left_idx);
                F_prob_right = F_prob(:,right_idx);

                n_temp = n_temp + 1;
                MFdiff_struc_list_all(n_temp) = JS_div(F_prob_left,F_prob_right);
            end
        end
    end
end

if is_func
    [~,F]=motif_func_cal_fun(A);
    F_prob = F + 1;
    F_prob = F_prob ./ sum(F_prob);
    
    if is_bilateral
        % % %
        MFdiff_func_list_bilateral = zeros(length(left_idx_list),1);
        for n_bi = 1:length(left_idx_list)
            left_idx = left_idx_list(n_bi);
            right_idx = right_idx_list(n_bi);

            F_prob_left = F_prob(:,left_idx);
            F_prob_right = F_prob(:,right_idx);

            MFdiff_func_list_bilateral(n_bi) = JS_div(F_prob_left,F_prob_right);
        end
    end
    
    if is_all
        % % %
        MFdiff_func_list_all = zeros(N*(N-1)/2,1);
        n_temp = 0;
        for left_idx = 1:N
            for right_idx = left_idx+1:N

                F_prob_left = F_prob(:,left_idx);
                F_prob_right = F_prob(:,right_idx);

                n_temp = n_temp + 1;
                MFdiff_func_list_all(n_temp) = JS_div(F_prob_left,F_prob_right);
            end
        end
    end
end



end

