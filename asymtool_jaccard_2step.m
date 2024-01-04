function [jaccard_index_in_out_list_bilateral,jaccard_index_in_out_list_all,jaccard_index_in_list_bilateral,jaccard_index_in_list_all] = asymtool_jaccard_2step(A,LRUL,is_bilateral,is_all,is_in,is_in_out)

%   2-step Jaccard Index
%
%   asymtool_jaccard_2step calculates the Jaccard similarity coefficient of
%   2step neighborhood. since the network is directed network, neighborhood
%   can be defined considering inward and outward directions. Here, we only
%   consider in neighborhood and in&out neighborhood.
%    
%   Inputs:     A,      nxn Adjacent matrix where n is the number of nodes.
%                       adjacency matrix for a graph where nodes are
%                       organized in a bilateral symmetric manner, with
%                       each left nodes in the front or right nodes.
%               LRU,    nx3 matrix indicating states of (L) left and (R) 
%                       right for bilaterally symmetric neurons and (U)
%                       unilateral.
%               is_bilateral,    false / true (default)
%               is_all,          false / true (default)
%               is_in,           false / true (default)
%               is_in_out,       false / true (default)
%
%   Outputs:    MFdiff_struc_list_bilateral, 
%               MFdiff_struc_list_all,
%               MFdiff_func_list_bilateral,
%               MFdiff_func_list_all,
%
%   Other m-files required: none                           
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
    is_in = true;
    is_in_out = true;
elseif nargin < 6
    is_in = true;
    is_in_out = true;
end
    

jaccard_index_in_out_list_bilateral = [];
jaccard_index_in_list_bilateral = [];
jaccard_index_in_out_list_all = [];
jaccard_index_in_list_all = [];

N = size(A,1);
    
left_idx_list = find(LRUL(:,1));
right_idx_list = find(LRUL(:,2));

% % % 
if is_bilateral
    if is_in_out
        jaccard_index_in_out_list_bilateral = zeros(length(left_idx_list),1);
    end
    if is_in
        jaccard_index_in_list_bilateral = zeros(length(left_idx_list),1);
    end
    for n_bi = 1:length(left_idx_list)
        left_idx = left_idx_list(n_bi);
        right_idx = right_idx_list(n_bi);

        conn_left_1_in = A(:,left_idx);
        conn_left_2_in = sum(A(:,logical(conn_left_1_in)),2)>=1;
        conn_right_1_in = A(:,right_idx);
        conn_right_2_in = sum(A(:,logical(conn_right_1_in)),2)>=1;

        conn_left_1_out = A(left_idx,:)';
        conn_left_2_out = sum(A(logical(conn_left_1_out),:)',2)>=1;
        conn_right_1_out = A(right_idx,:)';
        conn_right_2_out = sum(A(logical(conn_right_1_out),:)',2)>=1;

        if is_in_out
            conn_left = [(conn_left_1_in+conn_left_2_in)>=1;(conn_left_1_out+conn_left_2_out)>=1];
            conn_right = [(conn_right_1_in+conn_right_2_in)>=1;(conn_right_1_out+conn_right_2_out)>=1];
            
            jaccard_index_in_out_list_bilateral(n_bi) = sum(conn_left.*conn_right) / sum((conn_left+conn_right)>=1);
        end
        
        if is_in
            in_conn_left =(conn_left_1_in+conn_left_2_in)>=1;
            in_conn_right = (conn_right_1_in+conn_right_2_in)>=1;

            jaccard_index_in_list_bilateral(n_bi) = sum(in_conn_left.*in_conn_right) / sum((in_conn_left+in_conn_right)>=1);
        end
    end
end

% % %
if is_all
    if is_in_out
        jaccard_index_in_out_list_all = zeros(N*(N-1)/2,1);
    end
    if is_in
        jaccard_index_in_list_all = zeros(N*(N-1)/2,1);
    end
    n_temp = 0;
    for left_idx = 1:N
        for right_idx = left_idx+1:N
            
            n_temp = n_temp + 1;

            conn_left_1_in = A(:,left_idx);
            conn_left_2_in = sum(A(:,logical(conn_left_1_in)),2)>=1;
            conn_right_1_in = A(:,right_idx);
            conn_right_2_in = sum(A(:,logical(conn_right_1_in)),2)>=1;

            conn_left_1_out = A(left_idx,:)';
            conn_left_2_out = sum(A(logical(conn_left_1_out),:)',2)>=1;
            conn_right_1_out = A(right_idx,:)';
            conn_right_2_out = sum(A(logical(conn_right_1_out),:)',2)>=1;
            
            if is_in_out
                
                conn_left = [(conn_left_1_in+conn_left_2_in)>=1;(conn_left_1_out+conn_left_2_out)>=1];
                conn_right = [(conn_right_1_in+conn_right_2_in)>=1;(conn_right_1_out+conn_right_2_out)>=1];
                
                if sum((conn_left+conn_right)>=1) == 0
                    jaccard_index_in_out_list_all(n_temp) = nan;
                else
                    jaccard_index_in_out_list_all(n_temp) = sum(conn_left.*conn_right) / sum((conn_left+conn_right)>=1);
                end
            end

            if is_in
                in_conn_left =(conn_left_1_in+conn_left_2_in)>=1;
                in_conn_right = (conn_right_1_in+conn_right_2_in)>=1;

                if sum((in_conn_left+in_conn_right)>=1) == 0
                    jaccard_index_in_list_all(n_temp) = nan;
                else
                    jaccard_index_in_list_all(n_temp) = sum(in_conn_left.*in_conn_right) / sum((in_conn_left+in_conn_right)>=1);
                end
            end
        end
    end
end

jaccard_index_in_out_list_bilateral(isnan(jaccard_index_in_out_list_bilateral)) = 0;
jaccard_index_in_list_bilateral(isnan(jaccard_index_in_list_bilateral)) = 0;
jaccard_index_in_out_list_all(isnan(jaccard_index_in_out_list_all)) = 0;
jaccard_index_in_list_all(isnan(jaccard_index_in_list_all)) = 0;

end

