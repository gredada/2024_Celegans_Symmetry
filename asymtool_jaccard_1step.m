function [jaccard_index_in_out_list_bilateral,jaccard_index_in_out_list_all,jaccard_index_in_list_bilateral,jaccard_index_in_list_all] = asymtool_jaccard_1step(A,LRUL,is_bilateral,is_all,is_in,is_in_out)
% Calculate 2-step jaccard index
%   A: (N by N) adjacency matrix (either directed or undirected)

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
        conn_right_1_in = A(:,right_idx);

        conn_left_1_out = A(left_idx,:)';
        conn_right_1_out = A(right_idx,:)';

        if is_in_out
            conn_left = [(conn_left_1_in)>=1;(conn_left_1_out)>=1];
            conn_right = [(conn_right_1_in)>=1;(conn_right_1_out)>=1];
            
            jaccard_index_in_out_list_bilateral(n_bi) = sum(conn_left.*conn_right) / sum((conn_left+conn_right)>=1);
        end
        
        if is_in
            in_conn_left =(conn_left_1_in)>=1;
            in_conn_right = (conn_right_1_in)>=1;

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
            conn_right_1_in = A(:,right_idx);

            conn_left_1_out = A(left_idx,:)';
            conn_right_1_out = A(right_idx,:)';
            
            if is_in_out
                
                conn_left = [(conn_left_1_in)>=1;(conn_left_1_out)>=1];
                conn_right = [(conn_right_1_in)>=1;(conn_right_1_out)>=1];
                
                if sum((conn_left+conn_right)>=1) == 0
                    jaccard_index_in_out_list_all(n_temp) = nan;
                else
                    jaccard_index_in_out_list_all(n_temp) = sum(conn_left.*conn_right) / sum((conn_left+conn_right)>=1);
                end
            end

            if is_in
                in_conn_left =(conn_left_1_in)>=1;
                in_conn_right = (conn_right_1_in)>=1;

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

