function [Mirror_Symmetry,Mirror_Symmetry_null] = asymtool_mirror_symmetry(Adj,cont_ind)
%   Mirror Symmetric Index
%
%   [Mirror_Symmetry,Mirror_Symmetry_null] = asymtool_mirror_symmetry(Adj, 
%   cont_ind) calculates the mirror symmetric index of Adj. Mirror
%   symmetric index 
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
%               is_bilateral,    0 : no / 1 : yes (default)
%               is_all,          0 : no / 1 : yes (default)
%
%   Outputs:    redundancy_index_list_bilateral, 
%               redundancy_index_list_all,
%
%   Other m-files required: distance_bin from BCT
%   Subfunctions: none
%   MAT-files required: none
Adj_sym=Adj.*Adj(cont_ind,cont_ind);
Mirror_Symmetry=sum(sum(Adj_sym))/sum(sum(Adj));

null_n=1000;
Mirror_Symmetry_null=zeros(null_n,1);

parfor i=1:null_n
    null_Adj=randmio_dir(Adj,1000);
    null_Adj_sym=null_Adj.*null_Adj(cont_ind,cont_ind);
    Mirror_Symmetry_null(i) = sum(sum(null_Adj_sym))/sum(sum(null_Adj));
end

end

