function [Mirror_Symmetry,Mirror_Symmetry_null] = asymtool_mirror_symmetry(Adj,cont_ind)

%   Mirror Symmetric Index
%
%   [Mirror_Symmetry,Mirror_Symmetry_null] = asymtool_mirror_symmetry(Adj, 
%   cont_ind) calculates the mirror symmetric index of Adj. Mirror
%   symmetric index 
%
%   Inputs:  A,         nxn Adjacent matrix where n is the number of nodes.
%                       adjacency matrix for a graph where nodes are
%                       organized in a bilateral symmetric manner, with
%                       each left nodes in the front or right nodes.
%            cont_ind,  nx1 matrix indicating index of contralateral
%                       neuron.
%
%   Outputs:    Mirror_Symmetry, 
%               Mirror_Symmetry_null,
%
%   Other m-files required: randmio_dir from BCT
%   Subfunctions: none
%   MAT-files required: none
% 
%   ____________________________________________________________________
%

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

