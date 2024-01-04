function D = KL_div(p,q)

%   Kullback–Leibler divergence
%
%   Inputs:     p,q     proabability distributions.
%
%   Outputs:    D, 
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   ____________________________________________________________________
%
temp_p = p .* log(p./q);
temp_p(p==0) = 0;
D = sum(temp_p);
end