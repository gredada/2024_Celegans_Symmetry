function D = JS_div(p,q)

%   Jenson-Shannon diversity
%
%   Inputs:     p,q     proabability distributions.
%
%   Outputs:    D,
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: KL_div
%
%   ____________________________________________________________________
%
m = 0.5 * (p + q);
D = 0.5 * KL_div(p,m) + 0.5 * KL_div(q,m);
end