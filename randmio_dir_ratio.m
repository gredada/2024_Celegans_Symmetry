function [A_shuffled, eff]=randmio_dir_ratio(A, cont_ind, p_asym, p_sym, ITER)
%RANDMIO_DIR    Random graph with preserved in/out degree distribution
%
%   R = randmio_dir(W, ITER);
%   [R eff] = randmio_dir(W, ITER);
%
%   This function randomizes a directed network, while preserving the in-
%   and out-degree distributions. In weighted networks, the function
%   preserves the out-strength but not the in-strength distributions.
%
%   Input:      W,      directed (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%               eff,    number of actual rewirings carried out
%
%   References: Maslov and Sneppen (2002) Science 296:910
%
%
%   2007-2012
%   Mika Rubinov, UNSW
%   Olaf Sporns, IU

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Mar 2012: Limit number of rewiring attempts, count number of successful
%             rewirings (Olaf Sporns)
n=size(A,1);

A_sym=A.*A(cont_ind,cont_ind);
A_asym=A-A_sym;

R = zeros(n);

R(randsample(find(A_sym(:)),floor(length(find(A_sym(:)))*p_sym)))=1;
R(randsample(find(A_asym(:)),floor(length(find(A_asym(:)))*p_asym)))=1;

A_rest=A-R;

[i j]=find(R);
K=length(i);
ITER=K*ITER;

% maximal number of rewiring attempts per 'iter'
maxAttempts= round(n*K/(n*(n-1)));
% actual number of successful rewirings
eff = 0;

for iter=1:ITER
    att=0;
    while (att<=maxAttempts)                                     %while not rewired
        while 1
            e1=ceil(K*rand);
            e2=ceil(K*rand);
            while (e2==e1),
                e2=ceil(K*rand);
            end
            a=i(e1); b=j(e1);
            c=i(e2); d=j(e2);

            if all(a~=[c d]) && all(b~=[c d]);
                break           %all four vertices must be different
            end
        end

        %rewiring condition
        if ~(A(a,d) || A(c,b))
            R(a,d)=R(a,b); R(a,b)=0;
            R(c,b)=R(c,d); R(c,d)=0;

            A(a,d)=A(a,b); A(a,b)=0;
            A(c,b)=A(c,d); A(c,d)=0;
            
            
            j(e1) = d;          %reassign edge indices
            j(e2) = b;
            eff = eff+1;
            break;
        end %rewiring condition
        att=att+1;
    end %while not rewired
    
    
end %iterations

A_shuffled=A_rest+R;