function D = KL_div(p,q)
    temp_p = p .* log(p./q);
    temp_p(p==0) = 0;
    D = sum(temp_p);
end