function D = JS_div(p,q)
    m = 0.5 * (p + q);
    D = 0.5 * KL_div(p,m) + 0.5 * KL_div(q,m);
end