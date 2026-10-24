function G = make_disk(P, r)
    
    n = size(P,1);
    s = zeros(n,1);
    
    for i = 1:n
        s(i) = P(i,1)^2 + P(i,2)^2;
    end

    delta = s * ones(n,1)' + ones(n,1) * s' - 2 * (P * P');
    G = (delta <= r^2);
    G = G - eye(n,n);



end