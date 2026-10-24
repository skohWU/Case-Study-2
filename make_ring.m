function [G] = make_ring(n,k)

    G = zeros(n,n);
    S = circshift(eye(n),1);
    
    if k > n/2
        error('make_ring:k_too_large', 'k exceeds the maximum allowed value of n/2');
    else
        for j = 1:k
            G = G + S^j + (S^j)';
        end
    end

end

