function [G] = make_ring(n,k)

G = zeros(n,n);

for j = 1:k
    S = circshift(eye(n),j);
    G = G + S + (S)';
end

end

