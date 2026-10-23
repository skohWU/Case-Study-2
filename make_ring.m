function [G] = make_ring(n,k)

G = zeros(n,n);
S = circshift(eye(n),1);

for j = 1:k

    G = G + S^j + (S^j)';
end

end

