function [W_A,W_B] = makeProtocols(d,G)
    D = diag(d);
    L = D - G;

    W_A =  (eye(size(G)) + D) \ (eye(size(G)) + G);

    epsilon = 1/max(d);

    W_B = eye(size(G)) - epsilon * L;
end