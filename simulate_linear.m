function X = simulate_linear(W, x0, K)
    X = zeros(size(W,1),K+1);
    X(:,1) = x0;

    for i = 1:K
        X(:,i+1) = W * X(:,i);
    end
end