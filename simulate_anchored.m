function X = simulate_anchored(W, lamda, a, x0, K)
    Lamda = diag(lamda);

    X = zeros(size(W,1),K+1);
    X(:,1) = x0;

    for i = 1:K
        X(:,i+1) = Lamda * W * X(:,i) + (eye(size(W)) - Lamda) * a;
    end

    round(norm(X(:,K+1) - W^K * x0));
end