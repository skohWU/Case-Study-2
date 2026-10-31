function [J] = score_mission(D, S)
    T = 125;


    n = size(D.G,1);
    
    % Leaders
    lambda = ones(n,1);
    lambda(S) = 0.5;
    lambda(D.faulty) = 0;

    a = D.P0(:,1);
    a(S) = D.g(1);

    b = D.P0(:,2);
    b(S) = D.g(2);

    [degree, ~, ~] = graph_stats(D.G);
    [W, ~] = makeProtocols(degree, D.G);

    X_new = simulate_anchored(W, lambda, a, D.P0(:,1), T);
    Y_new = simulate_anchored(W, lambda, b, D.P0(:,2), T);

    % Baseline
    lambda_base = ones(n,1);
    lambda_base(D.faulty) = 0;

    X_base = simulate_anchored(W, lambda_base, ...
                               D.P0(:,1), D.P0(:,1), T);

    Y_base = simulate_anchored(W, lambda_base, ...
                               D.P0(:,2), D.P0(:,2), T);

    working_robots = true(n,1);
    working_robots(D.faulty) = false;


    e_new = zeros(T,1);
    e_base = zeros(T,1);
    total = 0;

    for k = 1:T
        e_new(k) = mean( ...
            (X_new(working_robots,k) - D.g(1)).^2 + ...
            (Y_new(working_robots,k) - D.g(2)).^2);

        e_base(k) = mean( ...
            (X_base(working_robots,k) - D.g(1)).^2 + ...
            (Y_base(working_robots,k) - D.g(2)).^2);

        total = total + (1 - e_new(k)/e_base(k));
    end
    
    J = (100/T) * total;

end
