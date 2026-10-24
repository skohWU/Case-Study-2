function [degree, edge_count, diameter] = graph_stats(G)
    
    n = size(G,1);
    degree = G * ones(n,1);
    
    edge_count = 1/2 * ones(n,1)' * G * ones(n,1);

    I = eye(n,n);
    R = I + G;
    k = 1;

    while true
        if any(R == 0,"all")
            R = (I + G) * R;
            R(R ~= 0) = 1;
            k = k + 1;
        else
            diameter = k;
            break
        end

        if k > 1000
            diameter = inf;
            break
        end
    end

end