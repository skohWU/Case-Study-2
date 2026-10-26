clearvars;
close all;

n = 100;

rng(105); 
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];

% Finding r_min
r = 0;

while true
    G_disk_rmin = make_disk(P0,r);
    [~, ~, diameter_rmin] = graph_stats(G_disk_rmin);
    r = r + 0.1;

    if diameter_rmin ~= inf
        r_min = r;
        break
    end
end

G_disk = make_disk(P0,r_min);
[degree_disk, edge_count_disk, diameter_disk] = graph_stats(G_disk);



% Finding ring's k.
k = 1;

while true
    G_ring_k0 = make_ring(n,k);
    [~, edge_count_ring_k0, ~] = graph_stats(G_ring_k0);

    G_ring_k1 = make_ring(n,k+1);
    [~, edge_count_ring_k1, ~] = graph_stats(G_ring_k1);

    if abs(edge_count_disk - edge_count_ring_k1) < abs(edge_count_disk - edge_count_ring_k0)
        k = k + 1;
    else
        break
    end
end

% Getting the stats for each network.

G_ring = make_ring(n,k);
[degree_ring, edge_count_ring, diameter_ring] = graph_stats(G_ring);

G_star = make_star(n);
[degree_star, edge_count_star, diameter_star] = graph_stats(G_star);

G_complete = make_complete(n);
[degree_complete, edge_count_complete, diameter_complete] = graph_stats(G_complete);

% Make a table comparing different networks.

stats = table([degree_complete' * ones(n,1) * 1/100, degree_star' * ones(n,1) * 1/100 , degree_ring' * ones(n,1) * 1/100 , degree_disk' * ones(n,1) * 1/100]',[max(degree_complete), max(degree_star), max(degree_ring), max(degree_disk)]', [edge_count_complete, edge_count_star, edge_count_ring, edge_count_disk]', [diameter_complete, diameter_star, diameter_ring, diameter_disk]', VariableNames = ["Mean Degree", "Max Degree", "Edges", "Diameter"], RowNames=["Complete", "Star", "Ring", "Disk"]);

% Graph disk network with connections.

figure('Name',"Disk Graph Network");
plot(graph(G_disk), 'XData', P0(:,1), 'YData', P0(:,2), 'NodeColor','r', 'NodeLabel',{});
title("Disk Graph Network");
grid on;

% Make each protocol, A and B, for each network.

[W_A_complete, W_B_complete] = makeProtocols(degree_complete,G_complete);
[W_A_ring, W_B_ring] = makeProtocols(degree_ring,G_ring);
[W_A_star, W_B_star] = makeProtocols(degree_star,G_star);
[W_A_disk, W_B_disk] = makeProtocols(degree_disk,G_disk);

K = 100;

X = simulate_linear(W_A_disk,P0(:,1),K);
Y = simulate_linear(W_A_disk,P0(:,2),K);

P_final_disk = [X(:,K+1), Y(:,K+1)];

figure('Name',"Graph after Updating Positions");
plot(graph(G_disk), 'XData', P_final_disk(:,1), 'YData', P_final_disk(:,2), 'NodeColor','r', 'NodeLabel',{});
title("Graph after Updating Positions");
grid on;