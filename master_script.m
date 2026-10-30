clearvars;
close all;

n = 100;

% rng(105); 
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

C = eye(n) - 1/n * ones(n,1) * ones(n,1)';

% Simulate movement and create delta vector for T1 (Complete)

X_complete = simulate_linear(W_A_complete,P0(:,1),K);
Y_complete = simulate_linear(W_A_complete,P0(:,2),K);

delta_complete = zeros(n,K+1);

for k = 1:K+1
    delta_complete(k) = sqrt(norm(C * X_complete(:,k))^2 + norm(C * Y_complete(:,k))^2);
end

k_epsilon_complete = find(delta_complete <= 0.01*delta_complete(1),1);

% Simulate movement and create delta vector for T2 (Star)

X_star = simulate_linear(W_A_star,P0(:,1),K);
Y_star = simulate_linear(W_A_star,P0(:,2),K);

delta_star = zeros(n,K+1);

for k = 1:K+1
    delta_star(k) = sqrt(norm(C * X_star(:,k))^2 + norm(C * Y_star(:,k))^2);
end

k_epsilon_star = find(delta_star <= 0.01*delta_star(1),1);

% Simulate movement and create delta vector for T3 (Ring)

X_ring = simulate_linear(W_A_ring,P0(:,1),K);
Y_ring = simulate_linear(W_A_ring,P0(:,2),K);

delta_ring = zeros(n,K+1);

for k = 1:K+1
    delta_ring(k) = sqrt(norm(C * X_ring(:,k))^2 + norm(C * Y_ring(:,k))^2);
end

k_epsilon_ring = find(delta_ring <= 0.01*delta_ring(1),1);

% Simulate movement and create delta vector for T4 (Disk)

X_disk = simulate_linear(W_A_disk,P0(:,1),K);
Y_disk = simulate_linear(W_A_disk,P0(:,2),K);

delta_disk = zeros(n,K+1);

for k = 1:K+1
    delta_disk(k) = sqrt(norm(C * X_disk(:,k))^2 + norm(C * Y_disk(:,k))^2);
end

k_epsilon_disk = find(delta_disk <= 0.01*delta_disk(1),1);

figure('Name',"Graph after Updating Positions, Disk");
plot(graph(G_disk), 'XData', X_disk(:,K+1), 'YData', Y_disk(:,K+1), 'NodeColor','r', 'NodeLabel',{});
title("Graph after Updating Positions, Disk");
grid on;

% Leaders and Faulty Robots

lamda = ones(100,1);
a = P0(:,1);
b = P0(:,2);

% Simulate a Deployment with Leaders

num_leaders = 4;

test = make_deployment(rand*100);
leaders = choose_leaders(test.G, test.P0, test.faulty, test.g, num_leaders);

n = size(test.G,1);

lamda = ones(n,1);
lamda(leaders) = .5;
lamda(test.faulty) = 0;

a = test.P0(:,1);
a(leaders) = test.g(1);
b = test.P0(:,2);
a(leaders) = test.g(2);

figure('Name',"Test Network Start");
plot(graph(test.G), 'XData', test.P0(:,1), 'YData', test.P0(:,2), 'NodeColor','r', 'NodeLabel',{});
title("Test Network Start");
grid on;

[degree, ~, ~] = graph_stats(test.G);

[W, ~] = makeProtocols(degree, test.G);

X = simulate_anchored(W, lamda, a, test.P0(:,1), n);
Y = simulate_anchored(W, lamda, b, test.P0(:,2), n);

% Simple leader selection functions

function S = random_working(G, faulty, m)
    working_robots = true(size(G,1),1);
    working_robots(faulty) = false;
    working_indices = find(working_robots);

    S = working_indices(randperm(size(working_indices,1), m));
end

random_leaders = random_working(test.G,test.faulty,num_leaders);
J = score_mission(test, random_leaders);