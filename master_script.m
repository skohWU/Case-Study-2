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

figure('Name',"Disk Graph Network", 'Theme', "light");
plot(graph(G_disk), 'XData', P0(:,1), 'YData', P0(:,2), 'NodeColor','r', 'NodeLabel',{});
title("Disk Graph Network");
grid on;
exportgraphics(gcf, 'diskGraphNetwork.png', 'BackgroundColor', 'white');

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

delta_complete = zeros(K+1,1);

for k = 1:K+1
    delta_complete(k) = sqrt(norm(C * X_complete(:,k))^2 + norm(C * Y_complete(:,k))^2);
end

k_epsilon_complete = find(delta_complete <= 0.01*delta_complete(1),1);

% Simulate movement and create delta vector for T2 (Star)

X_star = simulate_linear(W_A_star,P0(:,1),K);
Y_star = simulate_linear(W_A_star,P0(:,2),K);

delta_star = zeros(K+1,1);

for k = 1:K+1
    delta_star(k) = sqrt(norm(C * X_star(:,k))^2 + norm(C * Y_star(:,k))^2);
end

k_epsilon_star = find(delta_star <= 0.01*delta_star(1),1);

% Simulate movement and create delta vector for T3 (Ring)

X_ring = simulate_linear(W_A_ring,P0(:,1),K);
Y_ring = simulate_linear(W_A_ring,P0(:,2),K);

delta_ring = zeros(K+1,1);

for k = 1:K+1
    delta_ring(k) = sqrt(norm(C * X_ring(:,k))^2 + norm(C * Y_ring(:,k))^2);
end

k_epsilon_ring = find(delta_ring <= 0.01*delta_ring(1),1);

% Simulate movement and create delta vector for T4 (Disk)

X_disk = simulate_linear(W_A_disk,P0(:,1),K);
Y_disk = simulate_linear(W_A_disk,P0(:,2),K);

delta_disk = zeros(K+1,1);

for k = 1:K+1
    delta_disk(k) = sqrt(norm(C * X_disk(:,k))^2 + norm(C * Y_disk(:,k))^2);
end

k_epsilon_disk = find(delta_disk <= 0.01*delta_disk(1),1);

% Simulate movement and create delta vector for T4 with r + .4

disk_plus = make_disk(P0, r_min + 0.4);

[degree_plus, ~, ~] = graph_stats(disk_plus);

[W_A_plus, ~] = makeProtocols(degree_plus,disk_plus);

X_plus = simulate_linear(W_A_plus,P0(:,1),K);
Y_plus = simulate_linear(W_A_plus,P0(:,2),K);

delta_plus = zeros(K+1,1);

for k = 1:K+1
    delta_plus(k) = sqrt(norm(C * X_plus(:,k))^2 + norm(C * Y_plus(:,k))^2);
end

k_epsilon_plus = find(delta_plus <= 0.01*delta_plus(1),1);

figure(Theme="light");
hold on;
plot(delta_complete);
plot(delta_star);
plot(delta_ring);
plot(delta_disk);
plot(delta_plus)
legend("Complete", "Star", "Ring", "Disk", "Disk, r = r + 0.4");
yscale("log");
hold off;

exportgraphics(gcf, 'deltas.png');

% figure('Name',"Graph after Updating Positions, Disk");
% plot(graph(G_disk), 'XData', X_disk(:,K+1), 'YData', Y_disk(:,K+1), 'NodeColor','r', 'NodeLabel',{});
% title("Graph after Updating Positions, Disk");
% grid on;

% Simulate a Deployment with Leaders
seeds = randi([1, 250], 1, 30);

num_leaders = 4;

test = make_deployment(seeds(1));
leaders = choose_leaders(test.G, test.P0, test.faulty, test.g, num_leaders);

n = size(test.G,1);

lamda = ones(n,1);
%Put in code to iterate through lambda leaders
for i = leaders
    lamda(i) = 0.5;
end
lamda(test.faulty) = 0;

a = test.P0(:,1);
a(leaders) = test.g(1);
b = test.P0(:,2);
b(leaders) = test.g(2);

figure('Name',"Test Network Start");
plot(graph(test.G), 'XData', test.P0(:,1), ...
    'YData', test.P0(:,2), 'NodeColor','r', 'NodeLabel', ...
    {});
hold on;
plot(test.g(1), test.g(2), 'gp', 'MarkerSize',15, 'Color','g');
title("Test Network Start");
grid on;

[degree, ~, ~] = graph_stats(test.G);

[W, ~] = makeProtocols(degree, test.G);

X = simulate_anchored(W, lamda, a, test.P0(:,1), n);
Y = simulate_anchored(W, lamda, b, test.P0(:,2), n);

figure('Name',"Test Network End");
plot(graph(test.G), 'XData', X(:,n+1), ...
    'YData', Y(:,n+1), 'NodeColor','r', 'NodeLabel', ...
    {});
hold on;
plot(test.g(1), test.g(2), 'gp', 'MarkerSize',15, 'Color','g');
title("Test Network End");
grid on;

test_score = score_mission(test, leaders);

% Simple leader selection functions

function S = random_working(G, faulty, m)
    faulty_robots = true(size(G,1),1);
    faulty_robots(faulty) = false;
    working_indices = find(faulty_robots);

    S = working_indices(randperm(size(working_indices,1), m));
end


function S = highest_degree(G, faulty, m)
    [degree, ~, ~] = graph_stats(G);

    faulty_robots = false(size(G,1),1);
    faulty_robots(faulty) = true;
    broken_indices = find(faulty_robots);

    degree(broken_indices) = 0;

    [~,S] = maxk(degree,m);
end

function S = closest_to_goal(G, P0,faulty, g, m)
    faulty_robots = false(size(G,1),1);
    faulty_robots(faulty) = true;
    broken_indices = find(faulty_robots);

    distance = vecnorm(P0 - g, 2, 2);
    distance(broken_indices) = 10000;
    [~,S] = mink(distance,m);

end

function S = distances_to_goal(G, P0,faulty, g, m)
%Modification: Returns a n-m x m matrix containing
%different S 
n = size(G, 1);
S = zeros(n-m-length(faulty)+1, m);
faulty_robots = false(size(G,1),1);
faulty_robots(faulty) = true;
broken_indices = find(faulty_robots);

distance = vecnorm(P0 - g, 2, 2);
distance(broken_indices) = inf;

[~,newArr] = sort(distance);

for i = 1:size(S,1)
    S(i,:) = newArr(i:i+m-1);
end

end


% Comparing Simple Leader Selections

random_scores = zeros(size(seeds));
highest_scores = zeros(size(seeds));
closest_scores = zeros(size(seeds));

for i = 1:size(seeds,2)
    D = make_deployment(seeds(i));

    random_leaders = random_working(D.G,D.faulty,num_leaders);
    highest_degree_leaders = highest_degree(D.G, D.faulty, num_leaders);
    closest_leaders = closest_to_goal(D.G, D.P0, D.faulty, D.g, num_leaders);

    random_scores(i) = score_mission(D, random_leaders);
    highest_scores(i) = score_mission(D, highest_degree_leaders);
    closest_scores(i) = score_mission(D, closest_leaders);
end

figure('Theme', "light");
hold on;
plot(random_scores);
plot(highest_scores);
plot(closest_scores);
legend("Random Leaders", "Highest Degree Leaders", "Closest to g Leaders");
ylabel("Score");
xlabel("Simulation")
hold off;

exportgraphics(gcf, 'simpleSelection.png');

writematrix([random_scores', highest_scores', closest_scores'],"Simple_selection_scores.csv");

% Looking at one deployment, seed 19

seed19 = make_deployment(seeds(19));

figure('Name',"Seed 19", 'Theme', "light");
plot(graph(seed19.G), 'XData', seed19.P0(:,1), ...
    'YData', seed19.P0(:,2), 'NodeColor','r', 'NodeLabel', ...
    {});
hold on;
plot(seed19.g(1), seed19.g(2), '*', 'MarkerSize',15, 'Color',[1,0,1]);
title("Seed 19 Start");
grid on;

exportgraphics(gcf, 'seed19Start.png');

n = size(seed19.G,1);

[degree, ~, ~] = graph_stats(seed19.G);

[W, ~] = makeProtocols(degree, seed19.G);

% Random Leaders on 19

leaders = random_working(seed19.G, seed19.faulty, num_leaders);
<<<<<<< HEAD

lamda = ones(n,1);
for i = leaders
    lamda(i) = 0.5;
end
lamda(seed19.faulty) = 0;

a = seed19.P0(:,1);
a(leaders) = seed19.g(1);
b = seed19.P0(:,2);
b(leaders) = seed19.g(2);

X = simulate_anchored(W, lamda, a, seed19.P0(:,1), n);
Y = simulate_anchored(W, lamda, b, seed19.P0(:,2), n);

figure('Name',"Random Leaders 19 End");
plot(graph(seed19.G), 'XData', X(:,n+1), ...
    'YData', Y(:,n+1), 'NodeColor','r', 'NodeLabel', ...
    {});
hold on;
plot(seed19.g(1), seed19.g(2), 'gp', 'MarkerSize',15, 'Color','g');
title("Random Leaders 19 End");
grid on;
=======
export_for_viewer(seed19,leaders,"seed19random.json");
>>>>>>> main

% Highest degree Leaders on 19

leaders = highest_degree(seed19.G, seed19.faulty, num_leaders);
<<<<<<< HEAD

lamda = ones(n,1);
for i = leaders
    lamda(i) = 0.5;
end
lamda(seed19.faulty) = 0;

a = seed19.P0(:,1);
a(leaders) = seed19.g(1);
b = seed19.P0(:,2);
b(leaders) = seed19.g(2);

X = simulate_anchored(W, lamda, a, seed19.P0(:,1), n);
Y = simulate_anchored(W, lamda, b, seed19.P0(:,2), n);

figure('Name',"Highest Deg Leaders 19 End");
plot(graph(seed19.G), 'XData', X(:,n+1), ...
    'YData', Y(:,n+1), 'NodeColor','r', 'NodeLabel', ...
    {});
hold on;
plot(seed19.g(1), seed19.g(2), 'gp', 'MarkerSize',15, 'Color','g');
title("Highest Deg Leaders 19 End");
grid on;
=======
export_for_viewer(seed19,leaders,"seed19Highest.json");
>>>>>>> main

% Closest to g Leaders on 19

leaders = closest_to_goal(seed19.G, seed19.P0, seed19.faulty, seed19.g, num_leaders);
<<<<<<< HEAD

lamda = ones(n,1);
for i = leaders
    lamda(i) = 0.5;
end
lamda(seed19.faulty) = 0;

a = seed19.P0(:,1);
a(leaders) = seed19.g(1);
b = seed19.P0(:,2);
b(leaders) = seed19.g(2);

X = simulate_anchored(W, lamda, a, seed19.P0(:,1), n);
Y = simulate_anchored(W, lamda, b, seed19.P0(:,2), n);

figure('Name',"Closest Leaders 19 End");
plot(graph(seed19.G), 'XData', X(:,n+1), ...
    'YData', Y(:,n+1), 'NodeColor','r', 'NodeLabel', ...
    {});
hold on;
plot(seed19.g(1), seed19.g(2), 'gp', 'MarkerSize',15, 'Color','g');
title("Closest Leaders 19 End");
grid on;

% Distance to g Leaders on 19

leaders = distances_to_goal(seed19.G, seed19.P0, seed19.faulty, seed19.g, num_leaders);
count = 0;
for j = leaders
    lamda = ones(n,1);
    for i = j
        lamda(i) = 0.5;
    end
    lamda(seed19.faulty) = 0;
    
    a = seed19.P0(:,1);
    a(j) = seed19.g(1);
    b = seed19.P0(:,2);
    b(j) = seed19.g(2);
    
    X = simulate_anchored(W, lamda, a, seed19.P0(:,1), n);
    Y = simulate_anchored(W, lamda, b, seed19.P0(:,2), n);
    
    figure('Name', "#" + count + " Leaders 19 End");
    p = plot(graph(seed19.G), 'XData', X(:,n+1), ...
        'YData', Y(:,n+1), 'NodeColor','r', 'NodeLabel', ...
        {});
    hold on;
    highlight(p, j, 'NodeColor', 'b', 'MarkerSize', 3);
    plot(seed19.g(1), seed19.g(2), 'gp', 'MarkerSize',15, 'Color','g');
    title("Closest Leaders 19 End");
    grid on;
    count = count + 1;
end
=======
export_for_viewer(seed19,leaders,"seed19Closest.json");
>>>>>>> main
