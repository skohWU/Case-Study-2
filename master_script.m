clearvars;
close all;



n = 100;

rng(105); 
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];

%% Part I: Communication Networks

% Finding r_min: smallest radius (to the nearest 0.1) giving a connected
% disk graph. graph_stats returns diameter = inf when the network is
% disconnected 
r = 0;

while true
    G_disk_rmin = make_disk(P0,r);
    [~, ~, diameter_rmin] = graph_stats(G_disk_rmin);
    if diameter_rmin ~= inf
        r_min = r;
        break                           % r is the first connected radius 
    end
    r = r + 0.1;                        
end

G_disk = make_disk(P0,r_min);
[degree_disk, edge_count_disk, diameter_disk] = graph_stats(G_disk);



% Finding ring's k: increase k while k+1 gets the ring's edge count closer
% to the disk graphs edge count
k = 1;

while true
    G_ring_k0 = make_ring(n,k);
    [~, edge_count_ring_k0, ~] = graph_stats(G_ring_k0);

    G_ring_k1 = make_ring(n,k+1);
    [~, edge_count_ring_k1, ~] = graph_stats(G_ring_k1);

    if abs(edge_count_disk - edge_count_ring_k1) < ...
            abs(edge_count_disk - edge_count_ring_k0)
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
% Columns: mean degree, max degree, edge count, diameter
% (degree' * ones(n,1) * 1/100 is the sum of degrees divided by n = 100,
% i.e. the mean degree)

stats = table([degree_complete' * ones(n,1) * 1/100, degree_star' * ...
    ones(n,1) * 1/100 , degree_ring' * ones(n,1) * 1/100 , degree_disk' * ...
    ones(n,1) * 1/100]',[max(degree_complete), max(degree_star), ...
    max(degree_ring), max(degree_disk)]', [edge_count_complete, ...
    edge_count_star, edge_count_ring, edge_count_disk]', ...
    [diameter_complete, diameter_star, diameter_ring, diameter_disk]', ...
    VariableNames = ["Mean Degree", "Max Degree", "Edges", "Diameter"], ...
    RowNames=["Complete", "Star", "Ring", "Disk"]);

% Graph disk network with connections.

figure('Name',"Disk Graph Network", 'Theme', "light");
plot(graph(G_disk), 'XData', P0(:,1), 'YData', P0(:,2), 'NodeColor','r', 'NodeLabel',{});
title("Disk Graph Network");
grid on;
exportgraphics(gcf, 'diskGraphNetwork.png', 'BackgroundColor', 'white');


%% Part II: Rendezvous Protocols 

% Make each protocol, A and B, for each network.

[W_A_complete, W_B_complete] = makeProtocols(degree_complete,G_complete);
[W_A_ring, W_B_ring] = makeProtocols(degree_ring,G_ring);
[W_A_star, W_B_star] = makeProtocols(degree_star,G_star);
[W_A_disk, W_B_disk] = makeProtocols(degree_disk,G_disk);

K = 100;                                % number of simulation steps 

% Centering matrix --> C*x subtracts the average from x
C = eye(n) - 1/n * ones(n,1) * ones(n,1)';



% The next five blocks all do the same thing for one network:
% 1) simulate x and y with Protocol A (columns are x[0]...x[K])
% 2) compute delta[k], the spread of the swarm, for every step
% 3) find k_epsilon, the first step with delta <= 1% of delta[0]

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

% Plot delta[k] for all networks on one figure (log scale)

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



% k_epsilon values only
K_e1 = 5000;    % long enough for slow networks to converge
Ws    = {W_A_complete, W_A_star, W_A_ring, W_A_disk, W_A_plus};
names = ["Complete","Star","Ring","Disk (r_min)","Disk (r_min+0.4)"];

for i = 1:5
    X = simulate_linear(Ws{i}, P0(:,1), K_e1);
    Y = simulate_linear(Ws{i}, P0(:,2), K_e1);
    d = sqrt(sum((C*X).^2,1) + sum((C*Y).^2,1));   % delta[0] ... delta[K_e1]
    idx = find(d <= 0.01*d(1), 1);
    if isempty(idx)
        fprintf('%-18s k_eps not reached in %d steps\n', names(i), K_e1);
    else
        % -1 since column 1 is x[0]
        fprintf('%-18s k_eps = %d\n', names(i), idx - 1);   
    end
end


[degree_plus, edge_count_plus, diameter_plus] = graph_stats(disk_plus);
fprintf('r_min = %.1f\n', r_min);
fprintf('Disk r_min+0.4: edges = %d, mean deg = %.2f, max deg = %d, diameter = %d\n', ...
    edge_count_plus, mean(degree_plus), max(degree_plus), diameter_plus);
fprintf('cross-zone links: r_min = %d, r_min+0.4 = %d\n', ...
    sum(sum(G_disk(1:50,51:100))), sum(sum(disk_plus(1:50,51:100))));





%% Part E2 Validation

one   = ones(n,1);
names = ["Complete", "Star", "Ring", "Disk"];
Gs    = {G_complete, G_star, G_ring, G_disk};
ds    = {degree_complete, degree_star, degree_ring, degree_disk};
Kbig  = 2^20;

for i = 1:4
    d = ds{i};
    [~, E, ~]  = graph_stats(Gs{i});
    [W_A, W_B] = makeProtocols(d, Gs{i});

% Predicted meeting points 
    v  = one + d;
    cs = {[v'*P0(:,1), v'*P0(:,2)] / (n + 2*E), ... 
        [one'*P0(:,1), one'*P0(:,2)] / n};           

    endA = [mean(W_A^Kbig * P0(:,1)), mean(W_A^Kbig * P0(:,2))];
    endB = [mean(W_B^Kbig * P0(:,1)), mean(W_B^Kbig * P0(:,2))];

    fprintf('%-8s  error A = %.2e   error B = %.2e\n', names(i), norm(cs{1}-endA), norm(cs{2}-endB));
    fprintf('predicted A [%.4f %.4f]  simulated A [%.4f %.4f]\n', cs{1}, endA);
    fprintf('predicted B [%.4f %.4f]  simulated B [%.4f %.4f]\n', cs{2}, endB);
end



%% Part III: Simulate a deployment with leaders

% Pick 30 random practice seeds 
seeds = randi([1, 250], 1, 30);

num_leaders = 4;

% Build one test deployment and choose leaders with our strategy 
test = make_deployment(seeds(1));
leaders = choose_leaders(test.G, test.P0, test.faulty, test.g, num_leaders);

n = size(test.G,1);

% Responsiveness: followers 1, leaders 0.5, faulty 0
lamda = ones(n,1);
% Put in code to iterate through lambda leaders
for i = leaders
    lamda(i) = 0.5;
end
lamda(test.faulty) = 0;


% Anchors: landing positions, except leaders are anchored at g
a = test.P0(:,1);
a(leaders) = test.g(1);
b = test.P0(:,2);
b(leaders) = test.g(2);

% Plot the network at the start 
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


% Run the anchored simulation for n steps (n = number of robots here)
X = simulate_anchored(W, lamda, a, test.P0(:,1), n);
Y = simulate_anchored(W, lamda, b, test.P0(:,2), n);

% Plot the network at the end
figure('Name',"Test Network End");
plot(graph(test.G), 'XData', X(:,n+1), ...
    'YData', Y(:,n+1), 'NodeColor','r', 'NodeLabel', ...
    {});
hold on;
plot(test.g(1), test.g(2), 'gp', 'MarkerSize',15, 'Color','g');
title("Test Network End");
grid on;

test_score = score_mission(test, leaders);



% Steady state (no leaders)
% Solve (I - Lam*W) x* = (I - Lam) a with backslash, then compare to a 
% long simulation. err should be near 0. mean_dist_to_g shows how far 
% the swarm ends up from the rally point.


%% D.1: Steady state (no leaders)

D = make_deployment(seeds(1));
n = size(D.G,1);

[deg, ~, ~] = graph_stats(D.G);
[W, ~]      = makeProtocols(deg, D.G);

lam = ones(n,1);   lam(D.faulty) = 0;    % no leaders, faulty robots stuck
a = D.P0(:,1);     b = D.P0(:,2);        % anchors = landing positions
Lam = diag(lam);

% Steady state by backslash
M     = eye(n) - Lam*W;
xstar = M \ ((eye(n) - Lam) * a);
ystar = M \ ((eye(n) - Lam) * b);

res = norm(xstar - (Lam*W*xstar + (eye(n) - Lam)*a))

rho = max(abs(eig(Lam*W)))                    % must be less than 1
K   = min(ceil(log(1e-10) / log(rho)), 50000)

X = simulate_anchored(W, lam, a, D.P0(:,1), K);
Y = simulate_anchored(W, lam, b, D.P0(:,2), K);

% VALUE 1: solved steady state vs. long-run simulation
err = max(norm(xstar - X(:,end)), norm(ystar - Y(:,end)))

% VALUE 2: mean distance of the working robots from the rally point
work = setdiff(1:n, D.faulty);
mean_dist_to_g = mean(vecnorm([xstar(work) ystar(work)] - D.g(:)', 2, 2))

%% Simple leader selection functions


% Strategy 1: m random working (non-faulty) robots 
function S = random_working(G, faulty, m)
    faulty_robots = true(size(G,1),1);
    faulty_robots(faulty) = false;
    working_indices = find(faulty_robots);

    S = working_indices(randperm(size(working_indices,1), m));
end

% Strategy 2: the m working robots with the highest degree 
% (faulty robots get degree 0 so they are never picked)
function S = highest_degree(G, faulty, m)
    [degree, ~, ~] = graph_stats(G);

    faulty_robots = false(size(G,1),1);
    faulty_robots(faulty) = true;
    broken_indices = find(faulty_robots);

    degree(broken_indices) = 0;

    [~,S] = maxk(degree,m);
end


% Strategy 3: the m working robots that landed closest to g 
% (faulty robots get a huge distance so they are never picked)
function S = closest_to_goal(G, P0,faulty, g, m)
    faulty_robots = false(size(G,1),1);
    faulty_robots(faulty) = true;
    broken_indices = find(faulty_robots);

    distance = vecnorm(P0 - g, 2, 2);
    distance(broken_indices) = 10000;
    [~,S] = mink(distance,m);

end

% Exploration helper: returns every window of m consecutive robots in order
% of distance to g (not used in the comparison below)
function S = distances_to_goal(G, P0,faulty, g, m)
% Modification: Returns a n-m x m matrix containing
% different S
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

%% Comparing Simple Leader Selections
% Score each simple strategy on every seed

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

% Plot the scores of the three strategies on all 30 deployments 
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
export_for_viewer(seed19,leaders,"seed19random.json");

% Highest degree Leaders on 19

leaders = highest_degree(seed19.G, seed19.faulty, num_leaders);
export_for_viewer(seed19,leaders,"seed19Highest.json");

% Closest to g Leaders on 19

leaders = closest_to_goal(seed19.G, seed19.P0, seed19.faulty, seed19.g, num_leaders);
export_for_viewer(seed19,leaders,"seed19Closest.json");

grid on;



%% Testing our leader selection
% Goal: see how well our choose_leaders strategy does across many different
% practice deployments, not just one

% One slot per seed to store that deployment's score J. (size(seeds(), 2)
% is just the number of seeds, 30 here)
scores = zeros(size(seeds(),2),1);

% Go through every practice deployment, one seed at a time 
for i = 1:size(seeds(),2)

    % Build deployment i. The same seed always gives the same deployment
    % (network G, landing positions PO, faulty robots, and rally points g)
    D = make_deployment(seeds(i));

    % It only looks at the network and positions, and returns at most
    % num_leaders(4) working robots 
    leaders = choose_leaders(D.G, D.P0, D.faulty, D.g, num_leaders);

    % Run the mission with those leaders and record the score J for this
    % deployment (the average percent reduction in distance to g, compared
    % with having no leaders) 
    scores(i) = score_mission(D,leaders);
end

% Average score across all deployments. This is the number the competition
% uses (average over the hidden deployments) 
average_score = sum(scores) / length(scores);



%% Test Different Distance Based Leaders Across our 30 Seeds

% distanceScores = zeros(size(leaders,1),size(seeds,2));
% 
% for i = 1:size(seeds,2)
%     D = make_deployment(seeds(i));  
%     leaders = distances_to_goal(D.G, D.P0, D.faulty, D.g, num_leaders);
%     % count = 0;
% 
% 
%     for j = 1:size(leaders,1)
% 
%         distanceScores(j,i) = score_mission(D,leaders(j,:));
% 
%         % lamda = ones(n,1);
%         % for i = j
%         %     lamda(i) = 0.5;
%         % end
%         % lamda(seed19.faulty) = 0;
%         % 
%         % 
%         % 
%         % a = seed19.P0(:,1);
%         % a(j) = seed19.g(1);
%         % b = seed19.P0(:,2);
%         % b(j) = seed19.g(2);
%         % 
%         % X = simulate_anchored(W, lamda, a, seed19.P0(:,1), n);
%         % Y = simulate_anchored(W, lamda, b, seed19.P0(:,2), n);
%         % 
%         % figure('Name', "#" + count + " Leaders 19 End");
%         % p = plot(graph(seed19.G), 'XData', X(:,n+1), ...
%         %     'YData', Y(:,n+1), 'NodeColor','r', 'NodeLabel', ...
%         %     {});
%         % hold on;
%         % highlight(p, j, 'NodeColor', 'b', 'MarkerSize', 3);
%         % plot(seed19.g(1), seed19.g(2), 'gp', 'MarkerSize',15, 'Color','g');
%         % title("Closest Leaders 19 End");
%         % grid on;
% 
%         % count = count + 1;
%     end
% end

% writematrix(distanceScores,"distanceScores.csv");

%% Find weights for choose_leaders
% 
% step_size = 0.05; 
% max_score = -Inf;
% best_weights = [];
% 
% 
% for w1 = 0:step_size:1
%     for w2 = 0:step_size:(1 - w1) 
%         w3 = 1 - w1 - w2;
% 
%         current_weights = [w1, w2, w3];
% 
%         scores = zeros(size(seeds(),2),1);
% 
%         for i = 1:size(seeds(),2)
%             D = make_deployment(seeds(i));
%             scores(i) = score_mission(D,choose_leaders(D.G, D.P0, D.faulty, D.g, num_leaders, w1, w2, w3));
%         end
% 
%         current_score = sum(scores) / length(scores);
% 
%         % Check if this combination beats the previous best score
%         if current_score > max_score
%             max_score = current_score;
%             best_weights = current_weights;
%         end
% 
%     end
% end
% 
% % 4. Display the results
% fprintf('--- Optimization Results ---\n');
% fprintf('Maximum Score achieved: %.2f\n', max_score);
% fprintf('Optimal Weight 1: %.2f%%\n', best_weights(1) * 100);
% fprintf('Optimal Weight 2: %.2f%%\n', best_weights(2) * 100);
% fprintf('Optimal Weight 3: %.2f%%\n', best_weights(3) * 100);
