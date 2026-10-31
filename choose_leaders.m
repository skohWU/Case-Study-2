function S = choose_leaders(G, P0, faulty, g, m)
%CHOOSE_LEADERS  Choose which robots become leaders (ESE 105, Case Study II, Part III).
%
%   S = choose_leaders(G, P0, faulty, g, m)
%
%   Inputs
%     G       n-by-n adjacency matrix of the communication network
%     P0      n-by-2 landing positions; row i is robot i's (x, y)
%     faulty  vector of indices of robots with broken drives
%     g       1-by-2 rally point
%     m       maximum number of leaders
%
%   Output
%     S       vector of at most m distinct robot indices, none of them faulty
%
%   Rules (see Section D.4 of the handout)
%     - Use only the network's structure: anything you can compute from G,
%       P0, faulty and g with the tools from Part I (degrees, hop counts,
%       distances, products and powers of G, ...).
%     - Do NOT simulate the swarm, and do NOT call score_mission.
%     - Must run in under 2 seconds on one deployment.
%     - Do not change the function name or the order of the inputs.
%       Helper functions are allowed, but put them in this file, below
%       this function, so that this single file is all we need to run it.
%
<<<<<<< HEAD
%   Team: <Plswork>
=======
%   Team: Group 12
>>>>>>> main
%   Members: Spencer Koh, Yijia Mao, Dayanis Morales Salgado
%
%   Strategy (2-3 sentences describing your rule, in words):
%   <describe it here>
%   Prelim strategy 1: 
%   Select the most strongly connected robots as leaders, 
% given the connections aren't to a faulty one.
n = size(G, 1);
working = true(n, 1);
working(faulty) = false;
S = zeros(m, 1);
% ---- replace this placeholder with your strategy ----
%Strategy 1: Most connected

[degree, ~, ~] = graph_stats(G);

working_robots = false(size(G,1),1);
working_robots(faulty) = true;
broken_indices = find(working_robots);

degree(broken_indices) = 0;

[~,S] = maxk(degree,m);


end

