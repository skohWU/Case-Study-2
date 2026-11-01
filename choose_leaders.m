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
%   Team: Group 12
%   Members: Spencer Koh, Yijia Mao, Dayanis Morales Salgado
%
%   Strategy (2-3 sentences describing your rule, in words):
%   <describe it here>
%   Prelim strategy 1: 
%   Select the robots that have the highest weighted score. Score is
%   weighted between number of connections, distance to the goal, and
%   separation from other leaders.

working_robots = true(size(G,1), 1);
working_robots(faulty) = false;
working_indices = find(working_robots);

S = [];

degree = sum(G, 2);
degree = degree / (max(degree));

distance_to_g = sqrt(sum((P0 - g).^2, 2));
distance_score = distance_to_g / (max(distance_to_g));

for k = 1:m

    best_score = -inf;
    best_robot = -1;

    for i = working_indices'
        if ismember(i, S)
            continue;
        end
        if isempty(S)
            separation = 0;
        else
            distances = sqrt(sum((P0(S,:) - P0(i,:)).^2, 2));

            separation = min(distances);
            separation = separation/(max(distance_to_g));
            separation = min(separation, 1);
        end

        score = .5 * degree(i) + .2 * distance_score(i) + .3 * separation;

        if score > best_score
            best_score = score;
            best_robot = i;
        end
    end

    S(end+1) = best_robot;

end

end
