clearvars;
close all;

n = 100;

rng(105); 
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];

% Finding r_min
r = 0;

while true
    G_disk_rmin = make_disk(P0,r);
    [degree_rmin, edge_count_rmin, diameter_rmin] = graph_stats(G_disk_rmin);
    r = r + 0.1;

    if diameter_rmin ~= inf
        r_min = r;
        break
    end
end

G_disk = make_disk(P0,r_min);
[degree_disk, edge_count_disk, diameter_disk] = graph_stats(G_disk);


figure('Name',"Disk Graph Network");
plot(graph(G_disk), 'XData', P0(:,1), 'YData', P0(:,2), 'NodeColor','r', 'NodeLabel',{});
title("Disk Graph Network");
grid on;


%G_ring = make_ring();
%G_ring_stat = graph_stats(G_ring);