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


%G_ring = make_ring();
%G_ring_stat = graph_stats(G_ring);