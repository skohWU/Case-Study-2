n = 100;

rng(105); 
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];

G_disk = make_disk(P0,3);
[degree_disk, edge_count_disk, diameter_disk] = graph_stats(G_disk);


%G_ring = make_ring();
%G_ring_stat = graph_stats(G_ring);