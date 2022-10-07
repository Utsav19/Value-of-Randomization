

function fl_det = flows_test_det(G, deter_p, cap, N)
    G.Edges.Weight=(1-deter_p).*cap(:,k);
    fl_det=maxflow(G,1,size(N,1));