function fl_ran =  flows_test(G, supp_r, cap, N, K)
    for i=1:K
        G.Edges.Weight=(1-supp_r).*cap(:,i);
        fl_ran(i)=maxflow(G,1,size(N,1));
    end
end