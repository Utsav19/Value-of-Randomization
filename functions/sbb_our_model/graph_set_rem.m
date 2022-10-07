function [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r)

H = incidence(G);N=full(H);E=size(N,2);
suc_source=successors(G,1);
edges_source=[ones(size(suc_source,1),1) suc_source];
predec_sink=predecessors(G,p*r+2);
edges_sink=[predec_sink (p*r+2)*ones(size(predec_sink,1),1)];

set_non_rem=[];
for i=1:size(suc_source,1)-1
    if findedge(G,suc_source(i),suc_source(i+1))>0
        set_non_rem=[set_non_rem;[suc_source(i) suc_source(i+1)]];
    elseif findedge(G,suc_source(i+1),suc_source(i))>0
        set_non_rem=[set_non_rem;[suc_source(i+1) suc_source(i)]];
    end
end
for i=1:size(edges_sink,1)-1
    if findedge(G,predec_sink(i), predec_sink(i+1))>0
        set_non_rem=[set_non_rem;[predec_sink(i) predec_sink(i+1)]];
    elseif findedge(G,predec_sink(i+1), predec_sink(i))>0
        set_non_rem=[set_non_rem;[predec_sink(i+1) predec_sink(i)]];
    end
end
set_non_rem=[set_non_rem;edges_source;edges_sink];
a1= set_non_rem(:,1);
a2= set_non_rem(:,2);
set_non_rem=sort(findedge(G,a1,a2)');
set_non_rem= set_non_rem';
set_rem=[];
for e=1:E
    if ismember(e,set_non_rem)==0
        set_rem=[set_rem;e];
    end
end
end
