function [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub,F]=capacities_dirichlet(G,K,...
    E,N,set_non_rem, cap_non_rem)
size_set_non_rem=size(set_non_rem,1);
cap=zeros(E,K);
h=1;
for e=1:E
    if ismember(e,set_non_rem)
        cap(e,:)=100;      % infinite capacity of non-removable arcs
    else
        cap(e,:)=cap_non_rem(h,:);
        h=h+1;
    end
end

zeta_lb=0; % lowerbound on flow for any interdction plan and scenario
D_max=zeros(E,1);
for e=1:E
    D_max(e)=max(cap(e,:)');
end
G.Edges.Weight=D_max;
zeta_ub=maxflow(G,1,size(N,1));
diag_cap=zeros(E,E,K);
diag_cap_non_rem=zeros(E-size_set_non_rem,E-size_set_non_rem,K);
for j=1:K
    diag_cap(:,:,j)=diag(cap(:,j));
    diag_cap_non_rem(:,:,j)=diag(cap_non_rem(:,j));
end
end
