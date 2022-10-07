function [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub,F,cap_non_rem]=capacities(G,K,...
    E,N,set_non_rem)
size_set_non_rem=size(set_non_rem,1);
cap=zeros(E,K);
cap_non_rem = zeros(E-size_set_non_rem,K);
k=2;
F = rand(E-size_set_non_rem,k);
xi_mu = rand(k,1);
for it=1:K
    cap_non_rem(:,it) = F*expinv(rand(k,1),xi_mu);
end
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
