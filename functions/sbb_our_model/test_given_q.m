function [cvar_opt,zeta]=test_given_q(fl,K,alpha,q,u)
L=size(u,1);
zeta = sdpvar(1,1);
Del = sdpvar(L,K,'full');
objective = zeta+(1/(1-alpha))*sum(q*Del');
constraints = [Del(:)>=0];
for l=1:L
		constraints= constraints+(Del(l,:)>=u(l,1)*fl(l,:)-u(l,1)*zeta);
end
options=sdpsettings('solver','CPLEX','verbose',0);
optimize(constraints,objective, options);
cvar_opt = value(objective);
zeta=value(zeta);

% function [cvar_opt,zeta]=test_given_q(fl,K,alpha,q,u,C,supp,N,E)
% L=size(u,1);
% zeta = sdpvar(1,1);
% Del = sdpvar(L,K,'full');
% gam =sdpvar(L,size(N,1)-2,K,'full');%dual
% lam=sdpvar(L,E,K,'full');    %dual
% objective = zeta+(1/(1-alpha))*sum(q*Del');
% constraints = [Del(:)>=0,lam(:)>=0,lam(:)<=1];
% for l=1:L
% 	for k=1:K
% 		%constraints= constraints+(Del(l,:)>=u(l,1)*fl(l,:)-u(l,1)*zeta));
% 		constraints= constraints+(Del(l,k)>=u(l,1)*((ones(E,1)-supp(:,l))'*C(:,:,k)*lam(l,:,k)'-u(l,1)*zeta));
% 		constraints=constraints+(lam(l,:,k)'+N(2:end-1,:)'*gam(l,:,k)'-N(end,:)'>=0);
% 	end
% end
% options=sdpsettings('solver','cplex','verbose',0);
% optimize(constraints,objective, options);
% cvar_opt = value(objective);
% zeta=value(zeta);