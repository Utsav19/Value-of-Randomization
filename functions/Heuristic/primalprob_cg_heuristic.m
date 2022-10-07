function [u,cvar_opt,time_yalmip,time]=primalprob_cg_heuristic(fl,L,K,zeta,alpha,...
 Gamma,qhat,~,~,~,time_yalmip,time,l1b, cplx)

% For a given zeta and support set, determine the optimal solution
sdpvar mu;
u=sdpvar(L,1);
onei = ones(L,1);
t = sdpvar(1,1);
w_plus = sdpvar(K,1);
w_minus = sdpvar(K,1);
theta=sdpvar(K,1);
theta_minus=sdpvar(K,1);
chi = sdpvar(1,1);
Del = sdpvar(L,K,'full');
eta = sdpvar(L,1);
objective = t;
constraints = [Del(:)>=0,theta(:)>=0,theta_minus(:)>=0, theta(:)+theta_minus(:)==chi];

    for l=1:L
        for k=1:K
            constraints= constraints+(Del(l,k)>=fl(l,k)*u(l,1)-u(l,1)*zeta);
        end
    end
constraints =constraints+((Gamma)*(sum(w_plus(:))+sum(w_minus(:)))+l1b*chi+zeta+mu+ qhat*(w_plus-w_minus+theta-theta_minus)<=t);
for k=1:K
    constraints =constraints+(w_plus(k)-w_minus(k)+theta(k)-theta_minus(k)...
        -(1/(1-alpha))*sum(Del(:,k))+mu>=0);
end

constraints= constraints+[u(:)>=0, onei'*u==1];
constraints= constraints+[w_plus(:)>=0, w_minus(:)>=0];
if cplx
    options=sdpsettings('solver','CPLEX','verbose',0);
else
    options=sdpsettings('solver','GUROBI','verbose',0);
end
va=optimize(constraints,objective, options);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
cvar_opt = value(objective);
u = value(u);

