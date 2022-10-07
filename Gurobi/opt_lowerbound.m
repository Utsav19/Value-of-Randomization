function [u,cvar_opt,zeta,time_yalmip,time]=opt_lowerbound(fl,L,K,zeta_lb,...
    zeta_ub,alpha, Gamma,qhat,time_yalmip,time,l1b, cplex)

u=sdpvar(L,1);
onei = ones(L,1);
t = sdpvar(1,1);
w_plus = sdpvar(K,1);
w_minus = sdpvar(K,1);
theta=sdpvar(K,1);
theta_minus=sdpvar(K,1);
chi = sdpvar(1,1);
zeta = sdpvar(1,1);
Del = sdpvar(L,K,'full');
eta = sdpvar(L,1);

objective = t;
constraints = [Del(:)>=0,theta(:)>=0,theta_minus(:)>=0, theta(:)+theta_minus(:)==chi];

for l=1:L
    for k=1:K
        constraints= constraints+(Del(l,k)>=fl(l,k)*u(l,1)-eta(l,1));
    end
end
for k=1:K
    constraints =constraints+((Gamma)*(sum(w_plus(:))+sum(w_minus(:)))+l1b*chi+zeta+...
        -w_plus(k)+w_minus(k)-theta(k)+theta_minus(k)+(1/(1-alpha))*sum(Del(:,k))...
        + qhat*(w_plus-w_minus+theta-theta_minus)<=t);
end
for j=1:L
    constraints= constraints+[eta(j)>=u(j)*zeta_lb, eta(j)<=u(j)*zeta_ub,...
    eta(j)>=zeta+zeta_ub*(u(j)-1), eta(j)<=zeta+zeta_lb*(u(j)-1)];
end
constraints= constraints+[ u(:)>=0, onei'*u==1];
constraints = constraints + (onei'*eta==zeta); 
constraints= constraints+[w_plus(:)>=0, w_minus(:)>=0,zeta>=zeta_lb,zeta<=zeta_ub];

if cplex
    ops=sdpsettings('solver','CPLEX','verbose',0);
else
    ops=sdpsettings('solver','GUROBI','verbose',0);
end
va=optimize(constraints,objective, ops);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
cvar_opt = value(objective);
u = value(u);
zeta=value(zeta);
