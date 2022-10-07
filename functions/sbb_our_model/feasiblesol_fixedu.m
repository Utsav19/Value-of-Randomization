function [UB,zeta,time_yalmip,time] = feasiblesol_fixedu(fl_opt,u, Gamma,alpha,...
    K,qhat,L,zeta_lb,zeta_ub,time_yalmip,time,l1b, cplx)
t = sdpvar(1,1);
w_plus = sdpvar(K,1);
w_minus = sdpvar(K,1);
theta=sdpvar(K,1);
theta_minus =sdpvar(K,1);
chi = sdpvar(1,1);
zeta=sdpvar(1,1);
Del = sdpvar(L,K,'full');
objective = t;
constraints = [Del(:)>=0,theta(:)>=0,theta_minus(:)>=0, theta(:)+theta_minus(:)==chi];
for l=1:L
    for d=1:K
        constraints=constraints+(Del(l,d)>=u(l,1)*fl_opt(l,d)-u(l,1)*zeta);
    end
end

for k=1:K
    constraints =constraints+((Gamma)*(sum(w_plus(:))+sum(w_minus(:)))+l1b*chi+zeta...
        -w_plus(k)+w_minus(k)-theta(k)+theta_minus(k)+(1/(1-alpha))*sum(Del(:,k))...
        + qhat*(w_plus-w_minus+theta-theta_minus)<=t);
end
constraints= constraints+[w_plus>=0,w_minus>=0,zeta<=zeta_ub,zeta>=zeta_lb];
if cplx
    ops=sdpsettings('solver','CPLEX','verbose',0);
else
    ops=sdpsettings('solver','GUROBI','verbose',0);
end
va=optimize(constraints,objective, ops);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
UB=value(objective);
zeta=value(zeta);
end


