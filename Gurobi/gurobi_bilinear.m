% function to solve the interdictor's bilinear problem directly using gurobi's
% bilinear solver
function [u,UB,zeta,time_yalmip,time] = gurobi_bilinear(fl_opt, Gamma,alpha,...
    K,qhat,L,zeta_lb,zeta_ub,time_yalmip,time,l1b,~)
t = sdpvar(1,1);
w_plus = sdpvar(K,1);
w_minus = sdpvar(K,1);
theta=sdpvar(K,1);
theta_minus =sdpvar(K,1);
chi = sdpvar(1,1);
u=sdpvar(L,1);
zeta=sdpvar(1,1);
Del = sdpvar(L,K,'full');
objective = t;
constraints = [Del(:)>=0,theta(:)>=0,theta_minus(:)>=0, theta(:)+theta_minus(:)==chi];
U = repelem(u, 1, K);
constraints=constraints+(Del>=U.*(fl_opt-zeta));
constraints =constraints+((Gamma)*(sum(w_plus(:))+sum(w_minus(:)))+l1b*chi+zeta+...
    -w_plus+w_minus-theta+theta_minus+(1/(1-alpha))*sum(Del,1)'...
    + qhat*(w_plus-w_minus+theta-theta_minus)<=t);
constraints= constraints+[w_plus>=0,w_minus>=0,zeta<=zeta_ub,zeta>=zeta_lb,u(:)>=0,...
    ones(L,1)'*u==1];
options=sdpsettings('solver','gurobi','gurobi.NonConvex','2','Gurobi.TimeLimit','600','verbose',1);
va=optimize(constraints,objective, options);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
UB=value(objective);
zeta=value(zeta);
u=value(u);
end


