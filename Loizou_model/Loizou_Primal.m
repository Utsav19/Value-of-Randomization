function [u,cvar_opt,zeta,time_yalmip,time]=Loizou_Primal(fl,L,K,zeta_lb,...
    zeta_ub,alpha, Gamma,qhat,~,~,time_yalmip,time,l1b, cplx)
    t = sdpvar(1,1);
    w_plus = sdpvar(K,1);
    w_minus = sdpvar(K,1);
    theta=sdpvar(K,1);
    theta_minus =sdpvar(K,1);
    chi = sdpvar(1,1);
    zeta=sdpvar(1,1);
    Del = sdpvar(K,1);
    u = sdpvar(L,1);
    objective = t;
    onei = ones(L,1);
    constraints = [Del(:)>=0,theta(:)>=0,theta_minus(:)>=0, theta(:)+theta_minus(:)==chi];
    for k=1:K
        constraints=constraints+(Del(k,1)>=u'*fl(:,k)-zeta);
            constraints =constraints+((Gamma)*(sum(w_plus(:))+sum(w_minus(:)))+l1b*chi+zeta...
        -w_plus(k)+w_minus(k)-theta(k)+theta_minus(k)+(1/(1-alpha))*Del(k,1)...
        + qhat*(w_plus-w_minus+theta-theta_minus)<=t);
    end
    constraints= constraints+[ u(:)>=0, onei'*u==1];
    constraints= constraints+[w_plus>=0,w_minus>=0,zeta<=zeta_ub,zeta>=zeta_lb];
    if cplx
        options=sdpsettings('solver','CPLEX','verbose',0);
    else
        options=sdpsettings('solver','GUROBI','verbose',0);
    end
    va=optimize(constraints,objective, options);
    time_yalmip=time_yalmip+va.yalmiptime;
    time=time+va.solvertime;
    cvar_opt=value(objective);
    zeta=value(zeta);
    u = value(u);