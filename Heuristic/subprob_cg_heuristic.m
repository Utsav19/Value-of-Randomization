function [l_opt,time_yalmip,time] = subprob_cg_heuristic(E,N,B,...
    alpha, K, C,vphi_opt,p_opt,zeta,setcom,C_nocom,set_intrdictble,time_yalmip,time,~, cplx)
%identified interdiction action l_opt
%h=0 if no interdiction action identified
% N is node arc incidence  matrix
% B is the budget
tolerance = sqrt(eps);
scom=size(setcom,1);
l=binvar(E-scom,1);            % identify the support vector
z_lam=sdpvar(E-scom,K,'full');
gam =sdpvar(size(N,1)-2,K,'full');%dual
lam=sdpvar(E,K,'full');    %dual
Del=sdpvar(K,1);
objective = vphi_opt'*Del/(1-alpha)+p_opt;
constraints=[Del(:)>=0, lam(:)>=0,z_lam(:,:)>=0];
constraints=constraints+(ones(E-scom,1)'*l<=B);
for e=1:E-scom
    constraints=constraints+(z_lam(e,:)<=l(e)); 
    constraints=constraints+[z_lam(e,:)<=lam(set_intrdictble(e),:),...
            z_lam(e,:)>=lam(set_intrdictble(e),:)-(1-l(e))];
end
for k=1:K    
        constraints=constraints+(Del(k)>=ones(E,1)'*C(:,:,k)*lam(:,k)...
        -ones(E-scom,1)'*C_nocom(:,:,k)*z_lam(:,k)-zeta);
    constraints=constraints+(lam(:,k)+N(2:end-1,:)'*gam(:,k)-N(end,:)'>=0);
end
constraints=constraints+(lam(:)<=1);
if cplx
    options=sdpsettings('solver','CPLEX','verbose',0);
else
    options=sdpsettings('solver','GUROBI','verbose',0);
end
va=optimize(constraints,objective, options);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
if value(objective)<-tolerance 
    l_opt=round(value(l));
else
    l_opt=[];
end
end