function [l_opt,time_yalmip,time] = Loizou_subproblem(E,N,B,...
    ~, K, C,vphi_opt,p,setcom,C_nocom,set_intrdictble,time_yalmip,time,~, cplx)
tolerance = sqrt(eps);
scom=size(setcom,1);
l=binvar(E-scom,1);            
z_lam=sdpvar(E-scom,K,'full');
gam =sdpvar(size(N,1)-2,K,'full');
lam=sdpvar(E,K,'full');    
Del=sdpvar(K,1);
objective = vphi_opt'*Del+p;
constraints=[Del(:)>=0, lam(:)>=0,z_lam(:,:)>=0];
constraints=constraints+(ones(E-scom,1)'*l<=B);
for e=1:E-scom
    constraints=constraints+(z_lam(e,:)<=l(e)); 
    constraints=constraints+[z_lam(e,:)<=lam(set_intrdictble(e),:),...
            z_lam(e,:)>=lam(set_intrdictble(e),:)-(1-l(e))];
end
for k=1:K    
    constraints=constraints+(Del(k,1)>=ones(E,1)'*C(:,:,k)*lam(:,k)...
        -ones(E-scom,1)'*C_nocom(:,:,k)*z_lam(:,k));
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