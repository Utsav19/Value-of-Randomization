function [l_opt,time_yalmip,time] = subproblem_cg(E,N,B,alpha, K, C,vphi_opt,...
 pi_opt,p_opt,zeta_lb,zeta_ub,setcom,C_nocom,set_intrdictble,time_yalmip,time,l_lhat,~,cplx)
if B<=6
   tolerance = 1e-5;
else 
    tolerance = 1e-4; 
end
eta=sdpvar(1,1);
scom=size(setcom,1);
l=binvar(E-scom,1);            
z_lam=sdpvar(E-scom,K,'full');
gam =sdpvar(size(N,1)-2,K,'full');
lam=sdpvar(E,K,'full');  
Del=sdpvar(K,1);
objective = vphi_opt'*Del/(1-alpha)+p_opt+pi_opt*eta;
constraints=[eta>=zeta_lb,eta<=zeta_ub, Del(:)>=0, lam(:)>=0,z_lam(:,:)>=0];
constraints=constraints+(ones(E-scom,1)'*l<=B);
for e=1:E-scom
    constraints=constraints+(z_lam(e,:)<=l(e)); 
    constraints=constraints+[z_lam(e,:)<=lam(set_intrdictble(e),:),...
            z_lam(e,:)>=lam(set_intrdictble(e),:)-(1-l(e))];
end
for k=1:K    
    constraints=constraints+(Del(k)>=ones(E,1)'*C(:,:,k)*lam(:,k)...
        -ones(E-scom,1)'*C_nocom(:,:,k)*z_lam(:,k)-eta);
    constraints=constraints+(lam(:,k)+N(2:end-1,:)'*gam(:,k)-N(end,:)'>=0);
end
for i=1:size(l_lhat,2) % for every l in L_hat for which y_l does not equal zero
    m=l_lhat(:,i);
    index0=[];
    index1=[];
    for jj=1:E-scom
        if m(jj)==0 
            index0=[index0;jj];
        else
            index1=[index1;jj]; 
        end
    end
    constraints=constraints+(size(index0,1)-sum(l(index0))+sum(l(index1))<=E-scom-1);
end
constraints=constraints+(lam(:)<=1);
if cplx
    ops=sdpsettings('solver','CPLEX','verbose',0);
else
    ops=sdpsettings('solver','GUROBI','verbose',0);
end
va=optimize(constraints,objective, ops);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
if value(objective)<-tolerance 
    l_opt=round(value(l));
else
    l_opt=[];
end
end