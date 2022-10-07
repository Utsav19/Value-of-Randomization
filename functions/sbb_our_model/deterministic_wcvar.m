function [value_opt,l_opt,warm_ran] = deterministic_wcvar(E,N,B,alpha,Gamma, qhat, K,...
 C,setcom,C_nocom,set_rem,zeta_lb,zeta_ub,l1b,cplx)
scom=size(setcom,1);
l=binvar(E-scom,1);            
z_lam=sdpvar(E-scom,K,'full');
gam =sdpvar(size(N,1)-2,K,'full');
lam=sdpvar(E,K,'full');  
Del=sdpvar(K,1);
zeta=sdpvar(1,1);
t = sdpvar(1,1);
w_plus = sdpvar(K,1);
w_minus = sdpvar(K,1);
theta=sdpvar(K,1);
theta_minus =sdpvar(K,1);
chi = sdpvar(1,1);
objective = t;
constraints=[Del(:)>=0, lam(:)>=0,z_lam(:,:)>=0,theta(:)>=0,theta_minus(:)>=0,...
    theta(:)+theta_minus(:)==chi];
constraints=constraints+(ones(E-scom,1)'*l<=B);
for e=1:E-scom
    constraints=constraints+(z_lam(e,:)<=l(e)); 
    constraints=constraints+[z_lam(e,:)<=lam(set_rem(e),:),...
            z_lam(e,:)>=lam(set_rem(e),:)-(1-l(e))];
end
for k=1:K    
       constraints =constraints+((Gamma)*(sum(w_plus(:))+sum(w_minus(:)))+l1b*chi+zeta+...
        -w_plus(k)+w_minus(k)-theta(k)+theta_minus(k)+(1/(1-alpha))*Del(k,1)...
        + qhat*(w_plus-w_minus+theta-theta_minus)<=t);
         	constraints=constraints+(Del(k)>=ones(E,1)'*C(:,:,k)*lam(:,k)...
        -ones(E-scom,1)'*C_nocom(:,:,k)*z_lam(:,k)-zeta);
    constraints=constraints+(lam(:,k)+N(2:end-1,:)'*gam(:,k)-N(end,:)'>=0);
end
constraints=constraints+[lam(:)<=1,zeta>=zeta_lb,zeta<=zeta_ub,...
    w_plus(:)>=0, w_minus(:)>=0];
if cplx
    ops=sdpsettings('solver','CPLEX','verbose',0);
else
    ops=sdpsettings('solver','GUROBI','verbose',0);
end
optimize(constraints,objective, ops);
value_opt=value(objective);
l=value(l);
warm_ran=l;
l_opt=zeros(E,1);
inc=1;
for e=1:E
	if ~ismember(e,setcom)
	l_opt(e)=round(l(inc));
	inc=inc+1;
	end
end
end