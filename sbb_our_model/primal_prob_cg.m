function [ ell_yzero,time_yalmip,time,vphi,p,pi,LB_obj]=primal_prob_cg(fl,Lp,z_in,...
	K,zeta_lb,zeta_ub,alpha, Gamma,qhat,time_yalmip,time,l1b, cplex)
fl=fl';
A=[(Gamma)*(ones(K,1)*ones(K,1)')+ones(K,1)*qhat-eye(K) (Gamma)*(ones(K,1)*ones(K,1)')-ones(K,1)*qhat+eye(K)...
ones(K,1)*qhat-eye(K) -ones(K,1)*qhat+eye(K) ones(K,1) l1b*ones(K,1)  -ones(K,1);
zeros(K,K) zeros(K,K)  eye(K) eye(K) zeros(K,1) -ones(K,1)  zeros(K,1);
zeros(K,K) zeros(K,K)  -eye(K) -eye(K) zeros(K,1) ones(K,1) zeros(K,1);
zeros(2,4*K+3);
zeros(1,4*K) -1 0  0 ;
zeros(1,4*K) 1 0  0 ;
zeros(Lp,4*K) ones(Lp,1) zeros(Lp,1)  zeros(Lp,1);
zeros(Lp,4*K) -ones(Lp,1) zeros(Lp,1)  zeros(Lp,1);
zeros(1,4*K) 1 0  0;
zeros(1,4*K) -1 0  0;
-eye(K) zeros(K,3*K+3);
zeros(K,K) -eye(K) zeros(K,2*K+3);
zeros(K,2*K) -eye(K) zeros(K,K+3);
zeros(K,3*K) -eye(K) zeros(K,3)];
s=[zeros(3*K,1);1;-1;0;0;zeta_ub*ones(Lp,1);-zeta_lb*ones(Lp,1);zeta_ub;-zeta_lb;zeros(4*K,1)];
h=[zeros(4*K+2,1);1];
y=sdpvar(K+2,Lp,'full');
x=sdpvar(4*K+3,1);
objective= h'*x;
onel=zeros(Lp,1);
all_p=0;
for l=1:Lp
	onel(l)=1;
	B= [(eye(K))/(1-alpha) zeros(K,2); zeros(2*K,K+2); zeros(1,K) 1 0; zeros(1,K) -1 0; zeros(1,K) 0 1; zeros(1,K) 0 -1;zeros(Lp,K) zeta_ub*onel -onel; zeros(Lp,K) -zeta_lb*onel onel; zeros(4*K+2,K+2)];
	all_p=all_p+B*y(:,l);
	onel=zeros(Lp,1);
end
constraints=(A*x+all_p<=s);
for l=1:Lp
	W=[ -eye(K) fl(:,l) -ones(K,1); -eye(K) zeros(K,1) zeros(K,1);zeros(1,K) -1 0;zeros(1,K) zeta_lb -1; zeros(1,K) -zeta_ub 1];
	constraints=constraints+(W*y(:,l)<=0);
end
if cplex
    ops=sdpsettings('solver','CPLEX','verbose',0);
else
    ops=sdpsettings('solver','GUROBI','verbose',0);
end
va=optimize(constraints,objective, ops);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
LB_obj=value(objective);
x=value(x);
y=value(y);
ell_yzero=[];
for i=1:Lp
	hh=y(:,i);
	if norm(hh)>=sqrt(eps)
		ell=z_in(:,i);
		ell_yzero=[ell_yzero ell];
	end
end
ell_yzero=round(ell_yzero);
psi=dual(constraints(1));
vphi=psi(1:K,1);
p=psi(3*K+1)-psi(3*K+2);
pi=psi(3*K+3)-psi(3*K+4);









