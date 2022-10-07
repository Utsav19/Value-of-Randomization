function [ vpi, p,time_yalmip,time]= CgDual_heuristic(fl,Lp,K,zeta,...
  alpha, Gamma,qhat,~,~,~,time_yalmip,time,l1b, cplx)
fl=fl';
A=[(Gamma)*(ones(K,1)*ones(K,1)')+ones(K,1)*qhat-eye(K) (Gamma)*(ones(K,1)*ones(K,1)')-ones(K,1)*qhat+eye(K)...
ones(K,1)*qhat-eye(K) -ones(K,1)*qhat+eye(K)  l1b*ones(K,1)  -ones(K,1);
zeros(K,K) zeros(K,K)  eye(K) eye(K) -ones(K,1)  zeros(K,1);
zeros(K,K) zeros(K,K)  -eye(K) -eye(K)  ones(K,1) zeros(K,1);
zeros(2,4*K+2);
-eye(K) zeros(K,3*K+2);
zeros(K,K) -eye(K) zeros(K,2*K+2);
zeros(K,2*K) -eye(K) zeros(K,K+2);
zeros(K,3*K) -eye(K) zeros(K,2)];
s=[zeros(3*K,1);1;-1;zeros(4*K,1)];
h=[zeros(4*K+1,1);1];
y=sdpvar(K+1,Lp);
x=sdpvar(4*K+2,1);
objective= h'*x;
all_p=0;
for l=1:Lp
	B= [(eye(K))/(1-alpha) zeros(K,1); zeros(2*K,K+1); zeros(1,K) 1; zeros(1,K) -1; zeros(4*K,K+1)];
	all_p=all_p+B*y(:,l);
end
constraints=(A*x+all_p<=s);
for l=1:Lp
	W=[ -eye(K) fl(:,l)-zeta; -eye(K) zeros(K,1);zeros(1,K) -1];
	constraints=constraints+(W*y(:,l)<=0);
end
if cplx
    ops=sdpsettings('solver','CPLEX','verbose',0);
else
    ops=sdpsettings('solver','GUROBI','verbose',0);
end
va=optimize(constraints,objective, ops);
time_yalmip=time_yalmip+va.yalmiptime;
time=time+va.solvertime;
value(objective);
psi=dual(constraints(1));
for i=2:Lp+1
	sigma=dual(constraints(i));
end

vpi=psi(1:K,1);
p=psi(3*K+1)-psi(3*K+2);
