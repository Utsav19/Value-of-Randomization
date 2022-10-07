function [q_max, obj, time] = argmax_q(Gamma, us, flows_all, zeta, qhat, J_lists,...
	n_scenarios, i, K, time)
	q = sdpvar(n_scenarios, 1);
	z = sdpvar(n_scenarios,1);
	a = 0;
	for jj=1:size(J_lists{i},2)
	l=J_lists{i}(1,jj); k=J_lists{i}(2,jj);
	a = a+q(k)*(us(l)*(flows_all(l,k)-zeta));
	end;
	objective = -a;
	constraints = (q>=0);
	constraints = constraints + [sum(q) == 1, q == qhat' + z];
	constraints = constraints + [norm(z, inf) <= Gamma, norm(z, 1) <= Gamma*(K)^0.5];
	options  = sdpsettings('solver','gurobi', 'gurobi.TimeLimit','1800','verbose',0);
	va  = optimize(constraints, objective, options);
	time = time + va.solvertime;
	q_max = value(q);
	obj = -value(objective);
end







