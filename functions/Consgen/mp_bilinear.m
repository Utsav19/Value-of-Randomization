function [us, zeta,  t, obj, time] = mp_bilinear(zeta_lb, zeta_ub, alpha,...
    q_all,  flows_all, J_lists, time, rel_heur_ON)
    [n_plans,~]= size(flows_all);
    zeta = sdpvar(1);
    t = sdpvar(1);
    us = sdpvar(n_plans,1);
    objective = zeta + (1/(1-alpha))*t;
    constraints = (us>=0);   
    constraints = constraints + (zeta_lb<=zeta); 
    constraints = constraints + (zeta_ub>=zeta); 
    constraints = constraints + (sum(us) == 1); 
    constraints = constraints + (t>=0); 
    for j=1:length(J_lists)
     l=J_lists{j}(1,:); k=J_lists{j}(2,:);
    in = sub2ind(size(flows_all), l(:),k(:));
    constraints = constraints + (sum(us(l(:)).*(flows_all(in)-zeta).*q_all{j}(k(:))) <= t);
    end
    if rel_heur_ON
    ops=sdpsettings('solver','gurobi','gurobi.NonConvex','2','Gurobi.TimeLimit','1800','verbose',1);
    else
        ops=sdpsettings('solver','gurobi','gurobi.NonConvex','2','gurobi.NoRelHeurTime','0','Gurobi.TimeLimit','1800','verbose',1);
    end
    va=optimize(constraints,objective, ops);
    obj=value(objective);
    time = time + va.solvertime;
    us =value(us);
    zeta = value(zeta);
    t= value(t);
end


