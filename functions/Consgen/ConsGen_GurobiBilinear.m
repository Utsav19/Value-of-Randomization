function  [u, zeta,  UB, time, J_lists, flag,  t_all]  =  ConsGen_GurobiBilinear(flows_all, zeta_lb, zeta_ub,...
    qhat, eps_opt, alpha, Gamma, K, time, time_exit,rel_heur_ON,~)
    flag=0;
    J_lists = {};
    q_all = [];
    [n_plans,n_scenarios]= size(flows_all);
    u = zeros(n_plans,1);
    for i=1:n_plans
        if mod(i,10)==0
            flag =1;
        end
        u(i,1) = 1;
        J_lists{i} = [];
        ks = (1:n_scenarios);
        J_lists{i}(:,end+1:end+n_scenarios) = [i*ones(1,n_scenarios);ks];
        zeta = 0;
        [q_max, ~, time] = argmax_q(Gamma, u, flows_all, zeta, qhat,  J_lists, n_scenarios,i, K, time);
        q_all{i} = q_max;
    end
    gap =1;
    t_all =[];
    while(gap>=eps_opt)
        t1=time;
        if time > time_exit
            break;
        end
        [u, zeta,  t, LB, time] = mp_bilinear(zeta_lb, zeta_ub, alpha, q_all,  flows_all,...
            J_lists, time, rel_heur_ON);
        ls = find(u > 0);
        J_lists{i} = [];
        for j = 1:length(ls)
            l = ls(j);
            ks = find(flows_all(l, :) - zeta > 0);
            tmp = length(ks);
            J_lists{i}(:,end+1:end+tmp) = [l*ones(1,tmp);ks];
        end
        [q_max, obj, time] = argmax_q(Gamma, u, flows_all, zeta, qhat, J_lists,...
        n_scenarios, i, K, time);
        UB=zeta + obj/(1-alpha);
        q_all{i} = q_max;
        i = i + 1;
        gap = abs(100*(UB-LB)/UB);
        t_all(end+1) = time-t1;
    end
    sdpvar t;
    us = u;
    objective = zeta + (1/(1-alpha))*t;
    constraints = (t>=0);
    for j=1:length(J_lists)
        l=J_lists{j}(1,:); k=J_lists{j}(2,:);
        in = sub2ind(size(flows_all), l(:),k(:));
        constraints = constraints + (sum(us(l(:)).*(flows_all(in)-zeta).*q_all{j}(k(:))) <= t);
    end
    ops=sdpsettings('solver','gurobi','gurobi.method','1','gurobi.TimeLimit','1800','verbose',1);
    va=optimize(constraints,objective, ops);
    UB=value(objective);
    time = time + va.solvertime;
end