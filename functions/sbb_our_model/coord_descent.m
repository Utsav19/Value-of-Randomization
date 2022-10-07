function [ub,zeta,u,time_yalmip,time]=coord_descent(fl_opt,u,Gamma,alpha,...
    K,qhat,L,zeta_lb,zeta_ub,I,~,~,~,time_yalmip,time,l1b, cplex)
tol=1;
iter=1;
while (tol>sqrt(eps) && iter<=I)
    [UB1,zeta,time_yalmip,time] = feasiblesol_fixedu(fl_opt,u,Gamma,alpha,...
        K,qhat,L,zeta_lb,zeta_ub,time_yalmip,time,l1b, cplex);
    [UB2,u,time_yalmip,time] = feasiblesol_fixedzeta(fl_opt,zeta, Gamma,alpha, K,qhat,L,time_yalmip,time,l1b, cplex);
    tol=abs((UB2-UB1)/UB1);
    iter=iter+1;
end
ub=UB2;
end
