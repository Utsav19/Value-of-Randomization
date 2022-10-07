% load all_yy: for only random instances
% table 1 reports for random instances
% converge_ran11 reports for  instances VRS>0.1%
clc
clear all
close all
p=5;r=5;K=5;n_ins=10;alpha=0.05;B=floor(0.6*p);
Gamma=0.15*K;qhat=(1/K)*ones(1,K);n_iter_CD=2;partitions_sbb=1;
l1b=Gamma*(K^0.5);
time_exit=7200;tol4=0.01;tol2=0.01;
load ../Gurobi_SBB/instances_515
load ../Gurobi_SBB/flows_515
 s=rng(yy(1));
[G,dir1] = graph_generate_dir(p,r); % generate random graph instance
load ../Gurobi_SBB/SBB_515
for ij=1:10
    ij
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub,F]=capacities(G,K,...
        E,N,set_non_rem);   
    fl = fl_1{ij};
    time_cplex = 0;
    time_yalmip = 0;
    L = size(fl,1);
    [UB,zeta,ustar, time_yalmip,time_cplex] = Loizou_full_prob(fl, Gamma,alpha,...
        K,qhat,L,zeta_lb,zeta_ub,time_yalmip,time_cplex,l1b);
    index=find(ustar);
    ustar=ustar(index);
    fl = fl(index,:);
    L = size(ustar,1);
    [UB,zeta,time_yalmip,time_gurobi] = CG_feasible(fl,ustar, Gamma,alpha,...
    K,qhat,L,zeta_lb,zeta_ub,time_yalmip,time_cplex,l1b);
    % l_in=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
    % [fl_opt,z_in,time_yalmip,time_gurobi]= Loizou_support(alpha,Gamma,cap,diag_cap,...
    % l_in,E,N,B,K,G,qhat,set_non_rem,diag_cap_non_rem,set_rem,time_yalmip,time_cplex,l1b);
    % [u,cvar_opt,zeta,time_yalmip,time_gurobi]=Loizou_Primal(fl_opt,size(z_in,2),K,zeta_lb,...
    % zeta_ub,alpha, Gamma,qhat,set_non_rem,diag_cap_non_rem,time_yalmip,time_cplex,l1b);
    (cvar_ran_sbb(ij)-UB)
    ustar
    ustar_ins{ij}
    time_yalmip_Loizou(ij,1)=time_yalmip;
    time_cplex_Loizou(ij,1)=time_cplex;
    index=find(ustar);
    ustar=ustar(index);
    ustar_Loizou{:,ij}=ustar;
    UB_Loizou(ij)=UB;
    if ij+1> 10
        break;
    end
     s=rng(yy(ij+1));
    [G,dir1] = graph_generate_dir(p,r); % generate random graph instance

end
keyboard
% save('sBB_gurobi_515.mat','Gamma','alpha','K','B','p','r',...
%     'time_gurobi_all', 'UB_gurobi_all', 'ustar_gurobi','time_gurobi_sbb',...
%     'cvar_ran_sbb','ustar_sbb');

