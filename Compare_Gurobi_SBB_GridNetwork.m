clc
clear
close all
B=2;
alpha=0.05;
p=5; r=5; 
n_iter_CD=2;partitions_sbb=1;
K_all=[100;150;200];
time_exit=1800;
tol2=0.01;
num_ins=10;
cplx=0;
addpath(genpath('~/Desktop/CPLEX_SBB_git/'))
load instances_random
K = length(K_all);
convrge_gurobi_bilinear=zeros(num_ins,K);
convrge_sbb=zeros(num_ins,K);
convrge_sbb_withoutCG=zeros(num_ins,K);
convrge_consgen=zeros(num_ins,K);
time_gurobi_bilinear=zeros(num_ins,K);
time_sbb=zeros(num_ins,K);
cvar_sbb=zeros(num_ins,K);
time_sbb_CG=zeros(num_ins,K);
cvar_sbb_CG=zeros(num_ins,K);
time_consgen=zeros(num_ins,K);
cvar_consgen=zeros(num_ins,K);
cvar_gurobi_bilinear=zeros(num_ins,K);
time_last_consgen = zeros(num_ins,K);
iter_last_consgen = zeros(num_ins,K);

for ij=1:num_ins
    rng(seed_random(ij));
    [G,dir1] = graph_generate_dir(p,r);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r);
    size_set_non_rem=size(set_non_rem,1);
    K = K_all(end);
    [cap_all,diag_cap_all,diag_cap_non_rem_all,zeta_lb,zeta_ub,F, cap_non_rem_all]=capacities(G,K,...
        E,N,set_non_rem);
    [flow_all, supp_all, supp_rem] = compute_f_l_k(cap_all,set_non_rem,set_rem,B,E,G,N,K);
    for k = 1:length(K_all)
        K = K_all(k);
        cap = cap_all(:,1:K);
        diag_cap = diag_cap_all(:,:,1:K);
        diag_cap_non_rem = diag_cap_non_rem_all(:,:,1:K);
        cap_non_rem = cap_non_rem_all(:,1:K);
        flow_k = flow_all(:,1:K);
        Gamma=0.05*K;
        l1b=Gamma*(K^0.5);
        qhat=(1/K)*ones(1,K);
        % % Bilinear Gurobi
        rel_heur_ON=true;
        time_yalmip=0;
        time=0;
        flag=0;
        [~,UB,~,~,time] = gurobi_bilinear(flow_k, Gamma,alpha,...
            K,qhat,nchoosek(E-size_set_non_rem,B),zeta_lb,zeta_ub,time_yalmip,time,l1b, ~rel_heur_ON);
        convrge_gurobi_bilinear(ij,k)=flag;
        time_gurobi_bilinear(ij,k)=time;
        cvar_gurobi_bilinear(ij,k)=UB;
        % % Spatial Branch and Bound without CG
        l_0=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
        time_yalmip=0;
        time=0;
        flag = 0;
        [time_sbb, cvar_ran, ~,~,~,flag] = sbb_withoutCG(flow_k,K, Gamma,alpha,partitions_sbb,...
            flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
            diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
        time_sbb(ij,k)=time_sbb;
        cvar_sbb(ij,k)=cvar_ran;
        convrge_sbb_withoutCG(ij,k)=flag;
        % % Constraint Generation
        time =0;
        [u, zeta,   cvar, time, J_lists, flag, ts]  =  ConsGen_GurobiBilinear(flow_k,zeta_lb,zeta_ub,...
            qhat, tol2, alpha, Gamma, K, time, time_exit,~rel_heur_ON,tol2);
        convrge_consgen(ij,k) = flag;
        time_consgen(ij,k)=time;
        cvar_consgen(ij,k)=cvar;
        time_last_consgen(ij,k) = ts(end);
        iter_last_consgen(ij,k) = length(ts);

        % % Spatial Branch and Bound with CG
        time_yalmip = 0;
        time =0;
        flag=0;
        [time, cvar_random_policy, u_random_policy,~, flag]=sbb_CG(l_0,cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
            flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
            diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
        time_sbb_CG(ij,k)=time;
        cvar_sbb_CG(ij,k)=cvar_random_policy;
    end
end
avg_time_gurobi_bilinear = squeeze(mean(time_sbb_gurobi, 1));
avg_time_consgen = squeeze(mean(time_consgen, 1));
avg_time_sbb = squeeze(mean(time_sbb, 1));
avg_time_sbb_CG = squeeze(mean(time_sbb_CG, 1));


