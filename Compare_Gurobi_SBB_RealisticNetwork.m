clc
clear 
close all
alpha=0.05;
B =2; % B=1 for results in online supplement
B_sioux=1;
n_iter_CD=2;partitions_sbb=1;
K_all=[100;150;200];
time_exit=1800;
tol2=0.01;
num_ins=10;
cplx=0;
addpath(genpath('~/Desktop/ValueofRandomization/'));
load instances_random
G_sioux = readtable('SiouxFalls.csv');
G_nobelus=readtable('nobelus.csv');
K = length(K_all);
convrge_gurobi_bilinear=zeros(num_ins,K,2);
convrge_sbb=zeros(num_ins,K,2);
convrge_sbb_withoutCG=zeros(num_ins,K,2);
convrge_consgen=zeros(num_ins,K,2);
time_gurobi_bilinear=zeros(num_ins,K,2);
time_sbb=zeros(num_ins,K,2);
cvar_sbb=zeros(num_ins,K,2);
time_sbb_CG=zeros(num_ins,K,2);
cvar_sbb_CG=zeros(num_ins,K,2);
time_consgen=zeros(num_ins,K,2);
cvar_consgen=zeros(num_ins,K,2);
cvar_gurobi_bilinear=zeros(num_ins,K,2);
time_last_consgen = zeros(num_ins,K,2);
iter_last_consgen = zeros(num_ins,K,2);
for it1 = 1:2
    if it1==1
        G_data = G_sioux;
        B= B_sioux;
    else
        G_data = G_nobelus;
    end
    s=G_data{:,1};
    t=G_data{:,2};
    G=digraph(s,t);
    H = incidence(G);
    N=full(H);
    E=size(N,2);
    set_rem =1:E;
    set_non_rem=[];
    for ij=1:num_ins
        rng(seed_random(ij));
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
            time_yalmip=0;
            time=0;
            flag=0;
            [~,UB,~,~,time] = gurobi_bilinear(flow_k, Gamma,alpha,...
                K,qhat,nchoosek(E-size_set_non_rem,B),zeta_lb,zeta_ub,time_yalmip,time,l1b);
            convrge_gurobi_bilinear(ij,k,it1)=flag;
            time_gurobi_bilinear(ij,k,it1)=time;
            cvar_gurobi_bilinear(ij,k,it1)=UB;
            % % Spatial Branch and Bound without CG
            l_0=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
            time_yalmip=0;
            time=0;
            flag = 0;
            [time_sbb, cvar_ran, ~,~,~,flag] = sbb_withoutCG(flow_k,K, Gamma,alpha,partitions_sbb,...
                flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
                diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
            time_sbb(ij,k, it1)=time_sbb;
            cvar_sbb(ij,k, it1)=cvar_ran;
            convrge_sbb_withoutCG(ij,k, it1)=flag;
            % % Constraint Generation
            time =0;
            [u, zeta,   cvar, time, J_lists, flag, ts]  =  ConsGen_GurobiBilinear(flow_k,zeta_lb,zeta_ub,...
            qhat, tol2, alpha, Gamma, K, time, time_exit,tol2);
            convrge_consgen(ij,k,it1) = flag;
            time_consgen(ij,k)=time;
            cvar_consgen(ij,k)=cvar;
            time_last_consgen(ij,k, it1) = ts(end);
            iter_last_consgen(ij,k, it1) = length(ts);

            % % Spatial Branch and Bound with CG
            time_yalmip = 0;
            time =0;
            flag=0;
            [time, cvar_random_policy, u_random_policy,~, flag]=sbb_CG(l_0,cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
            flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
            diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
            time_sbb_CG(ij,k, it1)=time;
            cvar_sbb_CG(ij,k, it1)=cvar_random_policy;
        end    
    end
end
avg_time_gurobi_bilinear = mean(time_sbb_gurobi, 1);
avg_time_consgen = mean(time_consgen, 1);
avg_time_sbb = mean(time_sbb, 1);
avg_time_sbb_CG = mean(time_sbb_CG, 1);
mat = [avg_time_gurobi_bilinear avg_time_consgen(:)...
    avg_time_sbb(:) avg_time_sbb_CG(:)];
mat = reshape(mat, [],2);
filename = 'Compare_Gurobi_SBB_realnetworks.csv';
writematrix(mat, filename)


