clc
clear
close all
num_ins=10;alpha=0.05;
n_iter_CD=2; partitions_sbb=1;
set_discrete=[4;16;32;64;128];
time_exit=7200;tol2=0.01;tol4=0.01;
p=10;r=10;K=10;
qhat=(1/K)*ones(1,K); % center of ambiguity set
Gammas=[0.05*K;0.1*K;0.15*K];
B=floor(0.3*p);
addpath(genpath('~/Desktop/ValueofRandomization/'));
load instances_random % random instances
load instances_VRS % instances for which VRS>0.1
cplx = 1; % 1: solver is cplx, 0: solver gurobi
instance_notconvrg=zeros(num_ins,length(Gammas), 2);
time_all=zeros(num_ins,length(Gammas),2);
subopt_heuristic = zeros(num_ins, length(Gammas), length(set_discrete),2);
time_heuristic = zeros(num_ins, length(Gammas), length(set_discrete),2);
for ij=1:num_ins
    for ran =1:2
        rng(seed_random(ij)); 
        for kk=1:length(Gammas)
            if ran==2
                rng(seed_VRS(ij,kk,1))
            end
            Gamma=Gammas(kk);
            [G,dir1] = graph_generate_dir(p,r);
            [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r);
            size_set_non_rem=size(set_non_rem,1);
            [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub,F]=capacities(G,K,...
                E,N,set_non_rem);
            l_0=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
            time_yalmip=0;
            l1b=Gamma*(K^0.5);
            time=0;
            flag=0;
            [cvar_deterministic,deter_plan,l_deter] = deterministic_wcvar(E,N,B,alpha,Gamma,qhat, K,...
                diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b, cplx);
            [time, cvar_random_policy, u_random_policy,~, flag, least_lb]=sbb_CG(round(l_deter),cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
                flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
                diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
            time_all(ij,kk, ran)=time;
            instance_notconvrg(ij,kk, ran)=instance_notconvrg(ij,kk)+flag;
            for idx=1:length(set_discrete)
                time_s=cputime;
                tic;
                time_yalmip=0;
                time=0;
                [opt_zeta,~,~,time] = heuristic_zeta(alpha, Gamma,cap,diag_cap,...
                l_0,E,N,B,K,G,qhat,time_s,set_non_rem,diag_cap_non_rem,...
                set_discrete(idx),zeta_lb,zeta_ub,set_rem,...
                time_exit,time_yalmip,time,l1b,cplx);
                subopt_heuristic(ij, idx, kk, ran) =((opt_zeta-least_lb)/opt_zeta)*100;
                time_heuristic(ij,idx, kk, ran)=time;
             end
        end
    end
end

% generate table in the paper
avg_time = reshape(mean(time_all, 1),1,[]);
avg_time_heuristic = reshape(mean(time_heuristic, 1), length(set_discrete),[]);
avg_subopt = reshape(mean(subopt_heuristic,1), length(set_discrete),[]);
[m,n] = size(avg_time_heuristic);
mat = zeros(2*m,n);
mat(1:2:end,:) = avg_time_heuristic;
mat(2:2:end,:) = avg_subopt;
mat(end+1,:) = avg_time;

filename = 'Compare_SBB_Heuristic.csv';
writematrix(mat, filename)
