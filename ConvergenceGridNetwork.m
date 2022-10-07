% Experiments for Table 2 (random instances) and Table 3 (instances with
% VRS>0.1%)
clc
clear
close all
VRS =  input(['Enter 1 to run experiments for instances for which Value of ' ...
    'randomization is greater than 0.1%, otherwise 0:\n']);
num_ins=2;alpha=0.05;
size_net =[10; 15; 20];
n_iter_CD=2; partitions_sbb=1;
time_exit=7200;tol2=0.01;tol4=0.0001;
addpath(genpath('~/Desktop/ValueofRandomization/'));
load instances_random
load instances_VRS
cplx = 1;
instance_notconvrg=zeros(num_ins,length(size_net),3);
time_tol4_all=zeros(num_ins,length(size_net),3);
time_tol2_all = zeros(num_ins,length(size_net),3);
for ii=1:length(size_net)
    p=size_net(ii);
    r = size_net(ii);
    K=size_net(ii);
    Gammas=[0.05*K;0.1*K;0.15*K];
    qhat=(1/K)*ones(1,K);
    B=floor(0.3*p);
    for kk=1:length(Gammas)
        Gamma=Gammas(kk);
        for ij=1:num_ins
            ij
            if VRS==1
                rng(seed_VRS(ij,kk,ii));
            else
                rng(seed_random(ij));
            end
            [G,dir1] = graph_generate_dir(p,r);
            [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r);
            size_set_non_rem=size(set_non_rem,1);
            [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub,F]=capacities(G,K,...
                E,N,set_non_rem);
            l_0=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
            time_yalmip=0;
            time_gurobi=0;
            l1b=Gamma*(K^0.5);
            time=0;
            flag=0;
            [cvar_deterministic,deter_plan,l_deter] = deterministic_wcvar(E,N,B,alpha,Gamma,qhat, K,...
                diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b, cplx);
            [time_tol4, cvar_random_policy, u_random_policy,time_tol2, flag]=sbb_CG(round(l_deter),cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
                flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol4,time_exit,set_non_rem,...
                diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
            time_tol4_all(ij,ii,kk)=time_tol4;
            time_tol2_all(ij,ii,kk) =time_tol2;
            instance_notconvrg(ij,ii,kk)=instance_notconvrg(ij,ii,kk)+flag;
        end
    end
end

avg_time_tol4 = squeeze(mean(time_tol4_all, 1));
avg_time_tol2 = squeeze(mean(time_tol2_all, 1));
[m,n] = size(avg_time_tol2);
mat = zeros(2*m,n);
mat(1:2:end,:) = avg_time_tol2;
mat(2:2:end,:) = avg_time_tol4;

if VRS
    filename = 'Convergence_VRS.csv';
else
    filename = 'Convergence_random.csv';
end
writematrix(mat, filename)     