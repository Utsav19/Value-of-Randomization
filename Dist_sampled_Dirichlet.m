clc
clear
close all
p=15;r=15;K=15;num_ins=10;alpha=0.05;B=floor(0.3*p);
Gammas=[0.05*K;0.1*K;0.15*K];
n_iter_CD=2; partitions_sbb=1;
time_exit=7200;tol2=0.01;
addpath('./instances')
addpath('./functions')
addpath('./functions/sbb_our_model')
load instances_random
betas =[0.1,0.5,1,10000];
cplx = 1;
instance_notconvrg=zeros(num_ins,length(betas),length(Gammas));
time_all = zeros(num_ins,length(betas),length(Gammas));
for ij=1:num_ins
    rng(seed_random(ij));
    [G,dir1] = graph_generate_dir(p,r);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); 
    size_set_non_rem=size(set_non_rem,1);
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub,F]=capacities(G,K,...
        E,N,set_non_rem);
    l_0=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
    for ii=1:length(betas)
        flag=0;
        rng(seed_random(ij))
        in = 0.5*ones(K,1);
        q_true = gamrnd(in , 1);
        q_true = q_true./sum(q_true);
        qhat=q_true';
        for kk=1:length(Gammas)
            time_yalmip=0;
            time_gurobi=0;
            Gamma=Gammas(kk);
            l1b=Gamma*(K^0.5);
            time=0;
            [time, cvar_random_policy, u_random_policy,~, flag]=sbb_CG(l_0,cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
                flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
                diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
            time_all(ij,ii,kk)=time;
            instance_notconvrg(ij,ii,kk)=instance_notconvrg(ij,ii,kk)+flag;
        end
    end
end
avg_time = squeeze(mean(time_all, 1));
filename = 'dist_sampled_dirichlet.csv';
writematrix(avg_time, filename)
disp('Average computation time for distribution sampled from dirichlet', avg_time)
save('appendix_dirichlet_times.mat', 'avg_time');
