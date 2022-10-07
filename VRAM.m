clc
clear
close all
p=10;r=10;K=10;num_ins=10;B=0.3*p;
Gamma=0.1*K;qhat=(1/K)*ones(1,K);n_iter_CD=2;partitions_sbb=1;
l1b=Gamma*(K^0.5);
time_exit=7200;tol2=0.01;
alphas= [0;0.05;0.1;0.2;0.3;0.4];
ins_nconvrg=zeros(num_ins,length(alphas));
num_ins=10;
cplx=1;
addpath(genpath('~/Desktop/ValueofRandomization/'))
load instances_VRS
ins_seed= seed_VRS(:,2,1);
cvar_riskN= zeros(num_ins,length(alphas));
cvar_sbb= zeros(num_ins,length(alphas));
for ij=1:num_ins
    rng(ins_seed(ij));
    [G,~] = graph_generate_dir(p,r);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    l_0=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub,F]=capacities(G,K,...
        E,N,set_non_rem);
    time_yalmip = 0;
    alpha=0;
    time=0;
    [fl_opt,supp,time_yalmip,time]= Loizou_support(alpha,Gamma,cap,diag_cap,...
        l_0,E,N,B,K,G,qhat,set_non_rem,diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
    [u,cvar,zeta,time_yalmip,time]=Loizou_Primal(fl_opt,size(supp,2),K,zeta_lb,...
        zeta_ub,alpha, Gamma,qhat,set_non_rem,diag_cap_non_rem,time_yalmip,time,l1b, cplx);
    for ii =1:length(alphas)
        alpha=alphas(ii);
        fl=zeros(size(supp,2),K);
        for i=1:size(supp,2)
            fl(i,:) = flows_test(G, supp(:,i), cap, N,K);
        end
        L=length(u);
        time_yalmip=0;
        [UB,zeta,time_yalmip,~] = CG_feasible(fl,u, Gamma,alpha,...
            K,qhat,L,zeta_lb,zeta_ub,time_yalmip,time,l1b, cplx);
        cvar_riskN(ij,ii)=UB;
        flag=0;
        time=0;
        [cvar_deterministic,deter_plan,l_deter] = deterministic_wcvar(E,N,B,alpha,Gamma,qhat, K,...
            diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b, cplx);
        [~, cvar_random_policy, u_random_policy,~, flag]=sbb_CG(round(l_deter),cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
            flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
            diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
        cvar_sbb(ij,ii) = cvar_random_policy;
    end
end
vram1 =  100*(cvar_riskN-cvar_sbb)./cvar_sbb;
vram = vram1(:,2:end);
alphas = 1:1:length(alphas)-1;
alpha_copy=repelem(alphas,num_ins,1);
h1=figure;
boxchart(vram)
set(gca,'XTickLabel',{'0.05'; '0.1'; '0.2'; '0.3'; '0.4'});
ylabel('VRAM (in \%)', 'Interpreter','Latex', 'Fontsize',15);
set(gca,'fontsize',15)
xlabel('Risk aversion parameter ($\alpha$)', 'Interpreter','Latex','Fontsize',15);
saveas(h1,'boxplot_outsample_vram','pdf')
% save('VRAM.mat','vram');