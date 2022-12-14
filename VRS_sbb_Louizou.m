clc
clear
close all
K=10;
K_out=10000;
p=10;r=10;
num_ins=100;
alpha_set=[0.05;0.1; 0.2; 0.4];
Gamma_set=[0.05*K; 0.1*K; 0.15*K];
B_set=[floor(0.3*p); floor(0.6*p); floor(0.9*p)];
qhat=(1/K)*ones(1,K);
n_iter_CD=2;
partitions_sbb=1;
time_exit=7200;
cvar_loizou_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
cvar_sbb_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
cvar_det_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
tol2=0.01;
addpath(genpath('~/Desktop/ValueofRandomization/'))
load instances_VRS
ins_nconvrg=zeros(num_ins,1);
cplx=1;
%% Insample VRS from our model and Loizou
for ij=1:num_ins
    dir1=Direction_all{ij};
    G = graph_generate(p,r,dir1);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    cap_all =all_C{ij};
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub]=capacities_dirichlet(G,K,...
        E,N,set_non_rem, cap_all);
    for it1=1:size(alpha_set,1)
        for it2=1:size(B_set,1)
            for it3=1:size(Gamma_set,1)
                l1b=Gamma_set(it3)*(K^0.5);
                Gamma=Gamma_set(it3);
                B=B_set(it2);
                alpha=alpha_set(it1);
                flag=0;
                [cvar_D,deter_p,l_h] = deterministic_wcvar(E,N,B,alpha,Gamma,qhat, K,...
                    diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b,cplx);
                l_0=round(l_h);
                time_yalmip=0;
                time=0;
                [~, cvar_R, u_R,~,flag, leastLB,z_supp_R]...
                    = sbb_CG(l_0,cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
                    flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
                    diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
                l_in=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
                time_yalmip = 0;
                time =0;
                [fl_opt,z_supp_L,~,time]= Loizou_support(alpha,Gamma,cap,diag_cap,...
                    l_0,E,N,B,K,G,qhat,set_non_rem,diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
                [u_L,cvar_L,zeta,time_yalmip,time]=Loizou_Primal(fl_opt,size(z_supp_L,2),K,zeta_lb,...
                    zeta_ub,alpha, Gamma,qhat,set_non_rem,diag_cap_non_rem,time_yalmip,time,l1b,cplx);
                cvar_det_all(ij,it1,it2,it3) = cvar_D;
                cvar_sbb_all(ij,it1,it2,it3) = cvar_R;
                cvar_loizou_all(ij,it1,it2,it3) = cvar_L;
            end
        end
    end
end

% % generate table in the paper
cvar{1}=cvar_sbb_all;
cvar{2} =cvar_loizou_all;
mat=[];
for i=1:2
    rel_diff= 100*(cvar_det_all - cvar{i})./cvar{i};
    idx_1 = rel_diff>=1;
    idx_0 = rel_diff<0.1;
    VRS_high1 = sum(idx_1,1);
    VRS0 = sum(idx_0,1);
    VRS_less1 = num_ins-(VRS_high1+VRS0);
    avg_VRS = sum(rel_diff,1)./sum(idx_1,1);
    mat=[mat VRS0(:) VRS_less1(:)...
        VRS_high1(:) avg_VRS(:)];
end

filename = 'Insample_VRS_sbb_Loizou.csv';
writematrix(mat, filename);


%% out of sample test

% determine instances for which VRS>=1 for all Gammas
ij=1:1:100;
alpha=0.1;
B=9;
idx1 = find(alpha_set==alpha);
idx2 = find(B_set==B);
cvar_D= squeeze(cvar_det_all(:,idx1,idx2,:));
cvar_R= squeeze(cvar_sbb_all(:,idx1,idx2,:));
rel_SD_diff= 100*(cvar_D - cvar_R)./cvar_R;
idx = rel_SD_diff(:,1)>=1;
A=ij(idx);
for i=1:size(Gamma_set)-1
    idx_next = rel_SD_diff(:,i+1)>=1;
    Ap=ij(idx_next);
    inst_improv =intersect(A,Ap);
    A=inst_improv;
end

n_train=100;
instance_all=[];
alpha=0.1;
qhat=(1/K)*ones(1,K);
load seeds_out_sample.mat
rel_diff = zeros(length(inst_improv),length(betas));
betas = [0.1;0.3;0.5;0.8];
%% out of sample VRS
for i1=1:length(inst_improv)
    ij = inst_improv(i1);
    dir1=Direction_all{ij};
    G = graph_generate(p,r,dir1);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    cap_all =all_C{ij};
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub]=capacities_dirichlet(G,K,...
        E,N,set_non_rem, cap_all);
    for jj = 1:length(betas)
        beta = betas(jj);
        in =  beta*ones(K,1);
        dist_all = zeros(n_train,1);
        rng(seeds_dirichlet{i1});
        q_uni = ones(K,1)*(1/K);
        for i=1:n_train
            q_true = gamrnd(in , 1);
            q_true = q_true./sum(q_true);
            dist_qs = norm(q_true-q_uni,2);
            dist_all(i,1) = dist_qs;
        end
        sorted_gam = sort(dist_all);
        Gamma = sorted_gam(ceil(0.95*n_train));
        l1b=Gamma*(K^0.5);
        flag=0;
        [cvar_D,deter_p,l_h] = deterministic_wcvar(E,N,B,alpha,Gamma,qhat, K,...
            diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b, cplx);
        l_0=round(l_h);
        time_yalmip=0;
        time=0;
        ustar=[];
        [~, cvar_R, u_R,~,flag, leastLB,z_supp_R] = sbb_CG(l_0,cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
            flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol2,time_exit,set_non_rem,...
            diag_cap_non_rem,set_rem,time_yalmip,time,l1b, cplx);
        ustar_R=[];
        supp_R=[];
        for it=1:size(u_R,1)
            if u_R(it) >sqrt(eps)
                supp_R(:,end+1)=z_supp_R(:,it);
                ustar_R(end+1,1)=u_R(it,1);
            end
        end
        rel_diff(i1,jj)=rel_diff(i1,jj)+((cvar_D-cvar_R)/cvar_R)*100;
        flow_ran=zeros(size(supp_R,2),K);
        parfor i=1:size(supp_R,2)
            flow_ran(i,:) = flows_test(G, supp_R(:,i), cap, N,K);
        end
        flow_det= flows_test(G, deter_p, cap, N,K);
        parfor i=1:K_out
            q_true = gamrnd(in , 1);
            q_true = q_true./sum(q_true);
            [cvar_out_R,zeta_ran]=test_given_q(flow_ran,K,alpha,q_true',ustar_R);
            [cvar_det,zeta_det]=test_given_q(flow_det,K,alpha,q_true',1);
            cvar_out_ran_all(i,i1,jj) = cvar_out_R;
            cvar_out_det_all(i,i1,jj) =cvar_det;
        end
    end
end

% % generate table VRS for different values of concentration parameter
filename = 'insampleVRS_with_conc_param.csv';
avg_rel_diff =  squeeze(mean(rel_diff,1));
writematrix(avg_rel_diff, filename);

len_indx = length(inst_improv);
percentile_cvar = 95;
mat_prctile=[];F=[];
prctile_ran=zeros(len_indx,length(betas));
prctile_det=zeros(len_indx,length(betas));
for jj=1:4
    prctile_ran(:,jj)= prctile(cvar_out_ran_all(:,:,jj),percentile_cvar,1);
    prctile_det(:,jj) = prctile(cvar_out_det_all(:,:,jj),percentile_cvar,1);
    mat_prctile = [mat_prctile;prctile_det(:); prctile_ran(:)];
    F = [F;repelem(["deterministic plan"; "randomized strategy"],...
    len_indx,1)];
end
xs=[1;3;5;7];
beta_all = repelem(xs, 2*len_indx,1);
h1=figure;
varNames = {'Strategies','alphas','cvar'};
T= table(F,  beta_all, mat_prctile, 'VariableNames',varNames);
h=boxchart(T.alphas,T.cvar,'GroupByColor',T.Strategies, 'MarkerStyle','.');
ylabel('$95$th Percentile of the out-of-sample CVaR', 'Interpreter','Latex', 'Fontsize',15);
xlabel('Concentration parameter ($\beta$) of the Dirichlet distribution', 'Interpreter','Latex','Fontsize',15);
legend('Interpreter','Latex','Fontsize',15)
ax = h.Parent;  % axis handle
ax.XTickLabel([2,4,6,8]) = {'0.1','0.3','0.5','0.8'};
ax.XTickLabel([1,3,5,7,9,10]) = {''};
set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 15);
set(gca,'fontsize',15)
set(gca,'TickLength',[0 0])
saveas(h1,'boxplot_outsample_alpha_pointone','pdf')
