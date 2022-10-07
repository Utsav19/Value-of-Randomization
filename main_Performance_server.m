clc
clear
close all
K=10;
K_out=10000;
p=10;r=10;
num_ins=100;
alpha_set=[0.05;0.1; 0.2; 0.4];
alpha_set=0.1
Gamma_set=[0.05*K; 0.1*K; 0.15*K];
Gamma_set= 0.05*K
B_set=[floor(0.3*p); floor(0.6*p); floor(0.9*p)];
B_set=floor(0.9*p)
qhat=(1/K)*ones(1,K);
n_iter_CD=2;
partitions_sbb=1;
time_exit=7200;
cvar_L_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
cvar_R_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
cvar_D_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
tol2=0.01;
addpath(genpath('~/Desktop/ValueofRandomization/'))
load all_instances_performance
ins_nconvrg=zeros(num_ins,1);
cplx=1;
for ij=1:num_ins
    dir1=Direction_all{ij};
    G = graph_generate(p,r,dir1);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    cap_all =all_C{ij};
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub]=capacities_dirichlet(G,K,...
        E,N,set_non_rem, cap_all);
    for it1=1:size(Gamma_set,1)
        for it2=1:size(B_set,1)
            for it3=1:size(alpha_set,1)
                l1b=Gamma_set(it1)*(K^0.5);
                B=B_set(it2);
                alpha=alpha_set(it3);
                Gamma=Gamma_set(it1);
                flag=0;
                [cvar_D,deter_p,l_h] = deterministic_wcvar(E,N,B,alpha,Gamma,qhat, K,...
                    diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b,cplx);
                l_0=round(l_h);
                cvar_Det(ij,it1,it2,it3)=cvar_D;
                deter_plan(:,ij,it1,it2,it3)=deter_p;
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
                cvar_D_all(ij,it1,it2,it3) = cvar_D;
                cvar_R_all(ij,it1,it2,it3) = cvar_R;
                cvar_L_all(ij,it1,it2,it3) = cvar_L;
                rr=1;
                ustar_R=[];
                supp_R=[];
                for it=1:size(u_R,1)
                    if u_R(it) >sqrt(eps)
                        supp_R(:,rr)=z_supp_R(:,it);
                        ustar_R(rr,1)=u_R(it,1);
                        rr=rr+1;
                    end
                end
                rr=1;
                ustar_L=[];
                supp_L=[];
                for it=1:size(u_L,1)
                    if u_L(it) >sqrt(eps)
                        supp_L(:,rr)=z_supp_L(:,it);
                        ustar_L(rr,1)=u_L(it,1);
                        rr=rr+1;
                    end
                end
                supp_R_all{:,ij,it1,it2,it3}=supp_R;
                ustar_R_all{:,ij,it1,it2,it3}=ustar_R;
                supp_L_all{:,ij,it1,it2,it3}=supp_L;
                ustar_L_all{:,ij,it1,it2,it3}=ustar_L;
            end
        end
    end
end
cvar{1}=cvar_R_all;
cvar{2} =cvar_L_all;
for i=1:2
    rel_diff= 100*(cvar_D - cvar{i})./cvar{i};
    idx_1 = rel_diff>=1;
    idx_0 = rel_diff<0.1;
    idx_0to1= num_ins -(idx_0+idx_1);
    num_inst_VRS_higherthan1{i} = sum(idx_1,1);
    num_inst_VRS_0{i} = sum(idx_0,1);
    num_inst_VRS_0to1{i} = sum(idx_0to1,1);
    avg_VRS{i} = sum(rel_diff,1)./sum(idx_1,1);
end



n_ins=100; n_train=100;
instance_all=[];
alpha=0.1;
qhat=(1/K)*ones(1,K);
n_iter_CD=2;partitions_sbb=1;time_exit=7200;

cvar_D= cvar_D_all(:,:,3,:);
cvar_R= cvar_R_100(:,:,3,:);
rel_SD_diff= 100*(cvar_D - cvar_R)./cvar_R;
cvar_D= cvar_D_100(:,:,3,:);
cvar_R= cvar_R_100(:,:,3,:);
indx1 = rel_SD_diff(:,3,1,1)>=1;
indx2 = rel_SD_diff(:,3,1,2)>=1;
indx3= rel_SD_diff(:,3,1,3)>=1;
A=ij(indx1);
B=ij(indx2);
C=ij(indx2);
indxed = intersect(intersect(A,B),C);
B= 9;
conc_all = [0.1;0.3;0.5;0.8];
gg=1;
cplx=1;
for i1=1:length(indxed)
    ij = indxed(i1);
    dir1=Direction_all{ij};
    G = graph_generate(p,r,dir1);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    cap_all =all_C{ij};
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub]=capacities_dirichlet(G,K,...
        E,N,set_non_rem, cap_all);
    rng('shuffle')
    s{gg}=rng;
    for jj = 1:length(conc_all)
        param_conc = conc_all(jj);
        in =  param_conc*ones(K,1);
        a = zeros(n_ins,1);
        rng(s{gg})
        q_uni = ones(K,1)*(1/K);
        for i=1:n_train
            q_true = gamrnd(in , 1);
            q_true = q_true./sum(q_true);
            h = norm(q_true-q_uni,2);
            a(i,1) = h;
        end
        b = sort(a);
        Gamma = b(ceil(0.95*n_train));Gammas(jj,ij)=Gamma;
        l1b=Gamma*(K^0.5);
        flag=0;
        [cvar_D,deter_p,l_h] = minimize_wcvar(E,N,B,alpha,Gamma,qhat, K,...
            diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b, cplx);
        l_in=round(l_h);
        cvar_Det(jj,ij)=cvar_D;
        deter_plan(:,jj,ij)=deter_p;
        time_yalmip=0;
        time_gurobi=0;
        ustar=[];
        [rel_gap,gap,ub_iter,branch_iter,gap_tol4, time_tol4, cvar_R,zeta_ran,z_supp_R, u_R,time_tol2,...
            flag,gap_tol2,time_yalmip,time_cplex,least_lb] = sbb_CG(l_in,cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
            flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol4,time_exit,set_non_rem,...
            diag_cap_non_rem,set_rem,time_yalmip,time_gurobi,l1b);
        l_in=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
        time_yalmip = 0;
        time_cplex =0;
        rr=1;
        ustar_R=[];
        supp_R=[];
        for it=1:size(u_R,1)
            if u_R(it) >sqrt(eps)
                supp_R(:,rr)=z_supp_R(:,it);
                ustar_R(rr,1)=u_R(it,1);
                rr=rr+1;
            end
        end
        rr=1;
        supp_R_all{:,jj,ij}=supp_R;
        ustar_R_all{:,jj,ij}=ustar_R;
        if (((cvar_D-cvar_R)/cvar_R)*100>1e-1) && (size(ustar_R,1))>1
            rel_diff(jj,ij)=rel_diff(jj,ij)+((cvar_D-cvar_R)/cvar_R)*100;
            instance_diff(jj,ij)=instance_diff(jj,ij)+1;
        end
        fl_R=zeros(size(supp_R,2),K);
        parfor i=1:size(supp_R,2)
            fl_R(i,:) = flows_test(G, supp_R(:,i), cap, N,K);
        end
        fl_det=zeros(1,K);
        fl_det= flows_test(G, deter_p, cap, N,K);
        rel_improv(jj,ij)=rel_improv(jj,ij)+((cvar_D-cvar_R)/cvar_R)*100;
        instance_improv(jj,ij)=instance_improv(jj,ij)+1;
        parfor i=1:K_out
            q_true = gamrnd(in , 1);
            q_true = q_true./sum(q_true);
            [cvar_out_R,zeta_ran]=test_given_q(fl_R,K,alpha,q_true',ustar_R);
            [cvar_out_det,zeta_det]=test_given_q(fl_det,K,alpha,q_true',1);
            cvar_out_ran_all(jj,i1,i) = cvar_out_R;
            cvar_out_det_R(jj,i1,i) =cvar_out_det;
        end
    end
    gg=gg+1;
end
len_indx = length(indxed);
n=10000;
for jj=1:4
    Pctile_R(:,jj) = prctile(reshape(cvar_out_R_all(jj,:,:),length(indxed),n),95,2);
    Pctile_det(:,jj) = prctile(reshape(cvar_out_det_R(jj,:,:),length(indxed),n),95,2);
end
x =[];
g=[];
m=0;
colors=[];
in= 1:17;
params_diri =[];
F=[];
pm=[1;3;5;7];
for jj=1:length(conc_all)
    x_R = Pctile_R(:,jj);
    x_det = Pctile_det(:,jj);
    x = [x;x_det;x_R];
    g = [g;(m+1)*ones(size(x_det));  (m+3)*ones(size(x_R))];
    colors = [colors; 1 0 0; 0 0.5 0];
    params_diri = [params_diri; pm(jj)*ones(2*length(x_det),1)];
    len_indx= length(x_det);
    m=m+2;
    %     F = [F;repelem(["Deterministic plan "; "Randomized_{L}"; "randomized strategy"],len_indx,1)];
    F = [F;repelem(["deterministic plan"; "randomized strategy"],len_indx,1)];
end
h1=figure;
varNames = {'Strategies','alphas','cvar'};
T= table(F, params_diri, x, 'VariableNames',varNames);
h=boxchart(T.alphas,T.cvar,'GroupByColor',T.Strategies, 'MarkerStyle','.');
hold on
meandet = groupsummary(T.cvar,T.alphas,'mean');
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
