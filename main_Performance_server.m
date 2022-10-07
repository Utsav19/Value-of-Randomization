clc
clear
close all
K=10;K_out=10000;p=10;r=10;num_ins=100;
alpha_set=[0.05;0.1; 0.2; 0.4];Gamma_set=[0.05*K; 0.1*K; 0.15*K];
B_set=[floor(0.3*p); floor(0.6*p); floor(0.9*p)];
qhat=(1/K)*ones(1,K);n_iter_CD=2;partitions_sbb=1;time_exit=7200;
cvar_L_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
cvar_R_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
cvar_D_all = zeros(num_ins, size(Gamma_set,1),size(B_set,1),size(alpha_set,1));
tol2=0.01;
load all_instances_performance
ins_nconvrg=zeros(num_ins,1);
addpath('./instances')
addpath('./functions')
addpath('./functions/sbb_our_model')
addpath('./gen_network/')
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
                [cvar_D,deter_p,l_h] = minimize_wcvar(E,N,B,alpha,Gamma,qhat, K,...
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
load 'All Results'/sbb_Loizou_insample.mat
num_ins=100;
cvar{1}=cvar_R;
cvar{2} =cvar_L;
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
disp('Number of instances with VRS higher than 1')
disp(num_inst_VRS_higherthan1{1});
disp(num_inst_VRS_higherthan1{2});
disp('Number of instances with VRS in [0,0.1)')
disp(num_inst_VRS_0{1});
disp(num_inst_VRS_0{2});
disp('Number of instances with VRS in [0.1, 1]')
disp(num_inst_VRS_0to1{1});
disp(num_inst_VRS_0to1{2})
disp('Average VRS')
disp(avg_VRS{1});
disp(avg_VRS{2});

