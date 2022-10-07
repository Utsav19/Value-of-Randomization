% % convergence of spatial branch and bound algorithm
clc 
clear 
close all
K=10;K_out=1000;p=10;r=10;n_ins=100;
alpha_set=[0.025;0.05;0.1];Gamma_set=[0.05*K; 0.1*K; 0.15*K];
B_set=[floor(0.3*p); floor(0.6*p); floor(0.9*p)];
qhat=(1/K)*ones(1,K);n_iter_CD=2;partitions_sbb=1;time_exit=7200;
sum_cvar_ran=zeros(size(B_set,1),size(alpha_set,1),size(Gamma_set,1));
cvar_ran_sbb=zeros(size(sum_cvar_ran));
cvar_det=zeros(size(sum_cvar_ran));
gap_tol2=zeros(size(sum_cvar_ran));
instance_diff=zeros(size(sum_cvar_ran));
instance_diff_L=zeros(size(sum_cvar_ran));
instance_improv =zeros(size(sum_cvar_ran));
instance_improv_L =zeros(size(sum_cvar_ran));
cvar_out_R_all =zeros(size(sum_cvar_ran));
cvar_out_L_all =zeros(size(sum_cvar_ran));
cvar_out_det_L =zeros(size(sum_cvar_ran));
cvar_out_det_R =zeros(size(sum_cvar_ran));
rel_diff_out_L =zeros(size(sum_cvar_ran));
cvar_L_all = zeros(n_ins, size(B_set,1),size(alpha_set,1),size(Gamma_set,1));
cvar_R_all = zeros(n_ins, size(B_set,1),size(alpha_set,1),size(Gamma_set,1));
cvar_D_all = zeros(n_ins, size(B_set,1),size(alpha_set,1),size(Gamma_set,1));
rel_diff=zeros(size(sum_cvar_ran));
rel_diff_L=zeros(size(sum_cvar_ran));
rel_improv=zeros(size(sum_cvar_ran));
rel_improv_L=zeros(size(sum_cvar_ran));
tol2=0.01;
tol4=0.01;
load all_instances_performance
dir1=Direction_all{1};
G = graph_generate(p,r,dir1);
ins_nconvrg=zeros(n_ins,1);
addpath('./Loizou/')

for ij=1:n_ins
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    cap_all =all_C{ij};
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub]=capacities_dirichlet(G,K,...
    E,N,set_non_rem, cap_all);
    for it1=1:size(B_set,1)
        for it2=1:size(alpha_set,1)
            for it3=1:size(Gamma_set,1)
                l1b=Gamma_set(it3)*(K^0.5);
                B=B_set(it1);
                alpha=alpha_set(it2);
                Gamma=Gamma_set(it3);
                flag=0;
                [cvar_D,deter_p,l_h,zeta_D] = minimize_wcvar(E,N,B,alpha,Gamma,qhat, K,...
                    diag_cap,set_non_rem,diag_cap_non_rem,set_rem,zeta_lb,zeta_ub,l1b);
                
                l_in=round(l_h);
                cvar_Det(ij,it1,it2,it3)=cvar_D;
                deter_plan(:,ij,it1,it2,it3)=deter_p;
                time_yalmip=0;
                time_gurobi=0;
                ustar=[];
                
                [rel_gap,gap,ub_iter,branch_iter,gap_tol4, time_tol4, cvar_R,zeta_ran,z_supp_R, u_R,time_tol2,...
                    flag,gap_tol2,time_yalmip,time_cplex,least_lb] = spatial_BB_showconvergence(l_in,cap,diag_cap,E, N, G,B,K,Gamma,alpha,partitions_sbb,...
                    flag,qhat,n_iter_CD,zeta_lb,zeta_ub,tol2,tol4,time_exit,set_non_rem,...
                    diag_cap_non_rem,set_rem,time_yalmip,time_gurobi,l1b);
                
                l_in=[ones(B,1);zeros(E-size_set_non_rem-B,1)];
                time_yalmip = 0;
                time_cplex =0;
                [fl_opt,z_supp_L,time_yalmip,time_gurobi]= Loizou_support(alpha,Gamma,cap,diag_cap,...
                    l_in,E,N,B,K,G,qhat,set_non_rem,diag_cap_non_rem,set_rem,time_yalmip,time_cplex,l1b);
                [u_L,cvar_L,zeta,time_yalmip,time_gurobi]=Loizou_Primal(fl_opt,size(z_supp_L,2),K,zeta_lb,...
                    zeta_ub,alpha, Gamma,qhat,set_non_rem,diag_cap_non_rem,time_yalmip,time_cplex,l1b);
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
                if (((cvar_D-cvar_R)/cvar_R)*100>1e-1) && (size(ustar_R,1))>1
                    rel_diff(it1,it2,it3)=rel_diff(it1,it2,it3)+((cvar_D-cvar_R)/cvar_R)*100;
                    instance_diff(it1,it2,it3)=instance_diff(it1,it2,it3)+1;
                end
                if (((cvar_D-cvar_L)/cvar_L)*100>1e-1) && (size(ustar_L,1))>1
                    rel_diff_L(it1,it2,it3)=rel_diff_L(it1,it2,it3)+((cvar_D-cvar_L)/cvar_L)*100;
                    instance_diff_L(it1,it2,it3)=instance_diff_L(it1,it2,it3)+1;
                end
                if ((cvar_D-cvar_R)/cvar_R)*100>=1   
                    rel_improv(it1,it2,it3)=rel_improv(it1,it2,it3)+((cvar_D-cvar_R)/cvar_R)*100;
                    instance_improv(it1,it2,it3)=instance_improv(it1,it2,it3)+1;
                end
                if ((cvar_D-cvar_L)/cvar_L)*100>=1 
                    rel_improv_L(it1,it2,it3)=rel_improv_L(it1,it2,it3)+((cvar_D-cvar_L)/cvar_L)*100;
                    instance_improv_L(it1,it2,it3)=instance_improv_L(it1,it2,it3)+1;
                end
                
            end
        end
        
    end
       
     dir1=Direction_all{ij+1};
     G = graph_generate(p,r,dir1);
end
save('results_section_5_2_1.mat','alpha_set','Gamma_set','B_set','cvar_D_all', ...
    'cvar_R_all','cvar_L_all','rel_improv','instance_improv','instance_diff', ...
    'rel_improv_L','instance_improv_L','instance_diff_L', 'supp_R_all', ...
    'ustar_R_all', 'supp_L_all', 'ustar_L_all');