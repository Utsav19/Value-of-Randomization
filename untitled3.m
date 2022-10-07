% % out_of_sample performance Section 5.2.2
clc
clear 
close all


load data_out_section_5_2_2 
% contains all instances from Section 5.2.1 
% variable 'indxed' stores indices for instances for which there
% is in_sample improvement
% s is the seed used to obtain out-of-sample data
% optimal deteerministic plan: deter_plan
% optimal randomized plan: ustar_R_all: 
% support randomized plan: supp_R_all

K=10;
K_out=10000;
p=10;
r=10;
n_ins=100; % total number of instances in Section 5.2.1
n_train=100;
n=10;
instance_all=[];
B= 9;
alpha=0.05;
qhat=(1/K)*ones(1,K);n_iter_CD=2;partitions_sbb=1;time_exit=7200;
rel_diff(4,100)=0;
instance_diff(4,100)=0;
rel_improv(4,100)=0;
instance_improv(4,100) =0;
tol2=0.01;
tol4=0.01;
conc_all = [0.1;0.3;0.5;0.8];
gg=1; % index of the seed for out-of-sample data
for ij=1:length(indxed)
    dir1=Direction_all{ij};
    G = graph_generate(p,r,dir1);
    [E,N,set_non_rem,set_rem]= graph_set_rem(G,p,r); % generate set of removable and nonremovable arcs
    size_set_non_rem=size(set_non_rem,1);
    cap_all =all_C{ij};
    [cap,diag_cap,diag_cap_non_rem,zeta_lb,zeta_ub]=capacities_dirichlet(G,K,...
        E,N,set_non_rem, cap_all);
    rng(s{ij})
    for jj = 1:length(conc_all)
        param_conc = conc_all(jj);
        in =  param_conc*ones(K,1);
        a = zeros(n_ins,1);
        rng(s{ij});
        q_uni = ones(K,1)*(1/K);
        for i=1:n_train
            q_trued = gamrnd(in , 1);
            q_trued = q_trued./sum(q_trued);
            h = norm(q_trued-q_uni,2);
            a(i,1) = h;
        end
        q_true=[];
        for i=1:n
             q_true1 = gamrnd(in , 1);
            q_true{i} = q_true1./sum(q_true1);
        end
        b = sort(a);
        Gamma = b(ceil(0.95*n_train));
        Gammas(jj,ij)=Gamma;
        l1b=Gamma*(K^0.5);
        deter_p= deter_plan(:,jj,ij);
        supp_R = supp_R_all{:,jj,ij};
        ustar_R=ustar_R_all{:,jj,ij};
        fl_R=zeros(size(supp_R,2),K);
        parfor i=1:size(supp_R,2)
            fl_R(i,:) = flows_test(G, supp_R(:,i), cap, N,K);
        end
        fl_R_all{ij,jj} = fl_R;
        
        fl_det=zeros(1,K);
        fl_det= flows_test(G, deter_p, cap, N,K);
        fl_det_all{ij,jj} =fl_det;
        parfor i=1:n
            [cvar_out_R,zeta_ran]=test_given_q(fl_R,K,alpha,q_true{i}',ustar_R);
            [cvar_out_det,zeta_det]=test_given_q(fl_det,K,alpha,q_true{i}',1);
            cvar_out_R_all(jj,ij,i) = cvar_out_R;
            cvar_out_det_R(jj,ij,i) =cvar_out_det;
        end
    end
end 

