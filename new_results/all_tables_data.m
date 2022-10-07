clc
clear all
n_ins=100
cvar_L = zeros(n_ins, 3,3,2);
cvar_R = zeros(n_ins, 3,3,2);
cvar_D= zeros(n_ins, 3,3,2);
load insample_VRS_apr25_2rows_first50
cvar_L(1:50,:,:,:) = cvar_L_all(1:50,:,:,:);
cvar_R(1:50,:,:,:) = cvar_R_all(1:50,:,:,:);
cvar_D(1:50,:,:,:) = cvar_D_all(1:50,:,:,:);
load insample_VRS_apr25_2rows_last50
cvar_L(51:100,:,:,:) = cvar_L_all(51:100,:,:,:);
cvar_R(51:100,:,:,:) = cvar_R_all(51:100,:,:,:);
cvar_D(51:100,:,:,:) = cvar_D_all(51:100,:,:,:);
% load in_sample_VRS_results_rev2
instance_SD_nobenefit =0;
instance_SD_improv = 0;
rel_SD_diff= 100*(cvar_D - cvar_R)./cvar_R;
instance_SD_nobenefit = sum((rel_SD_diff<0.1));
instance_SD_mid = sum(rel_SD_diff>=0.1 &rel_SD_diff<1 );
instance_SD_improv = sum((rel_SD_diff>=1));
for i=1:3
    for j=1:3
        for k=1:2
            rel_SD= 100*(cvar_D(:,i,j,k) - cvar_R(:,i,j,k))./cvar_R(:,i,j,k);
            ind_SD_improv = rel_SD>=1;
            if sum(ind_SD_improv)>0
            improve_SD_all(:, i,j,k) = sum(rel_SD(ind_SD_improv))/(sum(ind_SD_improv));
            end
        end
    end
end

improve_SD_all

instance_L_nobenefit =0;
instance_L_improv = 0;
rel_L_diff = 100*(cvar_D - cvar_L)./cvar_L;
instance_L_nobenefit = sum((rel_L_diff<0.1));
instance_L_improv = sum((rel_L_diff>=1),1);

for i=1:3
    for j=1:3
        for k=1:2
            rel_L= 100*(cvar_D(:,i,j,k) - cvar_L(:,i,j,k))./cvar_L(:,i,j,k);
            ind_L_improv = rel_L>=1;
            improve_L_all(:, i,j,k) = sum(rel_L(ind_L_improv))/(sum(ind_L_improv));
        end
    end
end
