clc
clear all
close all
load out_10000
% load /Users/utsav/Dropbox/Gurobi_SBB/Results/workspace_out_all17.mat
% load mat2_cvar
len_indx = length(indxed);
%     for it1=1:length(indxed)
%         m_R(it1,1) = prctile(reshape(cvar_out_R_all(4,indxed(it1),:),1000,1),95,1);
%         m_L(it1,1) = prctile(reshape(cvar_out_L_all(4,indxed(it1),:),1000,1),95,1);
%         m_det(it1,1) = prctile(reshape(cvar_out_det_R(4,indxed(it1),:),1000,1),95,1);
%     end
R_out_all = cvar_out_R_all(:,indxed(:),:);
% L_out_all = cvar_out_L_all(:,indxed(:),:);
det_out_all = cvar_out_det_R(:,indxed(:),:);

for jj=1:4
    Pctile_R(:,jj) = prctile(reshape(R_out_all(jj,:,:),length(indxed),n),95,2);
%     Pctile_L(:,jj) = prctile(reshape(L_out_all(jj,:,:),length(indxed),n),95,2);
    Pctile_det(:,jj) = prctile(reshape(det_out_all(jj,:,:),length(indxed),n),95,2);
end
x =[];
g=[];
m=0;
colors=[];
%     aa=1:4;
%     bb=6:11;
%     cc=13:15;
%     in=[aa';bb';cc'];
in= 1:17;
params_diri =[];
F=[];
pm=[1;3;5;7];
for jj=1:length(conc_all)
%     if jj==4
%         x_R = Pctile_R(in(:),jj);
% %         x_L = Pctile_L(in(:),jj);
%         x_det = Pctile_det(in(:),jj);
%     else
        x_R = Pctile_R(:,jj);
%         x_L = Pctile_L(:,jj);
        x_det = Pctile_det(:,jj);
%     x = [x;x_det;x_L;x_R];
x = [x;x_det;x_R];
%     g = [g;(m+1)*ones(size(x_det)); (m+2)*ones(size(x_L)); (m+3)*ones(size(x_R))];
g = [g;(m+1)*ones(size(x_det));  (m+3)*ones(size(x_R))];
%     colors = [colors; 1 0 0; 0 0 1; 0 0.5 0];
  colors = [colors; 1 0 0; 0 0.5 0];
%     params_diri = [params_diri; pm(jj)*ones(3*length(x_det),1)];
params_diri = [params_diri; pm(jj)*ones(2*length(x_det),1)];
    len_indx= length(x_det);
%     m=m+3;
m=m+2;
%     F = [F;repelem(["Deterministic plan "; "Randomized_{L}"; "randomized strategy"],len_indx,1)];
 F = [F;repelem(["deterministic plan"; "randomized strategy"],len_indx,1)];
end
h1=figure;
varNames = {'Strategies','alphas','cvar'};
T= table(F, params_diri, x, 'VariableNames',varNames);
h=boxchart(T.alphas,T.cvar,'GroupByColor',T.Strategies, 'MarkerStyle','.')
meandet = groupsummary(T.cvar,[ T.alphas T.Strategies],'mean');
hold on
% meandet = reshape(meandet, 3,4)
% plot([1;3;5;8], meandet,'-o')
ylabel('$95$th Percentile of the out-of-sample CVaR', 'Interpreter','Latex', 'Fontsize',12);
xlabel('Concentration parameter ($\beta$) of the Dirichlet distribution', 'Interpreter','Latex','Fontsize',12);
legend('Interpreter','Latex','Fontsize',12)

ax = h.Parent;  % axis handle
% ax.XTick =[];
ax.XTickLabel([2,4,6,8]) = {'0.1','0.3','0.5','0.8'};
ax.XTickLabel([1,3,5,7,9,10]) = {''};
set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 12);
set(gca,'TickLength',[0 0])
% writetable(T,'myData.csv','WriteVariableNames', true)
% x = [x;m_R;m_L;m_det];
% g = [g;(m+1)*ones(size(x_det)); (m+2)*ones(size(x_L)); (m+3)*ones(size(x_R))];
% colors = [colors; 1 0 0; 0 0 1; 0 0.5 0];
% boxplot(x,g)
saveas(h1,'boxplot_out_of_samplle_dirichlet','pdf')
