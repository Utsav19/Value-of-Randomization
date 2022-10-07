function [optimal,zeta,time_yalmip,time] = heuristic_zeta(alpha, Gamma,cap,diag_cap,...
            l_0,E,N,B,K,G,qhat,time_s,set_non_rem,diag_cap_non_rem,...
            Np,zeta_lb,zeta_ub,set_rem,...
            time_exit,time_yalmip,time,l1b,cplx)
h=(zeta_ub-zeta_lb)/Np;
d=zeta_lb:h:zeta_ub;
d=d';
opt = [];
for i=1:Np+1
	zeta=d(i);
	[fl_opt,z_supp,time_yalmip,time]=optsupport_cg_heuristic(alpha,Gamma,...
		cap,diag_cap,l_0,E,N,B,K,zeta,G,qhat,time_s,set_non_rem,diag_cap_non_rem,set_rem,...
		time_exit,time_yalmip,time,l1b, cplx);
	[~,opt(end+1),time_yalmip,time]=primalprob_cg_heuristic(fl_opt,...
		size(z_supp,2),K,zeta, alpha,Gamma,qhat,time_s,set_non_rem, diag_cap_non_rem,time_yalmip,...
		time,l1b, cplx);
end
optimal=opt(1);
for i=2:Np+1
	if opt(i)<optimal
		optimal=opt(i);
		zeta=d(i);
	end
end