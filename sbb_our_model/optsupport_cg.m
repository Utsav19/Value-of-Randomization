function [fl_opt,z_in,time_yalmip,time]= optsupport_cg(alpha,Gamma,D,C,...
    l_in,E,N,B,K,zeta_lb,zeta_ub,G,qhat,~,setcom,C_nocom,set_intrdictble,time_yalmip,time,time_exit,l1b, cplx)
if time>time_exit
    return
end
NL=1;
l_opt=l_in;
z_in=zeros(E,1);
inc=1;
z_lhat=l_opt;
for e=1:E
    if ~ismember(e,setcom)
        z_in(e)=l_in(inc);
        inc=inc+1;
    end
end
while (~isempty(l_opt))
    for k=1:K
        G.Edges.Weight=(1-z_in(:,NL)).*D(:,k);
        fl_opt(NL,k)=maxflow(G,1,size(N,1));
    end
    [~,time_yalmip,time,vphi,p,pim,~]=primal_prob_cg(fl_opt,NL,z_lhat,K, zeta_lb,zeta_ub,alpha,...
        Gamma,qhat,time_yalmip,time,l1b, cplx);

    [l_opt,time_yalmip,time] =subproblem_cg(E,N,B,alpha, K, C,vphi, pim,p,...
        zeta_lb,zeta_ub,setcom,C_nocom,set_intrdictble,time_yalmip,time,z_lhat,l1b, cplx);
    if time>time_exit
        return
    end
    if isempty(l_opt)
        break
    else
        NL=NL+1;
        zp_in=zeros(E,1);
        inc=1;
        for e=1:E
            if ~ismember(e,setcom)
                zp_in(e)=l_opt(inc);
                inc=inc+1;
            end
        end
        l_in=[l_in l_opt];
        z_lhat=[z_lhat l_opt];
        z_in=[z_in zp_in];
    end

end