function [fl_opt,z_in,time_yalmip,time]= optsupport_cg_heuristic(alpha,Gamma,D,C,...
    l_in,E,N,B,K,zeta,G,qhat,time_s,setcom,C_nocom,set_intrdictble,time_exit,time_yalmip,time,l1b, cplx)

NL=1;
l_opt=l_in;
z_in=zeros(E,1);
inc=1;
for e=1:E
    if ~ismember(e,setcom)
        z_in(e)=l_in(inc);
        inc=inc+1;
    end
end
while (~isempty(l_opt))
    if time>time_exit
        return
    end
    for k=1:K
        G.Edges.Weight=(1-z_in(:,NL)).*D(:,k);
        fl_opt(NL,k)=maxflow(G,1,size(N,1));
    end
    [vphi,p,time_yalmip,time]=dual_cg_heuristic(fl_opt,NL,K, zeta,alpha, Gamma,qhat,...
        time_s,setcom,C_nocom,time_yalmip,time,l1b, cplx);
    [l_opt,time_yalmip,time] =subprob_cg_heuristic(E,N,B,alpha, K, C,vphi,p,zeta,setcom,C_nocom,...
        set_intrdictble,time_yalmip,time,l1b, cplx);
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
        z_in=[z_in zp_in];
    end
end