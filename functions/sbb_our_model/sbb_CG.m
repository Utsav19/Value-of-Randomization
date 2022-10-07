function [time_tol4, bestUB, ustar,time_tol2,flag, leastLB,opt_supp] = sbb_CG(l_in, D,C,E, N, G,B,K,Gamma,alpha,...
    ~,flag,qhat,n_iter_CD,zeta_lb,zeta_ub,...
    tol2,tol4,time_exit,set_notintrdictble,diag_capnotintrdct,set_intrdictble,time_yalmip,time_gurobi,l1b, cplx)
time_s=cputime;
time_tol4=100;
time_tol2=100;
import java.util.ArrayDeque
que = ArrayDeque();
bestUB=100;
z_lb=zeta_lb;
[fl_opt,z_supp,time_yalmip,time_gurobi]=optsupport_cg(alpha,Gamma,D,C,l_in,E,N,B,K,z_lb,zeta_ub,G,qhat,time_s,...
    set_notintrdictble,diag_capnotintrdct,set_intrdictble,time_yalmip,time_gurobi,time_exit,l1b,cplx);
if time_gurobi>time_exit
    return
end
[u,LB,zeta,time_yalmip,time_gurobi]=opt_lowerbound(fl_opt,size(z_supp,2),K,z_lb,zeta_ub, alpha, Gamma,qhat,...
    time_yalmip,time_gurobi,l1b,cplx);
vec=[ z_lb; zeta_ub; LB; zeta;size(z_supp,2);size(fl_opt,1);u(:);fl_opt(:);z_supp(:)];
que.addFirst(vec);
n_iter=1;
while (~isEmpty(que))
    leastLB=1e6;
    s_q=size(que);
    quedum = ArrayDeque();
    for j=1:s_q % add the node with least lower bound to the front
        sl=que.removeFirst();
        if sl(3)<leastLB
            leastLB=sl(3);
            quedum.addFirst(sl)
        else
            quedum.addLast(sl);
        end
    end
    que=quedum;
    s=que.removeFirst();
    que=quedum;
    [UB, zeta_Rj,uf,time_yalmip,time_gurobi]=coord_descent(reshape(s(7+s(5):6+s(5)+s(6)*K),...
        s(6),K),s(7:6+s(5)),Gamma,alpha,K,qhat,s(5),s(1),s(2),...
        n_iter_CD,time_s,set_notintrdictble,diag_capnotintrdct,time_yalmip,time_gurobi,l1b,cplx);
    if UB<bestUB
        bestUB=UB;
        z_best=zeta_Rj;
        opt_supp=reshape(s(end-E*s(5)+1:end),E,s(5));
        ustar=uf;
    end
    % PRUNING
    s_que=size(que);
    for j=1:s_que
        g=que.removeFirst();
        if g(3)< bestUB || abs(g(3)-bestUB)<sqrt(eps)
            que.addLast(g);
        end
        g=[];
    end
    gap=((bestUB-leastLB)/bestUB)*100;
    if gap<tol2 && flag==0
        flag=1;
        time_tol2=time_gurobi;
    end
    if time_gurobi>time_exit
        flag=1;
        time_tol4=time_gurobi;
        if gap> tol2+sqrt(eps)
            time_tol2=time_gurobi;
        end
        return
    end
    if gap < tol4
        time_tol4=time_gurobi;
        return
    end

    if abs(((UB-s(3))/UB)*100) > tol4
        % BRANCHING
        [mid1,mid2]=sbb_branch(zeta_Rj,s(1),s(2));
        idx =[s(1); mid1;mid2; s(2)];
        for i=1:3
            [fl_opt,z_supp,time_yalmip,time_gurobi]=optsupport_cg(alpha,Gamma,D,C,l_in,...
                E,N,B,K,idx(i),idx(i+1),G,qhat,time_s,set_notintrdictble,diag_capnotintrdct,set_intrdictble,time_yalmip,time_gurobi,time_exit,l1b,cplx);
            [u,LB,zeta,time_yalmip,time_gurobi]=opt_lowerbound(fl_opt,size(z_supp,2),K,...
                idx(i),idx(i+1), alpha, Gamma,qhat,time_yalmip,time_gurobi,l1b,cplx);
            if time_gurobi>time_exit
                flag=1;
                time_tol4=time_gurobi;
                if gap(n_iter)> tol2+sqrt(eps)
                    time_tol2=time_gurobi;
                end
                return
            end
            que.addFirst([idx(i); idx(i+1); LB;zeta; size(z_supp,2);size(fl_opt,1);u(:);fl_opt(:);z_supp(:)]);
        end
    else
        que.addFirst(s);
    end
end
end
