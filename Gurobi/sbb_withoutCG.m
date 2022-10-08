function [time_tol4, bestUB, ustar,time_tol2,time_yalmip, flag]=sbb_withoutCG(fl_opt,K,Gamma,alpha,...
    ~,flag,qhat,n_iter_CD,zeta_lb,zeta_ub,...
    tol2,tol4,time_exit,set_notintrdictble,diag_capnotintrdct,~,time_yalmip,time,l1b, cplx)
    L=size(fl_opt,1);
    fh=0;
    time_s=cputime;
    import java.util.ArrayDeque
    que = ArrayDeque();
    bestUB=1e6;
    z_lb=zeta_lb;
    [u,LB,zeta,time_yalmip,time]=opt_lowerbound(fl_opt,L,K,z_lb,zeta_ub, alpha, Gamma,qhat,...
        time_yalmip,time,l1b, cplx);
    vec=[ z_lb; zeta_ub; LB; zeta;u(:)];
    que.addFirst(vec);
    n_iter=1;
    while (~isEmpty(que))
        leastLB=1e6;
        s_q=size(que);
        quedum = ArrayDeque();
        for j=1:s_q 
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
        [UB, zeta_Rj,uf,time_yalmip,time]=coord_descent(fl_opt,...
            s(5:end),Gamma,alpha,K,qhat,L,s(1),s(2),...
            n_iter_CD,time_s,set_notintrdictble,diag_capnotintrdct,time_yalmip,time,l1b, cplx);
        if UB<bestUB
            bestUB=UB;
            ustar=uf;
        end
        % PRUNING
        s_que=size(que);
        for j=1:s_que
            g=que.removeFirst();
            if g(3)< bestUB || abs(g(3)-bestUB)<sqrt(eps)
                que.addLast(g);
            end
        end
        gap=((bestUB-leastLB)/bestUB)*100;
        if gap<tol2 && fh==0
            fh=1;
            time_tol2=time;
        end
        if time>time_exit
            flag=1;
            time_tol4=time;
            if gap> tol2+sqrt(eps)
                time_tol2=time;
            end
            return
        end
        if gap < tol4
            time_tol4=time;
            return
        end

        if abs(((UB-s(3))/UB)*100) > tol4
            % BRANCHING
            [mid1,mid2]=sbb_branch(zeta_Rj,s(1),s(2));
            idx =[s(1); mid1;mid2; s(2)];
            for i=1:3
                [u,LB,zeta,time_yalmip,time]=opt_lowerbound(fl_opt,L,K,...
                idx(i),idx(i+1), alpha, Gamma,qhat,time_yalmip,time,l1b, cplx);
                que.addFirst([idx(i); mid1(i+1); LB;zeta; u(:)]);
                if time>time_exit
                    flag=1;
                    time_tol4=time;
                    if gap> tol2+sqrt(eps)
                        time_tol2=time;
                    end
                    return
                end
            end
        else
            que.addFirst(s);
        end
        n_iter=n_iter+1;
    end
end

