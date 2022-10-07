function [fl_opt, supp_full, supp_rem]  = compute_f_l_k(cap,set_non_rem,set_rem,B,E,G,N,K)
    size_set_non_rem=size(set_non_rem,1);
    iter=1;
    edgelist=zeros(B,nchoosek(E-size_set_non_rem,B));
    for it1=1:E-size_set_non_rem
        if B==1
            edgelist(:,iter)=it1;
            iter=iter+1;
        else
            if B==2
                for it2=it1+1:E-size_set_non_rem
                    %for it3=it2+1:E-size_set_non_rem
                    edgelist(:,iter)=[it1;it2];
                    iter=iter+1;
                    %end
                end
            else 
                for it2=it1+1:E-size_set_non_rem
                    for it3=it2+1:E-size_set_non_rem
                        edgelist(:,iter)=[it1;it2;it3];
                        iter=iter+1;
                    end
                end
            end
        end
    end
    parfor iter=1:nchoosek(E-size_set_non_rem,B)
        [fl_opt(iter,:), supp_full(:,iter), supp_rem(:,iter)] = flows_supp(cap,edgelist, E, set_non_rem, iter, size_set_non_rem, G, N, K)
    end
end
function [fl_opt1, supp_full1, supp_rem1] = flows_supp(cap,edgelist, E, set_non_rem, iter, size_set_non_rem, G, N, K)
    %      iter
    a= edgelist(:,iter);
        l_in=zeros(E-size_set_non_rem,1);
        l_in(a)=1;
        supp_rem1 = l_in;
        inc=1;
        z_in=zeros(E,1);      
        for e=1:E
            if ~ismember(e,set_non_rem)
                z_in(e)=l_in(inc);
                inc=inc+1;
            end
        end
        supp_full1= z_in;
        fl_opt1=flows_test(G, z_in, cap, N,K);
end