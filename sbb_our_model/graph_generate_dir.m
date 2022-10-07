function [G,dir1] = graph_generate_dir(p,q)
% generates network with arcs between column pointing toward sink and
% within the same column, the direction is random.

% p*q+2 is the total number of nodes in the network including s and t
% s: stores the nodes in the network
% n_instance is the index of the instance
% s_dir is the seed for the direction
s=zeros(p*q+2,1);
for i=1:p*q+2       
    s(i)=i;
end



% generate a digraph
G=digraph(s,s); 

% generate arcs from source to its neighbours
for i=2:p+1   
    G=addedge(G,1,s(i));
end

% generate arcs from neighbours of sink to itself
for i=p*(q-1)+2:p*q+1   
    G=addedge(G,s(i),p*q+2);
end

% generate arc between columns
for i=1:q-1              
    for j=p*(i-1)+2:p*i+1
        G=addedge(G,s(j),s(j+p)); 
    end
end

dir1=rand(1,q*(p-1))';
index=1;
for i=1:q                
    for j=p*(i-1)+2:p*i
        if dir1(index)>0.5
            G=addedge(G,s(j),s(j+1));
        else
            G=addedge(G,s(j+1),s(j));
        end
        index=index+1;
    end
end

%remove self-loops
G = rmedge(G, 1:numnodes(G), 1:numnodes(G));
