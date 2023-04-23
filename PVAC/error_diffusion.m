function [s ,query]=error_diffusion(G,m,t_th)


%Note: In this code the diffusion starts from just one node


% Initial nodes from which the error starts to diffuse


 
%% Calculating the normalization wp for diffusing the error and path
[~,Dlow]=sort(rand(G.N,1));
query=[Dlow(1);zeros(G.N-m,1)];
visited=zeros(G.N,1);
w_denpminator=zeros(G.N,1);
Nodes=1:1:G.N;
propagated = [];
for aux_v=1:1:G.N
    v=query(aux_v);
    visited(v)=1;
    non_visited_nodes=setdiff(Nodes,find(visited));
    propagated=setdiff(propagated,find(visited));
    %finding the neighbors of v
    flag=0;
    [~,Neigh_v,~]=find(G.W(v,:));
    Neigh_v=intersect(Neigh_v,non_visited_nodes); %non visited neighbors
    if ~isempty(Neigh_v)
        [~,idx]=sort(G.W(v,Neigh_v));
        Neigh_v=Neigh_v(idx);
        for ii=1:1:max(size(Neigh_v))
            w_denpminator(v)=w_denpminator(v)+G.W(v,Neigh_v(ii));
            propagated = [propagated Neigh_v(ii)];
            if flag==0
                query(sum(query>0)+1)=Neigh_v(ii);
                flag=1;
            end      
        end 
    else
        if aux_v<G.N && flag==0 && ~isempty(propagated)
            query(sum(query>0)+1)=propagated(1);
        elseif aux_v<G.N && flag==0
            query(sum(query>0)+1)=non_visited_nodes(1);
            
        end
    end
end    

s=zeros(G.N,1);
aux_step=1;
%% Diffusing the error in all the nodes in the graph
while sum(s)==0  % To guarantee that I am taking at least some samples. 

x=((m*aux_step)/G.N)*ones(G.N,1);
visited=zeros(G.N,1);
e=zeros(G.N,1);
for aux_v=1:1:G.N
    v=query(aux_v);
    visited(v)=1;
    u=x(v)-e(v);
    if u>t_th && sum(s(G.W(v,:)>0))==0
        s(v)=1;
    else
        s(v)=0;
    end     
    ep=s(v)-u;
    % finding the neighbors of v
    [~,Neigh_v,~]=find(G.W(v,:));
    
    for ii=1:1:max(size(Neigh_v))
        if visited(Neigh_v(ii))~=1
            e(Neigh_v(ii))=e(Neigh_v(ii))+ep*G.W(v,Neigh_v(ii))/w_denpminator(v);   
        else
        end      
    end     
    
end   

aux_step=aux_step+1;
end

%% Checking if the number of ones in s is equal to m
if sum(s)<m    
      s_discard=s;
      %discarding the sampling points in s and its neighbors
      [aux_var]=find(s);
      for rr=1:1:max(size(aux_var))
          [~,Neigh_v,~]=find(G.W(aux_var(rr),:));
          s_discard(Neigh_v)=1;
      end 
      [aux_var]=find(not(s_discard));
           if min(size(aux_var))>0
              aux_var=aux_var(randperm(max(size(aux_var))));% random permutation of vector elements 
                 for rr=1:1:min((m-sum(s)),max(size(aux_var)))
                     if sum(s(G.W(aux_var(rr),:)>0))==0
                        s(aux_var(rr))=1; 
                     end
                 end
                 if sum(s)<m
                    [aux_index,~,~]=find(not(s));
                    aux_index=aux_index(randperm(length(aux_index)));
                    s(aux_index(1:(m-sum(s))))=1;  
                 else
                 end    
           else
                 aux_var=find(not(s));
                 aux_var=aux_var(randperm(length(aux_var)));
                 for rr=1:1:m-sum(s)
                      s(aux_var(rr))=1; 
                 end
           end      
elseif sum(s)>m
      [aux_var]=find(s);
      % Permuting randomly the elements
      aux_var=aux_var(randperm(max(size(aux_var))));
      
      aux_var=aux_var(1:1:m);
      s=zeros(size(s));
      s(aux_var)=1; 
end    


s=s(:);
if sum(s)==m
   disp('Number of sampling nodes in error diffusion is m, OK');
   [sum(s), m]
else
   disp('Number of sampling nodes in error diffusion is different from m')   
end    
    



end

