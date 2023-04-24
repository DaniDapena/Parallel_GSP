function [s]=dot_diffusion_partition(G,m,t_th)

    % G graphs (G.N number of nodes) (G.W weight matrix)
    % m numebr of samples
    % t_th thresholf

Ad=G.W;
Ad(G.W>0)=1;
D=sum(Ad);
[ ~, idx_D]=sort(D,'descend');
[Dlow]=idx_D(1);
idx_D(1)=[];
W1=zeros(G.N,G.N);
%% Calculating the normalization wp for diffusing the error

qq=[Dlow;zeros(G.N-1,1)];
visited=zeros(G.N,1);
wp=zeros(G.N,1);
Nodes=1:1:G.N;
SP=[];
for aux_v=1:1:G.N
    v=qq(aux_v);
    visited(v)=1;
    non_visited_nodes=setdiff(Nodes,find(visited));
    %finding the neighbors of v
    [~,Neigh_v,~]=find(G.W(v,:));
    Neigh_v=intersect(Neigh_v,non_visited_nodes);
    if ~isempty(Neigh_v)
        [~,idx]=sort(D(Neigh_v),'descend');
        Neigh_v=Neigh_v(idx);
%          sigma=mean(G.W(v,Neigh_v));
        for ii=1:1:max(size(Neigh_v))
%           G.W(v, Neigh_v(ii))=exp(-G.W(v,Neigh_v(ii))^2/sigma^2);
             G.W(v, Neigh_v(ii))=G.W(v,Neigh_v(ii));
            wp(v)=wp(v)+G.W(v,Neigh_v(ii));
  
        end 
    end
    if aux_v<G.N  %% anternate high a low degree
        if mod(aux_v,2)~=0
            qq(sum(qq>0)+1)=idx_D(end);
            idx_D(end)=[];
            if G.W(qq(sum(qq>0)),qq(sum(qq>0)-1))>0
                W1(qq(sum(qq>0)),qq(sum(qq>0)-1))=1;
            end
        else
            qq(sum(qq>0)+1)=idx_D(1);
            idx_D(1)=[];
            neig= find(G.W(qq(sum(qq>0)),qq(qq>0)));
            if G.W(qq(sum(qq>0)),qq(sum(qq>0)-1))>0
                W1(qq(sum(qq>0)),qq(sum(qq>0)-1))=1;
            end
        end
    end
end    

s=zeros(G.N,1);
aux_step=1;
%% Diffusing the error in all the nodes in the graph
flag=0;
while sum(s)==0  % To guarantee that I am taking at least some samples. 

x=((m*aux_step)/G.N)*ones(G.N,1);
visited=zeros(G.N,1);
e=zeros(G.N,1);
for aux_v=1:1:G.N
    v=qq(aux_v);
    visited(v)=1;
    u=x(v)-e(v);
    [~,Neigh_v,~]=find(G.W(v,:));
    UU(v)=u;
    if u>t_th && sum(s(Neigh_v))==0
        s(v)=1;
    else
        s(v)=0;
    end     
    ep=s(v)-u;
    
    E(v)=ep;
    % finding the neighbors of v
    [~,Neigh_v,~]=find(G.W(v,:));
    
    for ii=1:1:max(size(Neigh_v))
        if visited(Neigh_v(ii))~=1
            e(Neigh_v(ii))=e(Neigh_v(ii))+ep*G.W(v,Neigh_v(ii))/wp(v);  
        else
        end      
    end     

end   
aux_step=aux_step+1;
end

%Checking if the number of ones in s is equal to m
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
      %Permuting randomly the elements
      aux_var=aux_var(randperm(max(size(aux_var))));
      
      aux_var=aux_var(1:1:m);
      s=zeros(size(s));
      s(aux_var)=1; 
end    


s=s(:);
% if sum(s)==m
%    disp('Number of sampling nodes in error diffusion is m, OK');
%    [sum(s), m]
% else
%    disp('Number of sampling nodes in error diffusion is different from m')   
% end    
    
% 
% 

end

