[samplig_nodes]=dot_diffusion_partition(G,m,partitions)
% Input
    % G graphs (G.N number of nodes) (G.W weight matrix)
    % m numnber of samples
    % partitions: vector length N with entries representing the nodes partition
% Ouput:
    % set of sampling nodes

    N =  G.N
    s = zeros(N);
    selected_nodes = {}
    for a=1:length(unique(partitons))
        nodes=find(map==a);   
        G_aux.W=G.W(nodes,nodes);  
        G_aux.N=length(nodes);
        m_temp=round(m*sum(G_aux.W(:))/sum(G.W(:)));
        [S_VC] =  dot_diffusion(G_aux,m_temp,0.5);
        selected_nodes{a}=pos(find(S_VC==1));
    end
    S = zeros(G.N,1);
    for ii=1:length(unique(map))
        S(selected_nodes{ii})=1;
    end
    samplig_nodes = find(S)
end