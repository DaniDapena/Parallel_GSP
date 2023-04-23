function sampling_nodes =  pVAC(W,n_samples, partitions)
    %  OUT: sampling_nodes: set of sampling nodes
    % IN:
    %  W weight matrix
    %  n_sampels: number of samples
    %  partitions: vector length N with entries representing the nodes partition

    
% Nodes
    N  = max(size(W))
%% Find the nodes on the border
    [i,j] = find(triu(W));
    aux_border = find(map(i) ~= map(j));
    border = unique([i(aux_border) j(aux_border)]);
    n_border = length(border);

%% Border Graphs
    W_border=W(border,border);
%% number of graphs 
    bins = conncomp(graph(W_p));
    unique_bins = unique(bins);

%% Initial pattern
    S_initial = zeros(N,1);
    selected_nodes={};
    S = zeros(N,1);
%% create pattern in each subgraph in the border]
    for k=1:length(unique(bins))
        Nodes = find(bins==unique_bins(k));
        Nodes = unique(map(border(Nodes)));
        Gbin.N = length(Nodes);
        Gbin.W = W(Nodes,Nodes);
        IP_bin = zeros(Gbin.N,1); % initial pattern in border
        m_temp = full(round(n_samples*sum(sum(Gbin.W))/sum(W(:))));
        Nodes_b = [];
        if m_temp>0
            [pattern_path_bin ,query] =e rror_diffusion(Gbin,m_temp,0.5)
            S_initial(Nodes) = pattern_path_bin;
        end
    end
    %% Sample path graph
    for partition=1:max(partitions)
        Nodes= find(partitions=partition);
        Nodes_border = intersect(Nodes,border);
        W_partition = W(Nodes,Nodes);  
        IP = S_initial(Nodes);
        N_aux = length(Nodes);
        m_temp=full(round(n_samples*sum(sum(W(Nodes,Nodes)))/sum(G.W(:))));
        S_VC_partition =  my_vc_v00(W_partition,m_temp,IP,0,Nodes_border);
        selected_nodes{partition}=Nodes(find(S_VC_partition==1));
    end
    for ii=1:length(unique(map))
        S(selected_nodes{ii})=1;
    end
    sampling_nodes = find(S)
end
