function x_rec = pbmr(x,sampling_nodes,L, epsilon, q, partitions)
    % x: noisy singnal vector length N with entries representing the nodes signal value
    % sampling_nodes: set of sampling nodes
    % L lamplacian matrix
    % epsilon: pertubation
    % q: laplacian power
    % partitions: vector length N with entries representing the nodes partition

    N = length(x) % number of nodes
    x_rec= zeros(N,1);
    c= zeros(N,1);

    for partition = 1:max(partitions)
        Nodes = find(partitions == partition);
        power = 1;
        Neighbor_nodes = Nodes % for initial serach
        while power<=q % incread to the  -hop
            Neighbor_nodes_aux  = []
            for node = 1:length(Neighbor_nodes)  % look for nodes neighbirs
                aux_nodes =  find(G.W(Nodes(node),:));
                aux_nodes = setdiff(aux_nodes,Nodes) % nodes that are not in the the set A yet
                Neighbor_nodes_aux  = union(Neighbor_nodes_aux , aux_nodes); % update aux variable
            end
            Neighbor_nodes = Neighbor_nodes_aux % update looking for neighbors
            Nodes  = union(Nodes, Neighbor_nodes); $ update A
            power = power+1;
        end
        c(Nodes) = c(Nodes) +1; % count how many times the nodes appears ina set A

        % Parameters and matris of the subgraph induced by A
        W_A = G.W(Nodes,Nodes);
        L_A = diag(sum(W_A)) - W_A;
        N_A = length(Nodes);
        L_A= (L_A^power_laplacian)+epsilon_set(i)*eye(N_A); 
        inverted_Laplacian_A = inv(L_A);

        % Sampling matrix m
        sampling_anis = sampling_patterns_anis(j,Nodes);
        M = zeros(length(sampling_nodes),N_A);
        ind_M = 1;
        for n=1:length(sampling_nodes)
            M(ind_M,sampling_nodes(n)) = 1;
            ind_M = ind_M + 1;
        end

        % Recovery Matrix
        recovery_matrix = inverted_Laplacian_A*M'*inv(M*inverted_Laplacian_A*M');

        % Signal in the partition
        x_A = x(Nodes);

        %% Sampling 
        x_A_sampled = M*x_A;

        % Recovery:
        x_rec_A = recovery_matrix*x_A_sampled;

        % Update signal
        x_rec(Nodes) =  x_rec(Nodes) + x_rec_A;
                
                
    end
    % divide signal
    x_rec =  x_rec./c;

end
