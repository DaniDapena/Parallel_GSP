function [S_0,S_in] = my_vc_v00(w,m,IP,a,Nodes_b)

%Input
    % w weight matrix
    % density_pattern density of the patter view as the ration of selected
    % IP initial patter (mask)
    % U control variable to use density


% Output
    % S final pattern vector quith value of 1 in the selected node
    % otherwise
    % S _in initial random pattern
    %% Inicial
    N=max(size(w));
    density_pattern=m/N;
    Ind_A=-1;
    Ind_B=-1;
    S_0=IP;
    max_r=5*N;
    c=zeros(N,1);
    optCost = 1e6;
    boundary =zeros(N,1);
    boundary(Nodes_b) = 1;
    %% Density mapping
    for i=1:N
        d(i)=(1/mean(w(i,:)));
    end
    
    

    if a~=0
        for i=1:N
            for j=1:N
                p=max(d(i),d(j));
                w(i,j)=w(i,j)*p^(a);
            end
        end
    end
    
    %% pricinpal path length
    G=graph(w);
    Geodesics=distances(G);
    
    WD = Geodesics;
    num_points=3000;
    a=0; b=max(nonzeros(WD(:)));
    step=(b-a)/num_points;
    dist=step*(1:1:num_points);
 

    ht=0;
    for rr=1:1:N
        v=WD(rr,:);
        v=nonzeros(v);
        for ii=1:1:max(size(dist))
            h(ii) = sum(v<dist(ii))+1;    
        end  
        ht=ht+h; clear h
    end
    ht=ht/N;
    aux_var=1./ht;
    [~,xlim,~]=find(aux_var<=0.8);
    aux_var=aux_var(xlim);
    dist=dist(xlim);
    
    V_1=abs(density_pattern-(aux_var));
    Den=find(V_1==min(V_1));
    lambda=dist(Den(1));
   
    sigma=lambda^2/log(10); % parameter to adjust the Kernel
    K_0 = exp(-(Geodesics).^2/(sigma)); % Mapping of the geodesig distace using a Kernel
%% Initial patterns 
   

    initial_values = randperm(N);


    Number_of_samples = m-sum(IP); % Number of selected nodes in the pattern
   
    [I] =randperm(N);
    I = setdiff(I,Nodes_b);
    I = I(randperm(length(I)));
    
    S_0(I(1:Number_of_samples))=1;
    S_in = S_0;
%% Sampling
    
    for r=1:max_r
        
    	c(find( S_0)) = sum(K_0(find(S_0),find( S_0)));
        c(find(~S_0)) = sum(K_0(find(S_0),find(~S_0))) - 1000;
        pos  = find(boundary==0);
        position_max = pos(find(c(pos) == max(c(pos))));
        position_max = position_max(1);
        position_min = pos(find(c(pos) == min(c(pos))));
        position_min = position_min(1);
        
        if (sum(c(find(S_0))) < optCost)
            optCost = sum(c(find(S_0)));
        else
            S_0(position_max) = 1;
            S_0(position_min ) = 0;
            break;
        end
%         position_max=position_max(randperm(length( position_max),1));
%         position_min = position_min(randperm(length(position_min),1));

        S_0(position_max) = 0;
        S_0(position_min ) = 1;

        if Ind_A ==  position_max && Ind_B == position_min
            break
        else
            Ind_A = position_min;
            Ind_B = position_max;
        end
  
    end
    
%     if sum(S)==round(N*density_pattern)
%        disp('The number of samples is OK')
%    else
%        if sum(S)>round(N*density_pattern)
%            disp('Error: The number of samples is greater than m')
%            [sum(S), round(N*density_pattern)]
%        else
%            disp('Error: The number of samples is smaller than m')
%            [sum(S), round(N*density_pattern)]
%        end
%    end
end

