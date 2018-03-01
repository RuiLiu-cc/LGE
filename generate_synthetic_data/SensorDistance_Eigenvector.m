function [Utop,Lap,W] = SensorDistance_Eigenvector(Num_nodes,prob,ampliphy,strange,Num_PC)
%This function is generate the adjacency matrix of a graph and based on the 
%graph to calculate the eigenvectors as the basic components for the data.
%Each pair of nodes in the graph has a probability p to connect together. 
%The weight of edge is uniformly random drew from 0 to 1.
%
% input: 
% num: number of nodes in the graph (number of feature)
% prob: probability of a node connecting to others
% ampliphy: increase amplitude of the edge weight
% strange: if strange = 1, then there will be 20% disturbance in the graph;
% if strange ==0, then there is no disturbance in the graph
% NoCp: number of components selected as the basic components.
% 
% output:
% Utop: selected eigenvectors
% L1: Laplacian matrix of the graph.
% W: adjacency matrix of the graph
% 
% Note: the adjacency matrix (W), normalized Laplacian matrix (L1) and Utop 
% will be saved in the "generate_synthetic_data" folder.
%
%==========================================================================
%
% Author: Liu Rui, SUTD, 4 May 2017
%
%==========================================================================

%% generate the graph

    Distance=zeros(Num_nodes);
    D = zeros(Num_nodes);
    while length(find(sum(D)==0))>0
        %number of pairs
        npairs=(1+Num_nodes)*Num_nodes/2;
        D = rand(1,npairs);
        for i=1:npairs
            R=binornd(1,prob);
            D(1,i)=ampliphy*D(1,i)*R;
        end

        if strange==1
            %choose some weight to increase
            Index= randperm(npairs,ceil(npairs/5));
            D (Index)= D(Index)*10;
        end

        t=1;
        for x=1:Num_nodes
            for y=x+1:Num_nodes
            Distance(x,y)=D(1,t);
            Distance(y,x)=D(1,t);
            t=t+1;
            end
        end
        figure(1);imagesc(Distance);
    
        %% calculate eigenvectors
        W=Distance;
        %creat matrix D
        D=zeros(Num_nodes);

        for i=1:Num_nodes
            for j=1:Num_nodes
                D(i,i)=D(i,i)+W(i,j);
            end
            D(i,i)=D(i,i)-W(i,i);
        end
    end
    
    %calculate Laplacian matrix L
    L=zeros(Num_nodes);
    L=D-W;

    %calculate normalized graph Laplacian L1
    D1=zeros(Num_nodes);
    for i=1: Num_nodes
        D1(i,i)=1/sqrt(D(i,i));
    end
    Lap=D1*L*D1;

    %calculate eigenvectors U and eigenvalues lamda
    [U,lamda1]=eig(Lap); 
    lamda2 = sum(lamda1,1);
    [lamda,index] = sort(lamda2,'ascend');

    %choose top eigenvectors to generate data and save
    Utop= U(:,index(1:Num_PC));
    save (sprintf('generate_synthetic_data/Graph_%sEigenvectors_NumNodes=%s_prob=%s_amplify=%s.mat',num2str(Num_PC),num2str(Num_nodes),num2str(prob),num2str(ampliphy)),'U','Utop','Lap','W');
    
end