function [Data,Lap,Adj,Lr_0] = Generate_synthetic_data (Num_samples,Num_nodes,prob,ampliphy,strange,Num_PC,Para)
% this function is to generate synthetic graph signal resident on a graph
% Get the eigenvectors of graph as principle components P (Num_nodes x Num_PC)
% Generate coefficient matrix Y (Num_samples x Num_PC) from N(mu,1/Num_samples) 
% distribution.
% Generate sparse noise matrix M from i.i.d Bernoulli distribution with k. (the
% elements that are not zero will have half +1 and half -1)
% Thus, the graph signal X (Num_nodes x Num_samples) is
%
%                 X = PY^T + M = Lr_0 + M
% 
%
% input:
% Num_samples: number of samples
% Num_nodes: number of nodes in the graph, or number of features.
% prob: probability of a node connecting to others (used for generate graph)
% ampliphy: increase amplitude of the edge weight(used for generate graph)
% strange: if strange = 1, then there will be 20% disturbance in the graph;
% if strange ==0, then there is no disturbance in the graph (used for generate graph)
% Num_PC: number of components selected as the basic components(used for generate eigenvector)
% Para.mu: mean value of normal distribution for Y.
% Para.k: ratio of disturbed entries over whole dataset. k =
% ||M||_0/(Num_samples*Num_nodes)
%
% output:
% Data: graph signal data 
% Lap: Lapalcian matrix as groundtruth of the graph
% Adj: Adjacency matrix of the graph
% Lr_0: low rank matrix as groundtruth
%==========================================================================
%
% Author: Liu Rui, SUTD, 4 May 2017
%
%==========================================================================
% generate the graph and its eigenvector as P (Num_nodes x Num_PC)
[P,Lap,Adj] = SensorDistance_Eigenvector(Num_nodes,prob,ampliphy,strange,Num_PC);

% generate coefficient matrix Y ~ N(mu,1/Num_nodes)
Y  = rand(Num_samples,Num_PC) * sqrt(1/ Num_nodes) + Para.mu;
% Y = normrnd(Para.mu,1/ Num_nodes,Num_samples,Num_PC);
Lr_0 = P*Y';
Lr_0 = normc(Lr_0); % normalize Lr_0 along each column



% generate sparse noise matrix M
M = zeros (Num_nodes, Num_samples);
Num_err = floor (Para.k * Num_samples * Num_nodes);
Idx = randperm (Num_samples * Num_nodes,Num_err); % generate the index of noisy elements in the dataset.
noise = randsrc (1,Num_err); % generate a (1,Num_err) vector with +1 or -1 equally
for i = 1:Num_err
    InX(i) = fix(Idx(i)/Num_samples)+1;
    InY(i) = mod(Idx(i),Num_samples);
    if InY(i) == 0
        InX(i)= InX(i) -1;
        InY(i) = Num_samples;
    end
    M(InX(i),InY(i)) = noise(i);
end

% get the graph signal
Data = Lr_0 + M;
figure(100);
subplot(1,3,1);imagesc(Lr_0);title('Lr_0');
subplot(1,3,2);imagesc(M);title('M');
subplot(1,3,3);imagesc(Data);title('Data');

save (sprintf('generate_synthetic_data/Data_%sEigenvectors_NumNodes=%s_NumSamples=%s_k=%s_Graph_prob=%s_amplify=%s.mat',...
    num2str(Num_PC),num2str(Num_nodes),num2str(Num_samples),num2str(Para.k),num2str(prob),num2str(ampliphy)),'Data','Lr_0');

end