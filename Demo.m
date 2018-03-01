% Small DEMO: test the LGE with some synthetic data
%==========================================================================
%
% Author: Liu Rui, SUTD, 27 Feb 2018
% 
% For more details, please refer to the paper: Joint Estimation of Low-Rank 
% Components and Connectivity Graph in High-Dimensional Graph Signals: 
% Application to Brain Imaging: https://arxiv.org/abs/1801.02303
%
%==========================================================================
clear; close all;
dbstop if error;
addpath('./generate_synthetic_data/');
addpath('./HalfVectorization/');
addpath('./related_function/');

%% initialization

% Initail the parameter for generating graph signal
Num_samples = 50; % number of samples
Num_nodes = 30; % number of features
prob = 0.3; % probability of a node connecting to others
ampliphy = 5; % amplitude of the edge weight
strange = 1; % if strange = 1, then there will be 20% disturbance in the graph, otherwise no disturbance.
Num_PC = 3; % number of components selected as the basic components
Para.mu = 0; % mean value of normal distribution for Y
Para.k = 0.1; % ratio of disturbed entries over whole dataset

% Choose the way to generate input data
% if New_data =1, then use the function to generate new data. 
% if New_data =2, then load the present data directly.
New_data = 2;

% Initial the parameter of LGE algorithm
para.sigma = 0.4; % control the noise sparsity
para.gamma = 1; % control the smoothness of the signal on the graph.
para.beta = 0.5; % control the off-diagonal elements distribution in the graph matrix.
para.max_Iter = 50; % maximum iteration number of the GLE
para.r1 = 0.1; % the parameter for the ADMM in step 1.
para.r2 = 0.1; % the parameter for the ADMM in step 1.

% Choose the way to generate initial graph
% if Initial_graph =1, then initial the graph from the data with gspbox
% function. Then we use k nearest neigberhood to build the graph: ParaG.k
% if Initial_graph =2, then initial the graph with random graph.
% if Initial_graph =3, then initial the graph with the groundtruth laplacian matrix.
Initial_graph = 1; 
ParaG.k =5; 

% Fix the randomness
rng(3);

%% generate graph and graph signal as groundtruth
if New_data == 1 
    [Data,Lap,Adj,Lr_0] = Generate_synthetic_data (Num_samples,Num_nodes,...
    prob,ampliphy,strange,Num_PC,Para);
    % plot the grounthruth low rank matrix, input data and groundtruth
    % graph.
    figure(1);
    subplot(2,3,1);imagesc(Data);title('Input noisy Data');
    subplot(2,3,2);imagesc(Lr_0);title('Low-rank matrix (Groundtruth)');
    subplot(2,3,3);imagesc(Lap);title('Graph (Groundtruth)');

elseif New_data == 2
    load Experiment2_Data_3Eigenvectors_NumNodes=30_NumSamples=50_k=0.1-1_Mode1.mat;
    Data = Data_all{1};
    load Graph_3Eigenvectors_NumNodes=30.mat;
    Adj = W;
    % plot the grounthruth low rank matrix, input data and groundtruth
    % graph.
    figure(1);
    subplot(2,3,1);imagesc(Data);title('Input noisy Data');
    subplot(2,3,2);imagesc(Lr_0);title('Low-rank matrix (Groundtruth)');
    subplot(2,3,3);imagesc(Lap);title('Graph (Groundtruth)');
end
%% calculate the initial graph for GLE
if Initial_graph ==1
    tempG = gsp_nn_graph(Data,ParaG); % graph generation function in gsp toolbox.
    initG = tempG.L;
    
elseif Initial_graph ==2
    tempG = rand (Num_nodes);
    mask = tempG>1-prob;
    initG = sparse(mask.*tempG);
    
elseif Initial_graph ==3
    initG = sparse (Lap);
    
end


%% estimate the low rank components and graph with proposed method LGE
[Lr,rank_Lr,OptG,Vn] = LGE(initG, Data, Lr_0, Lap, para);

%% Show the estimation results
% plot the estimation matrix
figure(1);
subplot(2,3,5);imagesc(Lr);title('Estimated low-rank matrix(Lr)');
subplot(2,3,6);imagesc(OptG);title('Estimated Graph');

% calculate estimation error
Lr_error = norm(Lr-Lr_0,'fro')/norm(Lr_0,'fro')
Graph_error = norm (OptG - Lap,'fro')/norm(Lap,'fro')
