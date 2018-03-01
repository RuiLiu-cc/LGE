function [Lr,rank_Lr,OptG,Vn] = LGE(initG, Data, Lr_0, Lap, para)
% This function is the implementation of Low rank matrix and graph estimation (LGE)
% algorithm. 
%==========================================================================
% Input:
% initG: intial graph as \phi_f (Laplacian matrix of graph) in the paper.
% Data: (pxn) noisy input data as X in the paper. p is the number of features
% and n is the number of samples.
% Lr_0: groudtruth low rank matrix.
% Lap: groudtruth Laplacian matrix of the graph.
% para
% para.sigma: control the noise sparsity.
% para.gamma: control the smoothness of the signal on the graph.
% para.beta: control the off-diagonal elements distribution in the graph
% matrix.
% para.max_Iter: maximum iteration number.
% para.r1: the parameter for the ADMM in step 1 (RPCAG).
% para.r2: the parameter for the ADMM in step 1 (RPCAG).
%
% Output:
% Lr: estimated low-rank matrix
% rank_Lr: the rank of Lr
% OptG: estimated graph (Laplacian matrix of graph)
% Vn: the selected principle compoents based on the rank k.
%==========================================================================
%
% Author: Liu Rui, SUTD, 27 Feb 2018
% 
% For more details, please refer to the paper: Joint Estimation of Low-Rank 
% Components and Connectivity Graph in High-Dimensional Graph Signals: 
% Application to Brain Imaging: https://arxiv.org/abs/1801.02303
%
%==========================================================================
%% initialization
[Num_dim,Num_sample]= size(Data);
update_Lr = ones(Num_dim,Num_sample);
update_M = ones(Num_dim,Num_sample);
update_G = initG;
obj_value = zeros(2,para.max_Iter);

%% start iterations until converge
disp('Algorithm starts');
disp('Iter|   Lr    |    M    |tr(LrGLr) |   G      | obj_value | Lr error  | G error  |Diff edge ratio(G)| rank |');

for i = 1: para.max_Iter
    %% step1: Given graph, update the low rank estimation
    [update_Lr,update_M,rank_Lr,Vn] = update_LowRank_UpdateLfirst(update_G, Data, para);
    % if update_Lr is NaN or Vn is 0, parameters need to be changed.
    if isnan(update_Lr(1,1)) || isinf(update_Lr(1,1))
        update_Lr = 0;
        sprintf('This setting does not work. \n');
        break
    end
    if  Vn == 0
        update_Lr = 0;
        sprintf('This setting does not work. \n');
        break
    end
    % print out the value of each term in the objective function
    obj_value(1,i) = print_out_each_part(update_Lr,update_M,update_G,Lr_0,Lap,para,rank_Lr,i);
    % plot Lr and graph after this step
    plot_Lr_Graph_eachStep(update_Lr,update_G);
      
    %% step2: Given low rank estimation, update graph
    [update_G] = update_Graph(update_Lr,para, Num_dim);
    % make the elements whose absolute value < 1e-3 in the G to be zero
    update_G(abs(update_G)<1e-3) = 0;
    % print out the value of each term in the objective function
    obj_value(2,i) = print_out_each_part(update_Lr,update_M,update_G,Lr_0,Lap,para,rank_Lr,i);
    % plot Lr and graph after this step
    plot_Lr_Graph_eachStep(update_Lr,update_G);
       
    %% stop criteria
    if i>=2 && abs(obj_value(2,i)-obj_value(2,i-1))/abs(obj_value(2,i-1))< 10^(-3)
        break
    end
    %update Data, use the updated Lr as the new X in the X = L+ M
    Data = update_Lr; 
    
end

OptG = update_G;
Lr = update_Lr;

end