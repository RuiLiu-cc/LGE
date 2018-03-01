function [Graph] = update_Graph(Lr,para,Num_dim)
% This function is to update graph G with a fixed low rank matrix. It is 
% revised from Dong Xiaowen's code for learning Laplacian Matrix in Smooth
% Graph Signal Representation Paper:http://web.media.mit.edu/~xdong/pub.html
%
% input:
% Lr: low rank matrix 
% para.sigma: The penalty term for graph regularization.
% para.beta: control the off-diagonal elements distribution in the graph
% matrix.
% Num_dim: number of features of data
%
% output:
% Graph: update laplacian matrix of the graph.
%==========================================================================
%
% Author: Liu Rui, SUTD, 27 Feb 2018
%
%==========================================================================

%% initialization
alpha = para.sigma;
beta = para.beta;

%% calculate Laplacian constrains
[A1,b1,A2,b2,mat_obj] = laplacian_constraint_vech(Num_dim);
p = vec(Lr*Lr')';

%% optimization with cvx toolbox
cvx_begin quiet %turn off the information showing
% cvx_solver gurobi
variable vechG(Num_dim*(Num_dim+1)/2,1)
minimize alpha*p*mat_obj*vechG + beta*sum_square_abs(mat_obj*vechG)
subject to
    A1*vechG == b1
    A2*vechG <= b2
cvx_end

%% convert from vector form to matrix form
Graph = reshape(mat_obj*vechG,Num_dim,Num_dim);

end