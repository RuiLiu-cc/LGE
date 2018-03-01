function [update_Lr,update_M,rank, Vn] = update_LowRank_UpdateLfirst(Graph,Data, para)
% This function is to update low rank matrix L and noise matrix M with a
% fixed graph. This step is based on RPCAG with ADMM method.
%==========================================================================
% Input:
% Graph: the laplacian matrix of input graph which is fixed.
% Data: noisy input data as X in the paper. It is (pxn).
% para.sigma: control the noise sparsity.
% para.gamma: control the smoothness of the signal on the graph.
% r1: The lagrangian penalty term
% r2: The second lagrangian penalty term.
%
% Output:
% update_Lr: estimated low rank matrix L.
% update_M: estimated noisy matrix M.
% rank: the rank of L.
% Vn: selected pricinple components.
%==========================================================================
%
% Author: Liu Rui, SUTD, 27, Feb, 2018
%
%==========================================================================
%% initialization
sigma = para.sigma;
gamma = para.gamma;
r1 = para.r1;
r2 = para.r2;
tau1 = 2/(r1 + r2);
[Num_feature,Num_sample]=size(Data);

% compute the inverse of the matrix for the gradient solution using the
% graph regularization. 
Temp = gamma*Graph + r2*eye(Num_feature);
Tempinv = inv(Temp);
% initialize all the matrices
J = []; 
K = rand(Num_feature,Num_sample); % initial L with randomness
L = K;
M = Data - L;
Z1 = Data - L - M;
Z2 = K - L;
Z1_old = Z1;
Z2_old = Z2;
K_old = K;
M_old = M;
L_old = L;
% intialize matrices to evaluate the termination of the algorithm.
norm_Z1 = [100];
norm_Z2 = [100];
L_NuclearNorm = [100 200];
M_l1Norm = [100 200];
Trace_Norm = [100 200];


%% Update Lr, M and K
for i = 1: 1500 
    % Update Lr. 
    E1 = Data - M + Z1/r1;
    E2 = K + Z2/r2;
    J = (r1*E1 + r2*E2)/(r1+r2);

    if isnan(J(1,1))
        rank = 1;
        Vn= 0;
        return
    elseif abs(J(1,1))>1e+100
        rank = 1;
        Vn = 0;
        return
    else
        [Usvd, Sigmasvd, Vsvd] = svdecon(J); 
    end
    
    Ssvd_thres = softThresh(Sigmasvd,tau1);
    L = Usvd*Ssvd_thres*Vsvd';
    
    % Update M
    J = Data - L + Z1/r1;
    M = softThresh(J,sigma/r1);
    
    % Update K
    K = r2*Tempinv*(L-Z2/r2);
    
    % Update the dual variables / lagrangian parameters.
    Z1 = Z1 + r1*( Data - L - M );
    Z2 = Z2 + r2*(K - L);
    
    % Relative difference
    if(i > 1)
        L_NuclearNorm = [L_NuclearNorm sum((diag(Ssvd_thres)))];
        M_l1Norm = [M_l1Norm sigma*sum(sum(abs(M)))];
        Trace_Norm = [Trace_Norm gamma*trace(L'*Graph*L)];
        norm_Z1 = [norm_Z1 sum(sum((Z1-Z1_old).^2))];
        norm_Z2 = [norm_Z2 sum(sum((Z2-Z2_old).^2))];
    end
    
    % stop criteria
    if abs(L_NuclearNorm(end)-L_NuclearNorm(end-1))/L_NuclearNorm(end) <0.001 && abs(M_l1Norm(end)-M_l1Norm(end-1))/M_l1Norm(end)<0.001...
        && abs(Trace_Norm(end)-Trace_Norm(end-1))/Trace_Norm(end)< 0.001 && norm_Z1(end)< 0.0005 && norm_Z2(end)< 0.0005
        break 
    end
    Z1_old = Z1;
    Z2_old = Z2;
    K_old = K;
    M_old = M;
    L_old = L;

end 

% assign the results to the output.
update_Lr = L;
update_M = M;

% determine the rank of the soltuion after the iteration has terminated.
indicator = find(diag(Ssvd_thres)>0);
rank = length(indicator);

% check if the matrix L is full-rank
if(max(indicator) == Num_feature)
    warning('Your low-rank matrix is full rank. The parameter sigma is quite big.')
    Y = Usvd(:,1:max(indicator));               
    Vn = Vsvd(:,1:size(Y,2));
else
    if(length(indicator) > 1)
        % if the matrix is not full rank or rank 1.
        Y = Usvd(:,1:max(indicator));               
        Vn = Vsvd(:,1:size(Y,2));        
    else
        % if the matrix is rank 1
        warning('Your low-rank matrix is rank 1. The parameter sigma is quite small.')
        Y = Usvd(:,1:max(indicator));               
        Vn = Vsvd(:,1:size(Y,2));
    end 
end


end