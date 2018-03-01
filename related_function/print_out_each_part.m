function [obj_value,Lr_error,G_error,Diff_edge_error] = print_out_each_part(update_Lr,update_M,update_G,Lr_0,Lap,para,rank,iteration)
% this function is to print out the value of each term in the objective 
% function and gives the value of the whole objective function, Lr
% accuracy and graph accuracy.
%==========================================================================
%
% Author: Liu Rui, SUTD, 4 May 2017
%
%==========================================================================
% calculate the value 
nuclear_term = norm_nuclear(update_Lr);
norm1_term = para.sigma*norm(update_M,1);
trace_term = para.gamma*trace(update_Lr'*update_G*update_Lr);
Fnorm_term = para.beta*(norm(update_G,'fro')^2);
obj_value = nuclear_term + norm1_term + trace_term + Fnorm_term;

Lr_error = norm((update_Lr-Lr_0),'fro')/norm(Lr_0,'fro');
G_error = norm((update_G-Lap),'fro')/norm(Lap,'fro');

%calculate whether it gives a good estimation there is a edge between two nodes or not.
Num_nodes = size(Lap,1);
Num_edge_updateG = abs(update_G)>0;
Num_edge_Lap = abs(Lap)>0;
Num_diff_edge = nnz(Num_edge_Lap - Num_edge_updateG)-Num_nodes;
Diff_edge_error = Num_diff_edge/((Num_nodes*Num_nodes-Num_nodes)/2);

% print out a list
fprintf('%g   |%.4f   |%.4f   |%.4f   |%.4f   |%.4f   |%.4f   |%.4f   |%.4f          |%g   | \n',...
    iteration,nuclear_term,norm1_term,trace_term,Fnorm_term,obj_value,Lr_error,G_error,Diff_edge_error,rank);

end