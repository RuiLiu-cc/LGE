function plot_Lr_Graph_eachStep(update_Lr,update_G)
% this function is to plot the matrix of Lr and graph matrix after each
% step finish to help better understand what is going on
%==========================================================================
%
% Author: Liu Rui, SUTD, 27 Feb 2018
%
%==========================================================================

figure(1);
subplot(2,3,5);imagesc(update_Lr);title('Step1: Updated low-rank matrix(Lr)');
subplot(2,3,6);imagesc(update_G);title('Step2: Updated Graph');

end