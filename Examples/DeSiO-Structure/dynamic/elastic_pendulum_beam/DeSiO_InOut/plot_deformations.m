% =========================================================================
% Plotting DeSiO Coordinate-Solution
% =========================================================================
clc;
clear all;
% close all;

addpath('G:\10_DeSiO\DeSiO-Structure\Pre_Postprocessing_files');
addpath('G:\10_DeSiO\DeSiO-Aero\pre_post_processing\versions\v03_082021');

% Loading of DeSiO result files
model = strucure_readmodel;
q    = load([model.strSimName '_q.dres']);
time = load([model.strSimName '_t.dres']);

% Function to extract dof-solution for node from solution.m file
node  = [20];
for i = 1:length(node)
    [u] = get_DeSiO_dof_solu(node(i),1,q);
    figure(); hold on; grid on;
    plot(time(:,1),u(:,1),'-.r');
    plot(time(:,1),u(:,2),'-.b');
    plot(time(:,1),u(:,3),'-.k');
    ylabel('displacement m'); xlabel('time s')
    title(['displacement vs. time of node' num2str(node(i))]);
    legend('ux','uy','uz');
end
% print('displ_time','-dpng', '-r500');
return