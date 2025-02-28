% =========================================================================
% Plotting DeSiO Coordinate-Solution
% =========================================================================
clc;
clear all;
close all;

addpath('..\..\files_for_prepostprocessing\DeSiO-Structure\');

% Loading of DeSiO result files
q    = load(['solution_fsi_q.dres']);
time = load(['solution_fsi_t.dres']);

% Function to extract dof-solution for node from solution.m file
node  = [39];
for i = 1:length(node)
    [u] = get_DeSiO_dof_solu(node(i),2,q);
    figure(); hold on; grid on;
    plot(time(:,1),u(:,1),'-r');
    plot(time(:,1),u(:,2),'-b');
    plot(time(:,1),u(:,3),'-k');
    ylabel('displacement m'); xlabel('time s')
    title(['displacement vs. time of node' num2str(node(i))]);
    legend('ux','uy','uz');
end
set(gca,'fontsize',12,'fontweight','bold');
% print('displ_time','-dpng', '-r500');
return