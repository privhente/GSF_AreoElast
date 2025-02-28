% =========================================================================
% Plotting DeSiO Energy-Solution
% =========================================================================
clc;
clear all;
close all;

addpath('..\..\files_for_prepostprocessing\DeSiO-Structure\');

% Loading DeSiO result files
time       = load(['solution_fsi_t.dres']);
invariants = load(['solution_fsi_e.dres']);

figure(); hold on; grid on;
plot(time(:,1),invariants(:,1),'-b');
plot(time(:,1),invariants(:,2),'-r');
plot(time(:,1),invariants(:,3),'-g');
leg = legend('kinetic energy', 'potential energy', 'total energy');
ylabel('energy Nm'); xlabel('time s');
set(gca,'fontsize',12,'fontweight','bold');
title(['energy vs. time']);
print('energy_time','-dpng', '-r500');
return