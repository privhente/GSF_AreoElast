% =========================================================================
% Plotting DeSiO DOF-Solution
% =========================================================================
clc;
clear all;
close all

addpath('G:\10_DeSiO\DeSiO-Structure\Pre_Postprocessing_files');
addpath('G:\10_DeSiO\DeSiO-Aero\pre_post_processing\versions\v03_082021');

% Loading of DeSiO result files
q    = load(['solution_fsi_q.dres']);
s    = load(['solution_fsi_v.dres']);
time = load(['solution_fsi_t.dres']);
step = load(['solution_fsi_steps.dres']);

% Function to extract dof-solution for node from solution.m file
node  = [51];
for i = 1:length(node)
    [u] = get_DeSiO_dof_solu(node(i),1,q);
    [v] = get_DeSiO_dof_solu(node(i),1,s);
    
    inz = [12*(node(i)-1)+1:12*(node(i)-1)+12];
    
    i10 = [0,0,0];
    i20 = [0,1,0];
    i30 = [0,0,1];
    
    d10 = q(1,inz(4:6));
    d20 = q(1,inz(7:9));
    d30 = q(1,inz(10:12));
    
    w1j = s(1,inz(4:6));
    w2j = s(1,inz(7:9));
    w3j = s(1,inz(10:12));
    
    alpha = 0.5*(cross(d10,d10-i10) + cross(d20,d20-i20) + cross(d30,d30-i30));
    alpha_dot(1,1:3) = 0.5*(cross(d10,w1j) + cross(d20,w2j) + cross(d30,w3j));
    for j = 1:size(time,1)-1
        d1j = q(j+1,inz(4:6));
        d2j = q(j+1,inz(7:9));
        d3j = q(j+1,inz(10:12));
        
        deltad1 = d1j-d10;
        deltad2 = d2j-d20; 
        deltad3 = d3j-d30; 
        
        w1j = s(j+1,inz(4:6));
        w2j = s(j+1,inz(7:9)); 
        w3j = s(j+1,inz(10:12));
        
        dalpha = 0.5*(cross(d1j,deltad1) + cross(d2j,deltad2) + cross(d3j,deltad3));
        alpha(j+1,1:3)     = alpha(j,1:3) + dalpha(1:3);
        alpha_dot(j+1,1:3) = 0.5*(cross(d1j,w1j) + cross(d2j,w2j) + cross(d3j,w3j));
        
        d10 = d1j; d20 = d2j; d30 = d3j;
    end
    
    % plot displacement and angle
    figure(1);
    subplot(2,2,1); hold on; grid on;
    tit = title(['$$\textrm{u vs. time }$$'])
    set(tit,'interpreter','latex','fontWeight','bold','fontsize',12);
    xlabel('time s'); ylabel('displacement ft');
    plot(time(:,1),u(:,1),'-r');
    plot(time(:,1),u(:,2),'-g');
    plot(time(:,1),u(:,3),'-b');
    legend('u_1','u_2','u_3');
    
    subplot(2,2,3); hold on; grid on;
    tit = title(['$$\phi_1 \textrm{ vs. time }$$'])
    set(tit,'interpreter','latex','fontWeight','bold','fontsize',12);
    xlabel('time s'); ylabel('rotation rad');
    plot(time(:,1),alpha(:,1),'-r');
    plot(time(:,1),alpha(:,2),'-g');
    plot(time(:,1),alpha(:,3),'-b');
    legend('\phi_1','\phi_2','\phi_3');
    
    subplot(2,2,2); hold on; grid on;
    tit = title(['$$u_3 \textrm{ vs. } v_3$$'])
    set(tit,'interpreter','latex','fontWeight','bold','fontsize',12);
    xlabel('displacement ft'); ylabel('velocity ft/s');
    plot(u(:,1),v(:,1),'-b');
    
    subplot(2,2,4); hold on; grid on;
    tit = title(['$$\phi_1 \textrm{ vs. } \dot{\phi_1}$$'])
    set(tit,'interpreter','latex','fontWeight','bold','fontsize',12);
    xlabel('rotation rad'); ylabel('angular velocity rad/s');
    plot(alpha(:,3),alpha_dot(:,3),'-r');
    
    % plot fft
    ft1 = fft_x(2,u(:,1),time(:,1),1);
    ft2 = fft_x(2,alpha(:,3),time(:,1),2);
    axis([0 10 0 1]);
    legend('bending', 'torsion');
end
figure(1); savefig(1,'state_res.fig');
figure(2); savefig(2,'fft.fig');
return
