% plot the aerodynamic forces and moments over time

clc; 
clear all; 
close all;

str_software_path = 'C:\Arbeit_lokal\DeSiO';

strpath = fileread([str_software_path '\DeSiO_path.yaml']);
strpath = textscan(strpath,'%s');
addpath(genpath(strpath{1}{8}));

currDir = cd;
strfilename = 'wt';

nx = [1;0;0];
ny = [0;1;0];
nz = [0;0;1];

% ADAPT FOR EACH SIMULATION!
% geometric parameters
theta = 6.0*pi/180;                         % tilt angle
n_axis_r = [cos(theta); 0; -sin(theta)];    % rotation axis
x0       = [-11.96539128; 0; 149.99307271]; % reference point

% calculate and plot lift coeficient
model = uvlm_readmodel();
t   = fun_load_file([model.strSimName '_uvlm_t.dres']);
t = t(:,1);
dp  = fun_load_file([model.strSimName '_uvlm_dps_cp.dres']);
qs  = fun_load_file([model.strSimName '_uvlm_qs_nodal.dres']);
[F] = uvlm_f_surface(model,t,qs,dp);
[Mi] = uvlm_m_surface_i(model,t,qs,dp,x0);
[M] = uvlm_m_surface(model,t,qs,dp,x0);

for i = 1:size(F,1)
   fx(i) = nx'*F(i,1:3)';
   fy(i) = ny'*F(i,1:3)';
   fz(i) = nz'*F(i,1:3)';
end

for i = 1:size(M,1)
   mx(i) = nx'*M(i,1:3)';
   my(i) = ny'*M(i,1:3)';
   mz(i) = nz'*M(i,1:3)';
end

M_n = M*n_axis_r;
F_n = F*n_axis_r;

figure(); hold on; grid on;
xlabel('step'); ylabel('F_x');
plot(t,F_n,'-k','linewidth',2); 
plot(t,fx,'--r','linewidth',2); 
plot(t,fy,'--g','linewidth',2); 
plot(t,fz,'--b','linewidth',2); 
legend('F_{ax}','F_x','F_y','F_z','box','off');
saveas(gcf, 'AeroForces.png');

figure(); hold on; grid on;
xlabel('step'); ylabel('M_x');
plot(t,M_n,'-k','linewidth',2); 
plot(t,mx,'--r','linewidth',2);
plot(t,my,'--g','linewidth',2);
plot(t,mz,'--b','linewidth',2);
legend('M_{ax}','M_x','M_y','M_z','box','off');
saveas(gcf, 'AeroMoments.png');


% function uvlm_plot_CL_time()
% =================================================================================================================
    clc; 
    clear all; 
    close all;

    str_software_path = 'C:\Users\hente\Desktop';

    strpath = fileread([str_software_path '\DeSiO_path.yaml']);
    strpath = textscan(strpath,'%s');
    addpath(genpath(strpath{1}{8}));
    
    currDir = cd;
    strfilename = 'wt';
    
    nx = [1;0;0];
    ny = [0;1;0];
    nz = [0;0;1];
    
    theta = 0.0;
    n_axis_r = [cos(theta);0;-sin(theta)];
    x0       = [0;0;0];

    % calculate and plot lift coeficient
    model = uvlm_readmodel();
    t   = fun_load_file([model.strSimName '_uvlm_t.dres']); t = t(:,1);
%     if not(exist('M.mat'))
        dp  = fun_load_file([model.strSimName '_uvlm_dps_cp.dres']);
        qs  = fun_load_file([model.strSimName '_uvlm_qs_nodal.dres']);
        [F] = uvlm_f_surface(model,t,qs,dp);
        [Mi] = uvlm_m_surface_i(model,t,qs,dp,x0);
        [M] = uvlm_m_surface(model,t,qs,dp,x0);
%         save('M');
%         save('F');
%     else
%         F = load('F.mat'); F = F.F;
%         M = load('F.mat'); M = M.M;
%     end
% for i = 1:size(Mi,2)
%     m = Mi(i).m;
%     figure(); hold on; grid on;
%     title(['surf' num2str(i)]);
%     xlabel('step'); ylabel('M_x');
%     plot(t,m(:,1),'-r','linewidth',2);
%     plot(t,m(:,2),'-g','linewidth',2);
%     plot(t,m(:,3),'-b','linewidth',2);
%     legend('M_x','M_y','M_z','box','off');    
% end

    for i = 1:size(F,1)
       fx(i) = nx'*F(i,1:3)';
       fy(i) = ny'*F(i,1:3)';
       fz(i) = nz'*F(i,1:3)';
    end
    
    for i = 1:size(M,1)
       mx(i) = nx'*M(i,1:3)';
       my(i) = ny'*M(i,1:3)';
       mz(i) = nz'*M(i,1:3)';
    end
    
    M_n = M*n_axis_r;
    F_n = F*n_axis_r;
    
    figure(); hold on; grid on;
    xlabel('step'); ylabel('F_x');
    plot(t,F_n,'-k','linewidth',2); 
    plot(t,fx,'-r','linewidth',2); 
    plot(t,fy,'-g','linewidth',2); 
    plot(t,fz,'-b','linewidth',2); 
    legend('F_[ax}','F_x','F_y','F_z','box','off');
    
    figure(); hold on; grid on;
    xlabel('step'); ylabel('M_x');
    plot(t,M_n,'-k','linewidth',2); 
    plot(t,mx,'-r','linewidth',2);
    plot(t,my,'-g','linewidth',2);
    plot(t,mz,'-b','linewidth',2);
    legend('M_{ax}','M_x','M_y','M_z','box','off');
return
    dlmwrite(['f_test.txt'],[fx',fy',fz'],'delimiter','\t','precision',17);
    dlmwrite(['m_test.txt'],[mx',my',mz'],'delimiter','\t','precision',17);
% =================================================================================================================return

