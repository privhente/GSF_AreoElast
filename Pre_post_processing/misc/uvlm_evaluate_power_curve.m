function uvlm_evaluate_power_curve()
clc;
clear all;
% close all;

currDir = cd;

str_software_path = 'C:\Users\hente\Desktop';
strpath = fileread([str_software_path '\DeSiO_path.yaml']);
strpath = textscan(strpath,'%s');
addpath(genpath(strpath{1}{8}));
 
theta = 0.10471;
omega = 7.527*2*pi/60;
X0 = [-12.0307;0;5.6135];
n_rotor = [cos(theta);0;-sin(theta)];
strname = 'IEA_15MW'

flag_plot_each = 1;

markersize = 8;
linewidth  = 2;
fontsize   = 14;
fontname   = 'times';
fontweight = 'bold';

% evaluation resutls
file    = dir(currDir); 
strpath = {}; arrfolder = {}; m = [];
for i = 1:size(file,1)
    if file(i).isdir == 1 && length(file(i).name)>2
        strpath = file(i).name;
        i_1 = findstr(strpath,'_p'); i_1 = i_1(end);
        i_2 = findstr(strpath,'_om');
        i_3 = findstr(strpath,'_v');
        if not(isempty(i_1))
            var1 = str2num(strpath(i_1+2:i_2-1));
            var2 = str2num(strpath(i_2+3:i_3-1));
            var3 = str2num(strpath(i_3+2:end));
            m(end+1,1:3) = [var3,var1,var2];
            arrfolder{end+1} = strpath;
        end        
    end
end

[val, inzm] = sort(m(:,1)); 
m = m(inzm,:);
arrfolder = arrfolder(inzm);

for i = 1:length(m)
   folder_name{i,:} = arrfolder{i};
end

if flag_plot_each == 1
    fig_m = figure(); grid on; hold on;
    fig_f = figure(); grid on; hold on;
end

arr_f = [];
arr_m = [];
str_leg = {};

rgbcolor = RGB_Color;
for i = 1:length(m)
    copyfile([currDir '\uvlm_get_m_f_time.m'],[currDir '\' folder_name{i,:} '\DeSiO-Aero\']);
    cd([currDir '\' folder_name{i,:} '\DeSiO-Aero']);
    disp(['... reading results in...\' folder_name{i,:}]);
    if not(exist('res.mat'))
        res = uvlm_get_m_f_time(X0);
    else
        res = load('res.mat'); res = res.res;
    end
    arr_f(i,1:3) = [res.F(end-2,1:3)];
    arr_m(i,1:3) = [res.M(end-2,1:3)];

    if flag_plot_each == 1
        figure(fig_m);
        xlabel('steps'); ylabel('M_{rotor} in MNm');
        plot(res.M*n_rotor,'-','linewidth',1.0,'color',rgbcolor(i,:)); 

        figure(fig_f);
        xlabel('steps'); ylabel('F_{thrust} in MNm');
        plot(res.F*n_rotor,'-','linewidth',1.0,'color',rgbcolor(i,:)); 

        str_leg{end+1} = ['vel ' num2str(m(i,1))];
    end
    
end
cd(currDir)

if flag_plot_each == 1
    figure(fig_f);
    legend(str_leg,'location','bestoutside','orientation','vertical','NumColumns',2,'fontsize',fontsize,'box','off');
    set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
    fig_f.Position = [100 100 1200 400]; pbaspect auto;
    % savefig(fig_f,['all_f' '.fig']); print(['all_f'],'-dpng', '-r500');    

    figure(fig_m);
    legend(str_leg,'location','bestoutside','orientation','vertical','NumColumns',2,'fontsize',fontsize,'box','off');
    set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
    fig_m.Position = [100 100 1200 400]; pbaspect auto;
    % savefig(fig_m,['all_f' '.fig']); print(['all_f'],'-dpng', '-r500');    
end
clc;

% read power curve from file
fid = fopen('PC_IEA_15MW.opt');
while ~feof(fid)
    tline = fgetl(fid);
    cell_line = textscan(tline,'%10.3f %s');
    num_wind_speed = cell_line{1}(1);
    for i = 1:num_wind_speed
        tline = fgetl(fid);
        arr_val(i,:) = str2num(tline);
    end
end
fclose(fid);

arr_omega  = m(:,3)*2*pi/60;
arr_torque = arr_m(:,1:3)*n_rotor;
arr_power  = arr_torque.*arr_omega;
arr_thrust = arr_f(:,1:3)*n_rotor;

fig_p = figure(); hold on; grid on;
xlabel('velocity in m/s')
ylabel('P_{rotor} in MW');
plot(m(:,1),arr_power/1e6,'--b','marker','x','markersize',markersize,'linewidth',linewidth);
plot(arr_val(:,1),arr_val(:,6),'--k','marker','x','markersize',markersize,'linewidth',linewidth);
legend('DeSiO','NREL');
set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname,'XLim',[m(1,1),m(end,1)+1]);%,'YTick',[0:1:14]);
fig_p.Position = [100 100 1200 400]; pbaspect auto;
savefig(fig_p,[strname '_power' '.fig']); print([strname '_power'],'-dpng', '-r500');    

fig_f = figure(); hold on; grid on;
xlabel('velocity in m/s')
ylabel('F_{thrust} in MN');
plot(m(:,1),arr_thrust/1e6,'--b','marker','x','markersize',markersize,'linewidth',linewidth);
plot(arr_val(:,1),arr_val(:,5),'--k','marker','x','markersize',markersize,'linewidth',linewidth);
legend('DeSiO','NREL');
set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname,'XLim',[m(1,1),m(end,1)+1]);%,'YTick',[0:0.2:2]);
fig_f.Position = [100 100 1200 400]; pbaspect auto;
savefig(fig_f,[strname '_thrust' '.fig']); print([strname '_thrust'],'-dpng', '-r500');    
return