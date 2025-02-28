function uvlm_evaluate_f_m()
clc;
clear all;
close all;

currDir = cd;
addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));

strname = '5mw';
theta = 0.10471;
omega = 5*2*pi/60;
X0 = [-11.9654;0;5.6071];
n_rotor = [cos(theta);0;-sin(theta)];

flag_plot_each = 0;

markersize = 8;
linewidth  = 2;
fontsize   = 14;
fontname   = 'times';
fontweight = 'bold';

pitch_nrel  = 3.86;
torque_nrel = 0.077;
power_nrel  = 0.037;
thrust_nrel = 0.205;

% evaluation resutls
file    = dir(currDir); 
strpath = {}; arrfolder = {}; m = [];
for i = 1:size(file,1)
    if file(i).isdir == 1 && length(file(i).name)>2
        strpath = file(i).name;
        i_1 = findstr(strpath,'_p');
        i_2 = findstr(strpath,'_om');
        if not(isempty(i_1))
            var1 =  str2num(strpath(i_1+2:i_2-1));
            m(end+1,1) = [var1];
            arrfolder{end+1} = strpath;
        end        
    end
end

[val, inzm] = sort(m(:,1)); 
m = m(inzm,:);
arrfolder = arrfolder(inzm);

% inzfind = find(m>=14 & m<=15);
% [inzfind,i1] = find(m==sort([10,12,14,14.46,15]));
% m = m(inzfind);
% arrfolder = arrfolder(inzfind);

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
        plot(res.M*n_rotor/1e6,'-','linewidth',1.0,'color',rgbcolor(i,:)); 

        figure(fig_f);
        xlabel('steps'); ylabel('F_{thrust} in MNm');
        plot(res.F*n_rotor/1e6,'-','linewidth',1.0,'color',rgbcolor(i,:)); 

        str_leg{end+1} = ['pitch ' num2str(m(i))];
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

% interpolate to get thrust, torque and power corresponding to pitch_nrel:
arr_torque = arr_m(:,1:3)*n_rotor/1e6;
arr_power  = arr_m(:,1:3)*n_rotor/1e6*omega;
arr_thrust = arr_f(:,1:3)*n_rotor/1e6;

torque_pitch_nrel = interp1(m,arr_torque,pitch_nrel,'spline')
power_pitch_nrel  = interp1(m,arr_power,pitch_nrel,'spline')
thrust_pitch_nrel = interp1(m,arr_thrust,pitch_nrel,'spline')

[m_val, torque_val] = fun_interpol_function(m,arr_torque,1000,2); [pitch_torque] = fun_find_xofy(m_val,torque_val,torque_nrel,2);
[m_val, power_val]  = fun_interpol_function(m,arr_power,1000,2); [pitch_power] = fun_find_xofy(m_val,power_val,power_nrel,2);
[m_val, thrust_val] = fun_interpol_function(m,arr_thrust,1000,2); [pitch_thrust] = fun_find_xofy(m_val,thrust_val,thrust_nrel,1);

fig = figure(); hold on; grid on;
xlabel('pitch in °')
yyaxis left; ylabel('M_{rotor} in MNm and P_{rotor} in MW');
h(1) = plot(m,arr_m(:,1:3)*n_rotor/1e6,'--b','marker','x','markersize',markersize,'linewidth',linewidth);
h(2) = plot(m,arr_m(:,1:3)*n_rotor/1e6*omega,'-b','marker','o','markersize',markersize,'linewidth',linewidth);
yyaxis right; ylabel('F_{thrust} in MN');
h(3) = plot(m,arr_f(:,1:3)*n_rotor/1e6,'--r','marker','.','markersize',markersize,'linewidth',linewidth);
yyaxis left;
str_leg = {'M_{rotor}','P_{rotor}','F_{thrust}'};

% additional points
for i = 1:length(pitch_torque); plot(pitch_torque(i),torque_nrel,'marker','x','color','k','markersize',markersize); end
for i = 1:length(pitch_power); plot(pitch_power(i),power_nrel,'marker','x','color','k','markersize',markersize); end
for i = 1:length(pitch_thrust); plot(pitch_thrust(i),thrust_nrel,'marker','x','color','k','markersize',markersize); end

legend(h(1:3),str_leg,'location','best','orientation','horizontal','fontsize',fontsize,'box','off');
set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
% fig.Position = [100 100 1200 400]; pbaspect auto;
savefig(fig,[strname '_f_m' '.fig']); print(['all_f_m'],'-dpng', '-r500');    

delta_pitch_torque = pitch_nrel - max(pitch_torque)
delta_pitch_power  = pitch_nrel - max(pitch_power)
delta_pitch_thrust = pitch_nrel - max(pitch_thrust)
return

function [x_vals] = fun_find_xofy(x_val, y_val, val,p)
    [val, inz_y] = sort(abs(y_val - val));
    x_vals = x_val(inz_y(1:p));
return

function [interpolatedX, interpolatedY] = fun_interpol_function(arr_valx,arr_valy,n_inter,p)
    coeffs = polyfit(arr_valx,arr_valy,p);
    interpolatedX = linspace(min(arr_valx), max(arr_valx), n_inter); 
    interpolatedY = polyval(coeffs, interpolatedX);
    return
