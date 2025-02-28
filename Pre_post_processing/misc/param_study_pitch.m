% =================================================================================================================
% Function to transform WindIO to DeSiO-input for modeling wind energy
% converter.
% 
% Author: Christian Hente
% Date: 05.05.2022
% updated: 27.10.2022
% =================================================================================================================
function param_study_pitch()
% =================================================================================================================
clc; clear all; close all; format short e;
currDir = cd;

str_dir_yaml = 'IEA-15-240-RWT';
camber_line_Factor = 0.6;

fid = fopen(['PC_IEA_15MW.opt']);
% Wind speed [m/s], Pitch [deg], Rot. speed [rpm]
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

for i = 1:size(arr_val,1)
    
    vel   = arr_val(i,1);
    pitch = arr_val(i,2);
    omega = arr_val(i,3);
    
    %
    copyfile([str_dir_yaml '.yaml'],[str_dir_yaml '_temp.yaml']);
    fid1 = fopen([str_dir_yaml '_temp.yaml'],'r');
    fid2 = fopen('test.yaml','w');
        while ~feof(fid1)
            
            strline = fgetl(fid1);
            TF_pitch = contains(strline,'pitch_angle');
            TF_omega = contains(strline,'rotor_velocity:','IgnoreCase',true);
            TF_flag  = contains(strline,'flag_rotor_velocity:','IgnoreCase',true);
            TF_vel = contains(strline,'vinf:','IgnoreCase',true);
            if TF_pitch == 1
                newstrline = ['pitch_angle: ' num2str(pitch)];
                fprintf(fid2,'  %s\n',newstrline);
            elseif TF_omega == 1 && TF_flag ~= 1
                newstrline = ['rotor_velocity: ' num2str(omega)];
                fprintf(fid2,'  %s\n',newstrline);
            elseif TF_vel == 1           
                newstrline = ['vinf: ' num2str(vel)];
                fprintf(fid2,'  %s\n',newstrline);
            else
                fprintf(fid2,'%s\n',strline);
            end
            
            TF = contains(strline,'jobname:','IgnoreCase',true);
            if TF==1
                strjobname = textscan(strline,'%s %s');
                strjobname = strjobname{2};
                strjobname = strjobname{1};
            end
            TF = contains(strline,'uptilt:','IgnoreCase',true);
            if TF==1
                strline_n = textscan(strline,'%s %10.5f');
                theta = strline_n{2}*180/pi;
            end
            
        end
    fclose(fid1); fclose(fid2);
    delete([str_dir_yaml '.yaml']);
    copyfile('test.yaml',[str_dir_yaml '.yaml']);
    delete('test.yaml'); delete([str_dir_yaml '_temp.yaml']);
    WindIO_to_DeSiO_WT();
    
    % create kinematics
    wio_jobname = [strjobname '_clf' num2str(camber_line_Factor) '_p'  num2str(pitch) '_om' num2str(omega) '_v' num2str(vel)];
    copyfile('write_kinematics.m',[wio_jobname '\DeSiO-Aero\']);
    copyfile('uvlm_plot_m_f_time.m',[wio_jobname '\DeSiO-Aero\']);
    
    % get coordinates of hub
    cd([wio_jobname '\DeSiO-Structure\initial_data\']);
    fid1 = fopen('rigidbodyinput.txt','r');
    X0   = [0;0;0];
    while ~feof(fid1)
        strline  = fgetl(fid1);
        if contains(strline,'!! rigid body')
            if contains(strline,'hub phi')
                strline  = fgetl(fid1);
                strline = strrep(strline,'d','e');
                strline_r = str2num(strline);
                X0 = strline_r(1:3);
                break
            end
        end
    end
    fclose(fid1);

    cd(currDir);
    cd([wio_jobname '\DeSiO-Aero']);
    write_kinematics(omega,theta,X0);

    % change input for uvlm_plot_m_f_time.m file
    copyfile('uvlm_plot_m_f_time.m','uvlm_plot_m_f_time_temp.m');
    fid1 = fopen('uvlm_plot_m_f_time_temp.m','r');
    fid2 = fopen('test.m','w');
        while ~feof(fid1)
            strline  = fgetl(fid1);
            TF_om    = contains(strline,'omega   =');
            TF_theta = contains(strline,'theta   =','IgnoreCase',true);
            TF_X0    = contains(strline,'X0      =','IgnoreCase',true);
            if TF_om==1
                newstrline = ['omega = ' num2str(omega) '*2*pi/60;'];
                fprintf(fid2,'    %s\n',newstrline);
            elseif TF_theta==1
                newstrline = ['theta = ' num2str(theta*pi/180) ';'];
                fprintf(fid2,'    %s\n',newstrline);
            elseif TF_X0==1
                newstrline = ['X0 = [' num2str(X0(1)) ';' num2str(X0(2)) ';' num2str(X0(3)) '];'];
                fprintf(fid2,'    %s\n',newstrline);
            else
                fprintf(fid2,'%s\n',strline);
            end
        end
    fclose(fid1); fclose(fid2);
    delete('uvlm_plot_m_f_time.m');
    copyfile('test.m','uvlm_plot_m_f_time.m');
    delete('test.m'); delete('uvlm_plot_m_f_time_temp.m');
    cd(currDir);
end

return