% =========================================================================
% Plotting controller outputs
% =========================================================================
clc;
clear all;
close all
addpath(genpath('C:\Users\maertins\Documents\DeSiO\Git\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));
model = fsi_readmodel;
values = load([model.strSimName '_servo.dres']);

time       = values(:,1); 
pitch1     = values(:,2);
pitch1     = pitch1*(180/pi);
pitch2     = values(:,3);
pitch3     = values(:,4);
mgenerator = values(:,5);
omega_hub  = values(:,6); %[rad/s]
omega_hub  = omega_hub*(60/(2*pi));

M_hub      = values(:,7);
thrust     = values(:,8);
v_wind     = values(:,9);


figure 
plot(time,omega_hub)
title('Rotation speed, v_w = ')
xlabel('time [s]')
ylabel('\omega_{hub} [RPM]')

figure 
plot(time,pitch1)
title('Pitch angle')
xlabel('time [s]')
ylabel('pitch [°]')

figure 
plot(time,mgenerator)
title('Generator torque')
xlabel('time [s]')
ylabel('generator torque [Nm]')

