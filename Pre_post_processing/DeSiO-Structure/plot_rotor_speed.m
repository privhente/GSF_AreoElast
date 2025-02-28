% script for calculating rotor speed from director-based formulation
clc;
clear; 
close all;
addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));

model = fsi_readmodel;
q     = load([model.strSimName '_q.dres']);
s     = load([model.strSimName '_v.dres']);
time  = load([model.strSimName '_t.dres']); 
t     = time(:,1);

nd    = -[-7.03233176e-01;7.03233176e-01;1.04528463e-01]; % axis in local director cos of node
node  = 2;

[res] = fun_extract_rotor_speed(q,s,t,node,nd); 
omega = res.val;

fig = figure(); hold on; grid on;
plot(t,omega,'-','lineWidth',2,'color','r');
ylabel('\omega [rad/s]'); xlabel('time [s]')
title(['rotor speed vs. time at node' num2str(node)]);
legend('\omega');
set(gca,'fontsize',12,'fontweight','bold');
return