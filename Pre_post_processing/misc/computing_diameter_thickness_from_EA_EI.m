function computing_diameter_thickness_from_EA_EI()
% =================================================================================================================
% Function to transform rotor data of openDAST model to WindIO-yaml for modeling wind energy
% 
% Author: Christian Hente
% Date: 30.03.2023
% =================================================================================================================
% =================================================================================================================
clc; 
clear all; 
close all; 
format short;

currDir = cd;

str_software_path = 'C:\Users\hente\Desktop';
strpath = fileread([str_software_path '\DeSiO_path.yaml']);
strpath = textscan(strpath,'%s');
addpath(genpath(strpath{1}{8}));

syms Da t rho A I E rhoA EI

f1 = rho* pi/4*(Da^2-(Da-2*t)^2)  - rhoA;
f2 = E*pi/64*(Da^4-(Da-2*t)^4) - EI;

so = solve([f1;f2],[Da;t]);

da = simplify(so.Da)
t = simplify(so.t)
return
