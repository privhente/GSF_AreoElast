function write_kinematics
% function convert_kinematic
close all;
clear all;
clc;

str_software_path = 'C:\Users\hente\Desktop';

strpath = fileread([str_software_path '\DeSiO_path.yaml']);
strpath = textscan(strpath,'%s');
addpath(genpath(strpath{1}{8}));

% =================================================================================================================
currDir = cd;
strfilename = 'solution'

% reading initial data from files:
deltat   = 0.025;
omega    = 7.25*2*pi/60;
nsteps   = 400;

theta = 0.0;
n_axis_r = [cos(theta);0;-sin(theta)];
X0       = [0;0;0];

% reading data from surfaceinput files
nn_g = 0;
fid = fopen('surfaceinput.txt','r');
    for i = 1:3
        tline = fgetl(fid);
    end
    numsurf = str2num(fgetl(fid));
    for k = 1:numsurf
        for i = 1:3
            tline = fgetl(fid);
        end
        tline = fgetl(fid);
        a = sscanf(tline,'%i %i %i %i',4);
        nn = a(1); ne = a(2); nnx = a(3); nny = a(4);
        for i = 1:3
            tline = fgetl(fid);
        end
        for i = 1:nn
            tline = fgetl(fid);
            tline = strrep(tline,'d','e');
            surface(k).coord(i,1:3) = sscanf(tline,'%f %f %f',3);
        end
        for i = 1:3
            tline = fgetl(fid);
        end
        for i = 1:ne
            tline = fgetl(fid);
            surface(k).inz(i,1:4) = sscanf(tline,'%i %i %i %i',4);
        end
    end
fclose(fid);

% calculate rotor rotation
tr_Step = 10;
for k = 1:tr_Step
    [HH] = Hermite (deltat*(k-1), 0, deltat*tr_Step);
    alpha_t(k) = HH.h3*(omega);
end

alpha_t(k+1:nsteps+1) = omega;
alpha_a(1:nsteps+1) = deltat*[0:nsteps].*alpha_t;

figure(); hold on; grid on;
yyaxis left; plot(alpha_t);
yyaxis right; plot(alpha_a);

% writing kinematics
vs = []; qs = []; alpha = 0;
for k = 1:nsteps+1
    alphat = alpha_t(k);
    alpha  = alpha_a(k);
    
%     n_axis_r = [1,0,0]';
    R_kin = cos(alpha)*eye(3) + sin(alpha)*skew(n_axis_r)+(1-cos(alpha))*n_axis_r*n_axis_r';
    R_kin_t = (-sin(alpha)*eye(3) + cos(alpha)*skew(n_axis_r)+(1+sin(alpha))*n_axis_r*n_axis_r');

    ie = 0; 
    for k1 = 1:size(surface,2)
        nn = size(surface(k1).coord,1);
        X_0 = [ones(nn,1)*X0(1),ones(nn,1)*X0(2),ones(nn,1)*X0(3)];
        X_C_n   = X_0 + (R_kin*(surface(k1).coord-X_0)')'; 
        X_C_n_t = alphat*(R_kin_t*(surface(k1).coord-X_0)')'; 
        ia = ie + 1;
        ie = ia + 3*nn - 1;
        qs(k,ia:ie) = reshape(X_C_n.',1,[]);
        vs(k,ia:ie) = reshape(X_C_n_t.',1,[]);
    end
end

% write kinematic position for nodes
fid1 = fopen([strfilename '_uvlm_qs.dres'],'w');
fid2 = fopen([strfilename '_uvlm_vs.dres'],'w');
    for i = 1:nsteps
        fprintf(fid1,'%10.5e\t',qs(i,:)); fprintf(fid1,'\n');
        fprintf(fid2,'%10.5e\t',vs(i,:)); fprintf(fid2,'\n');
    end
fclose(fid1);
fclose(fid2);

return

function [H] = Hermite (x, xi, xj)
    h = xj - xi;
    H.h1 = 1 - 3*(x - xi)^2 / (h^2) + 2*(x - xi)^3/(h^3);
    H.h2 = (x - xi) - 2/h*(x - xi)^2 + 1/(h^2)*(x - xi)^3;
    H.h3 = 3/h^2*(x - xi)^2 - 2/(h^3)*(x - xi)^3;
    H.h4 = -1/h*(x - xi)^2 + 1/(h^2)*(x - xi)^3;
return