function [e_rs, e_is, es, es_d0, Xs] = evaluate_for_flutter_eigenvalue_analysis()
% =========================================================================
% =========================================================================
% clc;
clear all;
% close all;

addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));
model      = fsi_readmodel;
model_uvlm = uvlm_readmodel;
deltat     = model_uvlm.simulationsettings(2);

% loading q
% q  = load([model.strSimName '_q.dres']);
% q0 = q(1,:);

% reading matrix in sparse format
csr_colums_s     = load('csr_columns_S.dres');
csr_values_s     = load('csr_values_S.dres');
csr_rowIndices_s = load('csr_rowIndices_S.dres');

csr_colums_kae     = load('csr_columns_Kae.dres');
csr_values_kae     = load('csr_values_Kae.dres');
csr_rowIndices_kae = load('csr_rowIndices_Kae.dres');

csr_colums_m     = load('csr_columns_M.dres');
csr_values_m     = load('csr_values_M.dres');
csr_rowIndices_m = load('csr_rowIndices_M.dres');

indicesq = load('indicesq.dres');
indicesv = load('indicesv.dres');
indicesg = load('indicesg.dres');

% number of coordinates
nq = length(indicesq);
nv = length(indicesv);
ng = length(indicesg);

% dense matrices
% system matrix only structural term considered
S   = csr2full(csr_values_s,csr_rowIndices_s,csr_colums_s,nq+nv+ng,nq+nv+ng); %
Kae = csr2full(csr_values_kae,csr_rowIndices_kae,csr_colums_kae,nq,nq+nv); %
M   = csr2full(csr_values_m,csr_rowIndices_m,csr_colums_m,nq,nq); %

Kqq  = S(indicesq,indicesq);
Kqv  = S(indicesq,indicesv);
H    = S(indicesg,indicesq);
Kaeq = Kae(indicesq,indicesq); 
Kaev = Kae(indicesq,indicesv);

M_g = (Kqv + 0.5*Kaev)*deltat;
err_m = norm(M_g-M);
disp(['error in mass: ' num2str(err_m)]);

% nullspace of constraint matrix
Nt = null(H);
N  = Nt';

% solving structural system without aerodynamic matrices
K_red_s = N*(2*Kqq+Kaeq)*Nt; 
M_red   = N*M_g*Nt;
[Vs,es] = eig(K_red_s,M_red);
es = diag(es);
[inz] = find(es>=0); es = es(inz);
[es,inzs] = sort(es);
es = sqrt(es);
Vs = Vs(:,inzs);

% solving generalized eigenvalue problem: u = a*exp(lam*t)
% format short: (K + lam*C + lam^2*M)
% solving for reduced system, using null-space projection
flag_aer_d = 1;
% K_red = N*(2*Kqq+Kaeq)*Nt; 
K_red = N*(2*Kqq)*Nt; 
M_red = N*M_g*Nt;
C_red = N*Kaev*Nt*flag_aer_d;

if flag_aer_d==0
    [X,e] = eig(K_red,M_red);
    e = diag(e);
    e_r        = real(e);
    e_i        = imag(e);
    [inz]      = find(e_r >= 0);
    e_r0       = e_r(inz);
    e_i0       = e_i(inz);
    [e_r0_s,inzs] = sort(e_r0,'ascend');
    e_i0_s        = e_i0(inzs);
    Xs         = X(:,inz(inzs));
    e_rs       = sqrt(e_r0_s);
    e_is       = e_i0_s;
else
    [X,e]      = polyeig(K_red,-C_red,M_red);
    e_r        = real(e);
    e_i        = imag(e);
    [inz]      = find(e_i >= 0);
    e_r0       = e_r(inz);
    e_i0       = e_i(inz);
    X0         = X(:,inz);
    [tem,inzs] = sort(abs(e_i0));
    e_rs       = e_r0(inzs);
%     e_is       = sqrt(e_i0(inzs));
    e_is       = (e_i0(inzs));
    Xs         = X0(:,inzs);
end

K_red_s = N*(2*Kqq)*Nt; 
M_red   = N*M_g*Nt;
[Vs,es] = eig(K_red_s,M_red);
es = diag(es);
[inz] = find(es>=0); es = es(inz);
[es,inzs] = sort(es);
es_d0 = sqrt(es);
Vs = Vs(:,inzs);

disp(['structural eigenvalues: rad/s:                            ' num2str(real(es(1:5))')]);
disp(['aeroelastic (structural) eigenvalues, w/o damping: rad/s: ' num2str(real(es_d0(1:5))')]);
disp(['str+ae: Im(eigenvalues) rad/s:                            ' num2str(e_is(1:5)')]);
disp(['str+ae: Re(eigenvalues) rad/s:                            ' num2str(e_rs(1:5)')]);
disp(['det(C_red):' num2str(det(C_red))]);
disp(['det(Kae_qq):' num2str(det(N*Kaeq*Nt))]);
disp([' ']);


return
figure(); hold on; grid on; bar(e_rs); title('real (damping) - part of eigenvalue');
fig = figure(); hold on; grid on;
title('imaginar (frequency) - part of eigenvalue');
plot(e_rs(1:5),e_is(1:5));

% write the eigenmodes to files:
nef    = 5;
T      = 1;
deltat = T/10;
asteps = 0;
arr_steps = []; arr_time = [];
for i = 1:nef
    dphi    = real(Xs(:,i)); dphi = dphi/norm(dphi);
    dq(i,:) = (N'*dphi)';
    omega   = real(e_is(i));
    j = 0; qs = [];
    for time = 0:deltat:4*T
        j             = j + 1;
        asteps        = asteps + 1;
        arr_qs(i).qs(j,:) = q0 + dq(i,:).*sin(pi/2*time/T);
        % set steps for visualizing
        % step, i_simuType, istep, iteration
        arr_steps(end+1,1:4) = [asteps, 0, j, j];
        arr_time(end+1,1:2)  = [asteps,time];
    end
end

currDir = cd;
cd('..\');
currDir_temp = cd;
[~,name,~]=fileparts(currDir_temp);

kinDir  = [currDir_temp '\DeSiO-Aero'];
cd(kinDir);

save('arr_steps','arr_steps');
save('arr_time','arr_time');

delete('*.dres');
delete('*.txt');

fid1 = fopen(['solution_qfem.dres'],'w');
for i = 1:nef
    qss = arr_qs(i).qs;
    for j = 1:size(qss,1)
        fprintf(fid1,'%10.5e\t',qss(j,:)); fprintf(fid1,'\n');
    end
end
fclose(fid1);

simu.structprevjobname = 'solution';
simu.simType = 'kinematic';
simu.currDir = kinDir;
simu.caseDir = kinDir;
fun_writeDeSiOFSIInput(simu,[],[],[]);

copyfile([currDir '\surfaceinput.txt'],[kinDir '\surfaceinput.txt']);
copyfile([currDir '\wakeinput.txt'],[kinDir '\wakeinput.txt']);
copyfile([currDir '\beaminput.txt'],[kinDir '\beaminput.txt']);
copyfile([currDir '\constraint12input.txt'],[kinDir '\constraint12input.txt']);
copyfile([currDir '\fsi_input.txt'],[kinDir '\fsi_input.txt']);
return