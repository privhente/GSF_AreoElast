function compare_eigenforms_mac()
clc
clear all
close all

addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));

strfilename = 'test';

scale  = 5.0;
inz_ef = [1:3];

% Read simulation 1
curDir = cd;
model1 = structure_readmodel;
[ef1,ev1,q0] = get_eigenforms(model1,1.0);
plot_eigenforms(ev1(inz_ef),ef1(inz_ef,:),q0,model1);

% Read simulation 2
curDir = cd;
model2 = structure_readmodel;
[ef2,ev2] = get_eigenforms(model2,1.0);

% Modal Assurance Criterion
mac=MAC(ef1,ef2);

% Plot MAC-Data
f1 = figure(); spy3D(mac,1,1); title('MAC');
colormap('jet');
colorbar;
axis([0 size(mac,1)+1 0 size(mac,2)+1 0 1]);
% savefig(f1,[strfilename '_eigenform_MAC.fig']);
% print([strfilename '_eigenform_MAC'],'-dpng','-r0');
return

function plot_eigenforms(ev,ef,q0,model)
    str_exp = [model.strSimName '_model;'];
    eval(str_exp);
    for j = 1:length(ev)
        figure(); grid on; hold on; axis equal; view(45,20);
        title(['Eigenform:' num2str(j) ' eigenvalue: ' num2str(ev(j))]);
        xlabel('x'); ylabel('y'); zlabel('z');
        n0 = 0;
        for i = 1:size(model.beams,2)
            for k = 1:size(model.beams(i).connectivity,1)
                coord1  = ef(j,[3*(model.beams(i).connectivity(k,1)+n0-1)+1:3*(model.beams(i).connectivity(k,1)+n0-1)+3]);
                coord2  = ef(j,[3*(model.beams(i).connectivity(k,2)+n0-1)+1:3*(model.beams(i).connectivity(k,2)+n0-1)+3]);
                coord10 = q0([3*(model.beams(i).connectivity(k,1)+n0-1)+1:3*(model.beams(i).connectivity(k,1)+n0-1)+3]);
                coord20 = q0([3*(model.beams(i).connectivity(k,2)+n0-1)+1:3*(model.beams(i).connectivity(k,2)+n0-1)+3]);
                plot3([coord1(1) coord2(1)],[coord1(2) coord2(2)],[coord1(3) coord2(3)],'-b','linewidth',2);
                plot3([coord10(1) coord20(1)],[coord10(2) coord20(2)],[coord10(3) coord20(3)],'-k','linewidth',1);
            end
            n0 = n0 + model.beams(i).nnodes;
        end
    end
    
function [ef,ev,q0] = get_eigenforms(model,scale)
    step  = fun_load_file([model.strSimName '_steps.dres']);
    inzm  = find(step(:,2)==3);
    t     = fun_load_file([model.strSimName '_t.dres']);
    q     = fun_load_file([model.strSimName '_q.dres']); 
    q0    = q(1,:); dq = q(inzm,:);
    
    inz   = [(1:12:size(dq,2))',(2:12:size(dq,2))',(3:12:size(dq,2))'];
    inz   = reshape(inz', [], 1);

    ev = t(inzm,1);
    ef = [];
    q0 = q0(inz);
    for j = 1:length(ev)
        ef(j,:) = q0 + dq(j,inz)*scale;
    end
    
   return
   
function mac=MAC(phi1,phi2)
    % phi: matrix of the identified mode shapes
    % mac: MAC matrix
    for I=1:size(phi1,1)
        for J=1:size(phi2,1)
            mac(I,J)=Mac(phi1(I,:),phi2(J,:));
        end
    end
return

function mAc=Mac(Phi1,Phi2)
    % This function calculates mac between phi1 and phi2
    mAc= (abs(dot(Phi1,Phi2)))^2/((dot(Phi1,Phi1))*(dot(Phi2,Phi2)));
return

