function checking_twist_chord()
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

% plot settings
markersize = 8;
linewidth  = 2;
fontsize   = 14;
fontname   = 'times';
fontweight = 'bold';

% filename for model data from yaml file for reference data
str_dir_yaml = 'DTU_10MW.yaml';

% filenames for models to compare
strpath0 = 'C:\Users\hente\Desktop\DeSiO-simulations\DTU_10MW\10mw_notwist_p0_om9.6_v12\DeSiO-Aero';
strpath1 = 'C:\Users\hente\Desktop\DeSiO-simulations\DTU_10MW\10mw_twist_p0_om9.6_v12\DeSiO-Aero';
strname  = 'cmp_twist';

% strpath0 = 'C:\Users\hente\Desktop\DeSiO-simulations\IEA_NREL_15MW\15mw_aero\newyaml_old_wioc_15mw_pc_tilt_notwist\DeSiO';
% strpath1 = 'C:\Users\hente\Desktop\DeSiO-simulations\IEA_NREL_15MW\15mw_aero\newyaml_old_wioc_15mw_pc_tilt_twist\DeSiO';
% strname  = 'newyaml_old_wioc_15mw_pc_tilt_twist';

% strpath0 = 'C:\Users\hente\Desktop\DeSiO-simulations\IEA_NREL_15MW\15mw_aero\oldyaml_old_wioc_15mw_pc_tilt_notwist\DeSiO';
% strpath1 = 'C:\Users\hente\Desktop\DeSiO-simulations\IEA_NREL_15MW\15mw_aero\oldyaml_old_wioc_15mw_pc_tilt_twist\DeSiO';
% strname  = 'oldyaml_old_wioc_15mw_pc_tilt_twist';

% strpath0 = 'C:\Users\hente\Desktop\DeSiO-simulations\IEA_NREL_15MW\15mw_aero\oldyaml_new_wioc_15mw_pc_tilt_notwist\DeSiO';
% strpath1 = 'C:\Users\hente\Desktop\DeSiO-simulations\IEA_NREL_15MW\15mw_aero\oldyaml_new_wioc_15mw_pc_tilt_twist\DeSiO';
% strname  = 'oldyaml_new_wioc_15mw_pc_tilt_twist';

% hub coordinates in global coordinate system
X0 = [-7.073;0;3.3688];

% blade length
L  = 86.4;

% get model from yaml file
model_ref    = YAML.read(str_dir_yaml);
pitch_angle  = model_ref.environment.pitch_angle;
twist_grid   = model_ref.components.blade.DeSiO.uvlm.twist.grid;
twist_values = model_ref.components.blade.DeSiO.uvlm.twist.values;
chord_grid   = model_ref.components.blade.DeSiO.uvlm.chord.grid;
chord_values = model_ref.components.blade.DeSiO.uvlm.chord.values;
phi_tilt     = model_ref.components.nacelle.DeSiO.drivetrain.uptilt*180/pi;
nblades      = model_ref.assembly.number_of_blades;

% read model data to compare
cd(strpath0);
model0 = uvlm_readmodel();

cd(strpath1);
model = uvlm_readmodel();
cd(currDir);

% get model data for both models to compare
res0 = get_modeldata2compare(model0,X0);
res  = get_modeldata2compare(model,X0);

% compare data of both model
n_rotor  = [cos(phi_tilt*pi/180);0;-sin(phi_tilt*pi/180)];

aa = 0; hp = []; hch = [];
strlegend = {'reference from yaml'};
f1 = figure(); hold on; grid on; title(['twist in span in rad']); xlabel('grid'); ylabel('twist in rad');
htref = plot(twist_grid,twist_values,'-k','linewidth',linewidth);

fch = figure(); hold on; grid on; title(['chord in m' ]); xlabel('grid'); ylabel('chord in m');
hcref = plot(chord_grid,chord_values,'-k','linewidth',linewidth);

rgbcolor = RGB_Color;
a_rotor = 0;
for i = 1:size(res0,2)
    chord0 = res0(i).chord;
    chord  = res(i).chord;
    s0     = res0(i).s;
    s      = res(i).s;
    
    grid0  = s0./L;
    gridm  = s./L;

    % compare chord length in span direction
    ii = round(sin(pi*(i+1)/2),3);
    if round(sin(pi*(i+1)/2),3) == 0
        aa = aa + 1;
        a_rotor = (aa-1)*2*pi/nblades*180/pi;
        strlegend{end+1} = ['blade ' num2str(aa)];
        color_i = rgbcolor(i,:);
    end
    figure(fch); hch(end+1) = plot(gridm,chord,'--','color',color_i,'linewidth',linewidth); 
%     hd = plot(gridm,chord,'xr','markersize',markersize);
    
    % calculate twist angle between model0 and model in blade coordinate
    % system
    twist = [];
    for j = 1:size(res0(i).unit_a,2)
        unit_a0 = res0(i).unit_a(:,j);
        unit_a  = res(i).unit_a(:,j);
        
        e10 = unit_a0; e10 = e10/norm(e10);
        e20 = cross(e10,n_rotor); e20 = e20/norm(e20);
        e30 = cross(e10,e20);
        
        unit_a_e = [unit_a'*e10;unit_a'*e20;unit_a'*e30];
        twist(j) = -sign(unit_a_e(3))*acos(unit_a'*e10);
    end
    figure(f1); hp(end+1) = plot(grid0,twist,'--','color',color_i,'linewidth',linewidth);
end
cd(currDir);

figure(f1);
legend([htref, hp(2:2:end)],strlegend,'location','best','orientation','vertical','fontsize',fontsize,'box','off');
set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
savefig(f1,[strname '_comparison_twist' '.fig']); print([strname '_comparison_twist'],'-dpng', '-r500');    

figure(fch); legend(hch(2:2:end),strlegend);
legend([hcref, hch(2:2:end)],strlegend,'location','best','orientation','vertical','fontsize',fontsize,'box','off');
set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
savefig(fch,[strname '_comparison_chord' '.fig']); print([strname '_comparison_chord'],'-dpng', '-r500');    
end

function res = get_modeldata2compare(model,X0)
    i_2   = [0;1;0];
    for i = 1:size(model.surfaces,2)
        aa = 0;
        nn  = model.surfaces(i).nn;
        ne  = model.surfaces(i).ne;
        nnx = model.surfaces(i).nnx;
        nny = model.surfaces(i).nny;
        coords = model.surfaces(i).coord;
        coords(:,3) = coords(:,3);
        if round(sin(pi*(i+1)/2),3) == 0
            flag1 = 2; flag2 = 1; 
        else
            flag1 = 1; flag2 = +1; 
        end
        for j = 1:nny
            aa  = aa + 1;
            inz = nnx*(j-1)+1:nnx*(j-1)+nnx;
            i1  = inz(1);
            i2  = inz(1) + fix((inz(end)-inz(1))/flag1+flag2) -1;

            a1 = coords(i1,:)' - X0;
            a2 = coords(i2,:)' - X0;

            vecref = sum(coords(inz,:))'./length(inz) - X0;
            veca   = a2-a1;
            unit_a = veca/norm(a1-a2);
            chord  = norm(a1-a2);
            twist  = sign(unit_a(1))*acos(unit_a'*i_2);
            
            res(i).vecref(1:3,aa) = vecref;
            res(i).veca(1:3,aa)   = veca;
            res(i).unit_a(1:3,aa) = unit_a;
            res(i).chord(aa)      = chord;
            res(i).twist(aa) = twist;
        end
        if round(sin(pi*(i+1)/2),3) == 0
            lref0 = norm(res(i).vecref(:,1));
        end        
        s = norm(res(i).vecref(:,1));
        for k = 1:size(res(i).vecref,2)-1
            vecrefi       = res(i).vecref(:,k);
            vecrefip1     = res(i).vecref(:,k+1);
            ds            = norm(vecrefip1-vecrefi);
            s(k+1)        = s(k) + ds;
        end
        s = s - lref0;
        res(i).s = s;
    end
end