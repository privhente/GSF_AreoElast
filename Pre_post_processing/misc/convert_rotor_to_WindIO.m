function model = convert_rotor_to_WindIO(info, model_inp)
% =================================================================================================================
% Function to transform rotor data of openDAST model to WindIO-yaml for modeling wind energy
% 
% Author: Christian Hente
% Date: 30.03.2023
% =================================================================================================================
% =================================================================================================================
currDir = cd;
% reading of old model
model = [];
if exist(info.strYAMLFileName,'file')
    disp(['reading information from existing yaml file: ' info.strYAMLFileName]);
    model = YAML.read(info.strYAMLFileName);
end

% merging two models
if exist('model_inp') && not(isempty(model_inp))
    disp(['merging information of existing model to new model']);
    inf_model_inp = whos('model_inp');
    [model] = recursive_asign_fields(model_inp,inf_model_inp.name,length(inf_model_inp.name),model);
end

% reading all files in director
wt_dir      = info.wt_dir;
airfoil_dir = info.airfoil_dir;

% openfast input filenames
str_wt_ElastoDynBlade = []; if isfield(info,'str_wt_ElastoDynBlade') str_wt_ElastoDynBlade = info.str_wt_ElastoDynBlade; end
str_wt_BeamDyn = [];        if isfield(info,'str_wt_BeamDyn')        str_wt_BeamDyn        = info.str_wt_BeamDyn; end
str_wt_AeroDynBlade = [];   if isfield(info,'str_wt_AeroDynBlade')   str_wt_AeroDynBlade   = info.str_wt_AeroDynBlade; end
str_wt_AeroDyn = [];        if isfield(info,'str_wt_AeroDyn')        str_wt_AeroDyn        = info.str_wt_AeroDyn; end
str_wt_ElastoDyn = [];      if isfield(info,'str_wt_ElastoDyn')      str_wt_ElastoDyn      = info.str_wt_ElastoDyn; end
str_wt_ElastoDyn_tow = [];  if isfield(info,'str_wt_ElastoDyn_tow')  str_wt_ElastoDyn_tow  = info.str_wt_ElastoDyn_tow; end
cd(wt_dir);

% checking which beam model shall be converted:
if ~isempty(str_wt_BeamDyn)
    str_wt_ElastoDyn = [];
end

% opening aerodyn file to obtain airfoil_data in model
if ~isempty(str_wt_AeroDyn)
    fid = fopen([str_wt_AeroDyn],'r');
    airfoil_data = []; 
    airfoil_ID = 0;
    if fid > 0
        disp(['reading AeroDyn information from ' str_wt_AeroDyn]);
        while ~feof(fid)
            tline = fgetl(fid);
            TF1 = contains(tline,'NumFoil','IgnoreCase',true);
            TF2 = contains(tline,'NumAFfiles','IgnoreCase',true);
            if (TF1 ~= 0 || TF2 ~= 0)
                str_cell_ne = textscan(tline,'%10.3f, %s');
                number_airfoils = str_cell_ne{1};
                for i = 1:number_airfoils
                    tline     = fgetl(fid);
                    inz_slash = 0;
                    if isempty(strfind(tline,'/'))
                        inz_slash = strfind(tline,'\');
                    else
                        inz_slash = strfind(tline,'/');
                    end
                    inz_ext = strfind(tline,'.dat');
                    if ~(isempty(inz_ext))
                        airfoil_ID = airfoil_ID + 1;
                        airfoil_data(airfoil_ID).name       = tline(inz_slash(end)+1:inz_ext-1);
                        airfoil_data(airfoil_ID).airfoil_ID = airfoil_ID;
                        airfoil_data(airfoil_ID).filename   = [tline(inz_slash(end)+1:inz_ext-1) '.dat'];
                        [airfoil_data] = fun_reading_airfoil_infos(info.airfoil_dir,airfoil_data,airfoil_ID);
                    end
                end
            end
        end
    fclose(fid);
    else
        disp('Warning: no AeroDyn.dat file found!');
    end
end

% opening aerodynBlade file to read:
% grid in span
% in-plane and out-of-plane distance from pitch axis to AC
% twist angle
% chord length
% airfoils
str_cell_ne = {}; ne = 0; grid = [];
if ~isempty(str_wt_AeroDynBlade)
    fid = fopen([str_wt_AeroDynBlade],'r');
    if fid > 0
        disp(['reading AeroDynBlade v15.00 information from ' str_wt_AeroDynBlade]);
        if ~isempty(str_wt_BeamDyn)
            disp(['  defining reference axis for aero model by local in-plane and out-of-plane offsets, .i.e. distances of AC to local (straight) blade axis  - see AeroDyn manual page 49']);
            disp(['  AC will be obtained from airfoil data in .coord file - see Aerodyn manual page 47, for AC']);
            disp(['  for using BeamDyn, pitch axis in WindIO_converter is equal to AC, i.e. the aifoil will be transfered from its AC to the aero reference line']);
        elseif ~isempty(str_wt_ElastoDynBlade)
            disp(['  defining reference axis for by local in-plane and out-of-plane offsets, .i.e. distances of AC to local (straight) blade axis  - see AeroDyn manual page 49']);
            disp(['  reference axis for beam model is the same as for aero model']);
            disp(['  pitch axis will be taken from ElastoDyn_Blade file']);
        end
        while ~feof(fid)
            tline = fgetl(fid);
            TF1 = contains(tline,'BlSpn','IgnoreCase',true);
            TF2 = contains(tline,'BlCrvAC','IgnoreCase',true);
            TF3 = contains(tline,'BlSwpAC','IgnoreCase',true);
            TF4 = contains(tline,'BlCrvAng','IgnoreCase',true);
            TF5 = contains(tline,'BlTwist','IgnoreCase',true);
            TF6 = contains(tline,'BlChord','IgnoreCase',true);
            TF7 = contains(tline,'BlAFID','IgnoreCase',true);
            if (TF1 ~= 0 && TF2 ~= 0 && TF3 ~= 0 && TF4 ~= 0 && TF5 ~= 0 && TF6 ~= 0 && TF7 ~= 0)
                tline = fgetl(fid);
                for j = 1:ne
                    tline = fgetl(fid);
                    arr_blade_values(j,:) = str2num(tline);
                end
                % define grid according to the span-discretization of openFAST-
                % input
                arr_span = arr_blade_values(:,1);
                grid = arr_span/max(arr_span);

                model.components.blade.DeSiO.uvlm.airfoil_position.grid = grid;
                for k = 1:ne
                    arr_airfoil_labels{k} = airfoil_data(arr_blade_values(k,7)).name;
                    arr_airfoil_AC(k)     = airfoil_data(arr_blade_values(k,7)).aerodynamic_center;
                end
                model.components.blade.DeSiO.uvlm.airfoil_position.labels = arr_airfoil_labels;

                model.components.blade.DeSiO.uvlm.chord.grid   = grid;
                model.components.blade.DeSiO.uvlm.chord.values = arr_blade_values(:,6)';

                model.components.blade.DeSiO.uvlm.twist.grid   = grid;
                model.components.blade.DeSiO.uvlm.twist.values = arr_blade_values(:,5)'.*pi/180;

                % here due to the converter, we put the pitch axis in the AC
                model.components.blade.DeSiO.uvlm.pitch_axis.grid = grid;
                model.components.blade.DeSiO.uvlm.pitch_axis.values = arr_airfoil_AC;

                % aerodynamic reference axis:
                % given by BlCrvAC, the local out-of-plane offset (when the
                % blade-pitch angle is zero) of the aerodynamic center normal
                % to the blade pitch axis (local blade one in blade root) and
                % by BlSwpAC, the local in-plane offset of AC to pitch axis
                arr_BlCrvAC  = arr_blade_values(:,2)'; % in flap-wise but normal to initial blade axis
                arr_BlSwpAC  = arr_blade_values(:,3)'; % in edge-wise but normal to initial blade axis
                arr_BlCrvAng = arr_blade_values(:,4)'; % specifies the local angle (in degrees) from the blade-pitch axis of a vector normal to the plane of the airfoil

    %             figure; hold on;
    %             i1 = [0;0;1];
    %             for k = 1:length(grid)-1
    %                 c1 = [arr_BlCrvAC(k),arr_BlSwpAC(k),arr_span(k)]';
    %                 c2 = [arr_BlCrvAC(k+1),arr_BlSwpAC(k+1),arr_span(k+1)]';
    %                 coord = c2-c1; 
    %                 e = coord/norm(coord);
    %                 alpha(k,1) = sign(coord(1))*acos(e'*i1)*180/pi;
    %             end
    %             plot(arr_span,[alpha;alpha(end)],'-k')
    %             plot(arr_span,arr_BlCrvAng,'-r')

                model.components.blade.DeSiO.uvlm.airfoil_orientation.grid   = grid;
                model.components.blade.DeSiO.uvlm.airfoil_orientation.values = arr_BlCrvAng;

%                 model.components.blade.DeSiO.beam.airfoil_orientation.grid   = grid;
%                 model.components.blade.DeSiO.beam.airfoil_orientation.values = arr_BlCrvAng;

                model.components.blade.DeSiO.uvlm.reference_axis.x.grid   = grid;
                model.components.blade.DeSiO.uvlm.reference_axis.x.values = arr_BlCrvAC;

                model.components.blade.DeSiO.uvlm.reference_axis.y.grid   = grid;
                model.components.blade.DeSiO.uvlm.reference_axis.y.values = arr_BlSwpAC;

                model.components.blade.DeSiO.uvlm.reference_axis.z.grid   = grid;
                model.components.blade.DeSiO.uvlm.reference_axis.z.values = arr_span;

                % set twist for beam - to have aerogrid and fe mesh same twist
    %             model.components.blade.DeSiO.beam.twist = [model.components.blade.DeSiO.uvlm.twist];
            end
            TF = contains(tline,'NumBlNds','IgnoreCase',true);
            if TF ~= 0
                str_cell_ne = textscan(tline,'%10.3f, %s');
                ne          = str_cell_ne{1};
                grid        = [0:1/(ne-1):1];
            end
        end
        fclose(fid);
    else
        disp('Warning: no AeroDyn_Blade.dat file found!');
    end
end

% open ElastDyn_Blade .dat to get pitch_axis, twist angle: only used if ElastoDyn is
% switched on
if ~isempty(str_wt_ElastoDynBlade)
    str_cell_ne = {}; ne = 0; grid = [];
    fid = fopen([str_wt_ElastoDynBlade],'r');
    if fid > 0
        disp(['reading ElastoDyn_Blade information from ' str_wt_ElastoDynBlade]);
        disp('   using ElastoDyn_Blade input - overwritting pitch axis in uvlm model');
        disp('   using ElastoDyn_Blade input - assigning reference axis fro uvlm model to beam model');
        while ~feof(fid)
            tline = fgetl(fid);
            TF1 = contains(tline,'BlFract','IgnoreCase',true);
            TF2 = contains(tline,'PitchAxis','IgnoreCase',true);
            TF3 = contains(tline,'StrcTwst','IgnoreCase',true);
            if (TF1 ~= 0 && TF2 ~= 0 && TF3 ~= 0)
                tline = fgetl(fid);
                for j = 1:ne
                    tline = fgetl(fid);
                    arr_values(j,:) = str2num(tline);
                end
%                 model.components.blade.DeSiO.uvlm.pitch_axis.grid = grid;
%                 model.components.blade.DeSiO.uvlm.pitch_axis.values = arr_values(:,2)';
                
                model.components.blade.DeSiO.beam.twist.grid = grid;
                model.components.blade.DeSiO.beam.twist.values = arr_values(:,3)';
                
                model.components.blade.DeSiO.beam.reference_axis = model.components.blade.DeSiO.uvlm.reference_axis;
            end
            TF = contains(tline,'NBlInpSt','IgnoreCase',true);
            if TF ~= 0
                str_cell_ne = textscan(tline,'%10.3f, %s');
                ne          = str_cell_ne{1};
                grid        = [0:1/(ne-1):1];
            end
        end
        fclose(fid);
    else
        disp('Warning: no ElastoDyn_Blade.dat file found!');
    end
end

% open BeamDyn.dat to get reference axis and initial twist for beam model
% and airfoil orientation but this from AeroDyn file!
str_cell_ne = {}; ne = 0; grid = [];
if ~isempty(str_wt_BeamDyn)
    fid = fopen([str_wt_BeamDyn],'r');
    if fid > 0
        disp(['reading BeamDyn information from ' str_wt_BeamDyn]);
        while ~feof(fid)
            tline = fgetl(fid);
            TF1 = contains(tline,'kp_xr','IgnoreCase',true);
            TF2 = contains(tline,'kp_yr','IgnoreCase',true);
            TF3 = contains(tline,'kp_zr','IgnoreCase',true);
            TF4 = contains(tline,'initial_twist','IgnoreCase',true);
            if (TF1 ~= 0 && TF2 ~= 0 && TF3 ~= 0 && TF4 ~= 0)
                tline = fgetl(fid);
                for j = 1:ne
                    tline = fgetl(fid);
                    arr_coord(j,:) = str2num(tline);
                end
                % set grid for reference axis
                for j = 1:ne-1
                    ds(j) = norm(arr_coord(j+1,1:3)-arr_coord(j,1:3));
                end
                res_ref = [0,cumsum(ds)];
                grid = res_ref./max(res_ref);
                model.components.blade.DeSiO.beam.reference_axis.x.grid = grid;
                model.components.blade.DeSiO.beam.reference_axis.y.grid = grid;
                model.components.blade.DeSiO.beam.reference_axis.z.grid = grid;

                model.components.blade.DeSiO.beam.reference_axis.x.values = arr_coord(:,1)';
                model.components.blade.DeSiO.beam.reference_axis.y.values = arr_coord(:,2)';
                model.components.blade.DeSiO.beam.reference_axis.z.values = arr_coord(:,3)';

                % set beam twist angle
                model.components.blade.DeSiO.beam.twist.grid = grid;
                model.components.blade.DeSiO.beam.twist.values = arr_coord(:,4)'.*pi/180;

            end
            TF = contains(tline,'Total number of key points','IgnoreCase',true);
            if TF ~= 0
                str_cell_ne = textscan(tline,'%10.3f, %s');
                ne          = str_cell_ne{1};
                grid        = [0:1/(ne-1):1];
            end
        end
        fclose(fid);
    else
        disp('Warning: no BeamDyn.dat file found!');
    end
end

% open ElastDyn for rotor informations, uptilt, precone, pitch, ...
arr_blade_pitch = [];
if ~isempty(str_wt_ElastoDyn)
    fid = fopen([str_wt_ElastoDyn],'r');
    if fid > 0
        disp(['reading ElastoDyn information from ' str_wt_ElastoDyn]);
        while ~feof(fid)
            tline = fgetl(fid);
            % pitch angle(s)
            TF = contains(tline,'BlPitch','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                arr_blade_pitch(end+1) = cell_str{1};
                model.environment.pitch_angle = [arr_blade_pitch(end)];
            end
            % number of blades
            TF = contains(tline,'NumBl','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.assembly.temp = [0,0,0];
                model.assembly.number_of_blades = [cell_str{1}];
            end
            % rotor diameter
            TF = contains(tline,'TipRad','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.assembly.rotor_diameter = [2*cell_str{1}];
            end
            % hub diameter
            TF = contains(tline,'HubRad','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.hub.DeSiO.temp = zeros(3,1);
                model.components.hub.DeSiO.diameter = [2*cell_str{1}];
            end
            % PreCone
            TF = contains(tline,'PreCone','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.hub.DeSiO.cone_angle = [cell_str{1}*pi/180]; % in rad
            end
            % Distance from rotor apex to hub mass (positiv upwind)
            TF = contains(tline,'HubCM','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.hub.DeSiO.Hub2Apex = [cell_str{1}];
            end
            % OverHang
            TF = contains(tline,'OverHang','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.nacelle.DeSiO.drivetrain.temp = zeros(3,1);
                model.components.nacelle.DeSiO.drivetrain.overhang = [-1.0*cell_str{1}];
            end
            % tilt angle
            TF = contains(tline,'ShftTilt','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.nacelle.DeSiO.drivetrain.uptilt = [-1.0*cell_str{1}*pi/180]; % in rad
            end
            % Twr2Shft
            TF = contains(tline,'Twr2Shft','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.nacelle.DeSiO.drivetrain.Twr2Shft = [cell_str{1}];
            end
            % TowerHt
            TF = contains(tline,'TowerHt','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.tower.DeSiO.towerHT = [cell_str{1}];
            end
            % center of mass of nacelle in tt-coordinates
            TF = contains(tline,'NacCMxn','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.nacelle.DeSiO.elastic_properties_mb.center_mass(1) = [cell_str{1}];
            end
            TF = contains(tline,'NacCMyn','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.nacelle.DeSiO.elastic_properties_mb.center_mass(2) = [cell_str{1}];
            end        
            TF = contains(tline,'NacCMzn','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.nacelle.DeSiO.elastic_properties_mb.center_mass(3) = [cell_str{1}];
            end         
            % yaw_mass
            TF = contains(tline,'YawBrMass','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.nacelle.DeSiO.drivetrain.yaw_mass = [cell_str{1}];
            end
            % platform mass
            TF = contains(tline,'PtfmMass','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                model.components.monopile.temp = zeros(3,1);
                model.components.monopile.transition_piece_mass = [cell_str{1}];
            end        
        end
        fclose(fid);
    else
        disp('Warning: no ElastoDyn.dat file found!');
    end
end

% open ElastDyn_tower for tower informations
if ~isempty(str_wt_ElastoDyn_tow)
    fid = fopen([str_wt_ElastoDyn_tow],'r');
    if fid >= 0
        disp(['reading ElastoDyn_Tower information from ' str_wt_ElastoDyn_tow]);
        model.components.tower.DeSiO.uvlm.M_aero  = 50;
        model.components.tower.DeSiO.uvlm.N_aero  = 10;
        model.components.tower.DeSiO.beam.M_struc = 50;
        model.components.tower.DeSiO.beam.dissipation.alpha_s = 0.2;
        model.components.tower.DeSiO.beam.dissipation.alpha_v = 0.2;
        while ~feof(fid)
            tline = fgetl(fid);
            TF = contains(tline,'NTwInpSt','IgnoreCase',true);
            if (TF ~= 0)
                cell_str = textscan(tline,'%10.3f, %s');
                number_tower_point = cell_str{1};
            end
            TF1 = contains(tline,'HtFract','IgnoreCase',true);
            TF2 = contains(tline,'TMassDen','IgnoreCase',true);
            TF3 = contains(tline,'TwFAStif','IgnoreCase',true);
            TF4 = contains(tline,'TwSSStif','IgnoreCase',true);
            if (TF1~=0 && TF2~=0 && TF3~=0 && TF4~=0)
                tline = fgetl(fid);
                for j = 1:number_tower_point
                    tline = fgetl(fid);
                    arr_tower_inf(j,:) = str2num(tline);
                    grid_tower(j) = arr_tower_inf(j,1);
                    x_coord_tower(j) = 0;
                    y_coord_tower(j) = 0;
                    if isfield(model.components,'tower') 
                        tower_height = model.components.tower.DeSiO.towerHT;
                    else
                        tower_height = 0;
                        disp('Warning, no tower defined');
                    end
                    z_coord_tower(j) = grid_tower(j)*model.components.tower.DeSiO.towerHT;
                    % assumed material: steel
                    [da(j),t(j),rho,E] = fun_get_diamter_thickness_from_EA_EI(arr_tower_inf(j,2),arr_tower_inf(j,3));
                end
                model.components.tower.DeSiO.uvlm.reference_axis.x.grid = grid_tower;
                model.components.tower.DeSiO.uvlm.reference_axis.y.grid = grid_tower;
                model.components.tower.DeSiO.uvlm.reference_axis.z.grid = grid_tower;

                model.components.tower.DeSiO.uvlm.reference_axis.x.values = x_coord_tower;
                model.components.tower.DeSiO.uvlm.reference_axis.y.values = y_coord_tower;
                model.components.tower.DeSiO.uvlm.reference_axis.z.values = z_coord_tower;

                model.components.tower.DeSiO.uvlm.outer_diameter.grid   = grid_tower;
                model.components.tower.DeSiO.uvlm.outer_diameter.values = da;
                model.components.tower.DeSiO.uvlm.thickness.grid        = grid_tower;
                model.components.tower.DeSiO.uvlm.thickness.values      = t;

                model.components.tower.DeSiO.beam.reference_axis  = model.components.tower.DeSiO.uvlm.reference_axis;
                model.components.tower.DeSiO.beam.outer_diameter  = model.components.tower.DeSiO.uvlm.outer_diameter;
                model.components.tower.DeSiO.beam.thickness       = model.components.tower.DeSiO.uvlm.thickness;
                model.components.tower.DeSiO.beam.material        = 'steel';
                model.components.tower.DeSiO.beam.shear_factor.k1 = 0.506;
                model.components.tower.DeSiO.beam.shear_factor.k2 = 0.506;
            end
        end
        fclose(fid);
    else
        disp('Warning: no ElastoDyn_tower.dat input file found');
    end
end

% set airfoil to model
if ~isempty(airfoil_data)
    model.airfoils = airfoil_data;
end

cd(currDir);
end

function [da,t,rho,E] = fun_get_diamter_thickness_from_EA_EI(rhoA,EI)
    % assumed material: steel
    rho = 7850; % kg/m^3
    E   = 2.1000e+11; % N/m^2

    da =[
      -((E*rhoA^2 - (16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + 4*pi*EI*rho^2)*(((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2))/(E*rhoA^2*pi^(1/2))
     -(((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + E*rhoA^2 + 4*pi*EI*rho^2)*(-((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) - 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2))/(E*rhoA^2*pi^(1/2))
       ((E*rhoA^2 - (16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + 4*pi*EI*rho^2)*(((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2))/(E*rhoA^2*pi^(1/2))
      (((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + E*rhoA^2 + 4*pi*EI*rho^2)*(-((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) - 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2))/(E*rhoA^2*pi^(1/2))];
    t = [
     -(((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2)/pi^(1/2)
     -(-((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) - 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2)/pi^(1/2)
       (((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) + 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2)/pi^(1/2)
      (-((16*pi^2*EI^2*rho^4 - E^2*rhoA^4)^(1/2) - 4*pi*EI*rho^2)/(E*rho*rhoA))^(1/2)/pi^(1/2)];

    da = max(abs(da));
    t  = min(abs(t));
    
    % check:
    A = pi/4*(da^2-(da-2*t)^2); I = pi/64*(da^4-(da-2*t)^4);
    rhoA_cal = rho*A; EI_cal = E*I;
end
