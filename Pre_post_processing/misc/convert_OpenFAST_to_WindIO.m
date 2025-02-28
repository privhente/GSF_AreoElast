function convert_OpenFAST_to_WindIO()
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

% reading of old model
info.strYAMLFileName = 'IEA-15-240-RWT_new.yaml';

% define order of variable name in .yaml
str_names = {'name'
             'description'
             'jobname'
             'assembly'
             'simulationparameter'
             'environment'
             'components'
             'airfoils'
             'materials'};

% set model parameter:
model.name                                            = 'NREL 15 MW reference wind turbine';
model.description                                     = 'created from openfast input files using converter';
model.jobname                                         = '15mw';

model.assembly.turbine_class                          = 'I';
model.assembly.turbulence_class                       = 'B';
model.assembly.drivetrain                             = 'direct_drive';
model.assembly.rotor_orientation                      = 'Upwind';

model.environment.init_rotor_velocity                 = 7.25;
model.environment.vinf                                = 10.21;
%model.environment.wind_field_file                    = 'mean.wnd'
model.environment.pitch_angle                         = 0;
model.environment.yaw_angle                           = 0;
model.environment.winddir                             = [1.0,0.0,0.0];
model.environment.air_density                         = 1.225;

model.simulationparameter.os                          = 'windows';
model.simulationparameter.time                        = 10.0;
model.simulationparameter.deltat_aero                 = 0.10;
model.simulationparameter.deltat_struc                = 0.10;
model.simulationparameter.cutoff                      = 0.01;
model.simulationparameter.nwakerows                   = 200;
model.simulationparameter.wake_diffusion              = 0.00;
model.simulationparameter.tolerance                   = 1.0e-6;
model.simulationparameter.niter                       = 50;
model.simulationparameter.gravity                     = 1;
model.simulationparameter.grav_vector                 = [0,0,-9.81];
model.simulationparameter.flag_init_self_weight       = 0;
model.simulationparameter.flag_init_rotor_velocity    = 1;
model.simulationparameter.fsi_type                    = 'weak';
    
model.components.blade.DeSiO.uvlm.M_aero              = [100];
model.components.blade.DeSiO.uvlm.N_aero              = [10];
model.components.blade.DeSiO.uvlm.mesh_type           = [2];
model.components.blade.DeSiO.uvlm.mesh_factor         = [0.8];
model.components.blade.DeSiO.uvlm.blade_root_position = [8];
model.components.blade.DeSiO.beam.M_struc             = [model.components.blade.DeSiO.uvlm.M_aero];
model.components.blade.DeSiO.beam.mesh_type           = [2];
model.components.blade.DeSiO.beam.mesh_factor         = [0.8];
model.components.blade.DeSiO.beam.dissipation.temp    = [0,0];
model.components.blade.DeSiO.beam.dissipation.alpha_s = [0.8];
model.components.blade.DeSiO.beam.dissipation.alpha_v = [0.8];

% reading all files in director
info.wt_dir      = ['C:\Users\hente\Desktop\15_MW_OpenFast_DeSiO\OpenFAST Inputdateien\IEA-15-240-RWT-Monopile-Onshore-Operating\IEA-15-240-RWT'];
info.airfoil_dir = ['C:\Users\hente\Desktop\15_MW_OpenFast_DeSiO\OpenFAST Inputdateien\IEA-15-240-RWT-Monopile-Onshore-Operating\IEA-15-240-RWT\Airfoils'];

% openfast input filenames
% info.str_wt_ElastoDynBlade = 'IEA-15-240-RWT_ElastoDyn_blade.dat';
info.str_wt_BeamDyn       = 'IEA-15-240-RWT_BeamDyn.dat';
info.str_wt_AeroDynBlade  = 'IEA-15-240-RWT_AeroDyn15_blade.dat';
info.str_wt_AeroDyn       = 'IEA-15-240-RWT_AeroDyn15.dat';
info.str_wt_ElastoDyn     = 'IEA-15-240-RWT-Monopile_ElastoDyn.dat';
info.str_wt_ElastoDyn_tow = 'IEA-15-240-RWT-Monopile_ElastoDyn_tower.dat';

% get and convert model data
model = convert_rotor_to_WindIO(info, model);

% sorting model
str_fieldnames = fieldnames(model);
inz_names = [1:length(str_names)]; a = 0;
for i = 1:length(inz_names)
    for j = 1:size(str_fieldnames,1)
        if strcmp(str_fieldnames{j},str_names{i})
            a = a + 1;
            inz_new(a) = j;
        end
    end
end
table_model = struct2table(model,'AsArray',true);
table_model = table_model(1,inz_new);
model = table2struct(table_model);

% write into yaml-file
YAML.write(info.strYAMLFileName, model);
end