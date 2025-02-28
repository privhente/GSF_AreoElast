% =================================================================================================================
% Function for writing DeSiO-Aero input files
% 
% Author: Christian Hente
% Date: 05.05.2022
% 
% input:
%   simu - struc variable that contains DeSiO-FSI simulation settings, i.e.
%   jobname, search radius
%   simu_struc - settings for structural simulation
%   simu_aero - settings for aerodynamic simulation
% =================================================================================================================
function fun_writeDeSiOFSIInput(simu,simu_struc,simu_aero,simu_timestep)
% =================================================================================================================
currDir = cd;

% creating uvlm input files and directories
strfilename = 'none';
if isfield(simu,'strfilename'); strfilename = simu.strfilename; end
if isfield(simu,'currDir'); currDir = simu.currDir; end
if isfield(simu,'caseDir'); caseDir = simu.caseDir; end
if strncmp(strfilename,'none',4)
    strfilename = 'solution';
else
    caseDir = [currDir '\' strfilename '\DeSiO'];
    mkdir(caseDir);
end

% if exist([currDir '\' 'DeSiO.bat'], 'file')
%    copyfile([currDir '\' 'DeSiO.bat'],[caseDir]);
% else
%   warningMessage = sprintf('Warning: file does not exist:\n%s', [currDir '\' 'DeSiO-FSI.bat']);
% end

cd(caseDir);
fid = fopen('simulationinput_fsi.txt','w');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');

    structprevjobname = 'none'; if isfield(simu,'structprevjobname'); structprevjobname = simu.structprevjobname; end
    aeroprevjobname   = 'none'; if isfield(simu,'aeroprevjobname');   aeroprevjobname   = simu.aeroprevjobname; end

    fprintf(fid,'!! filename (1), prev filename aero (2), prev filename structure (3)\n');
    fprintf(fid,'%s\t%s\t%s\n',strfilename,aeroprevjobname,structprevjobname);
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! total time (1) deltat aero (2) cutoff (3) time for load factor (4) type of load factor function (5)\n');
    
    lf_duration = 0.0;
    if isfield(simu,'lf_duration'); lf_duration = simu.lf_duration; end
    
    lf_type = 'constant';
    if isfield(simu,'lf_type'); lf_type = simu.lf_type; end
    
    time_aero = 1.0; 
    if isfield(simu_aero,'time'); time_aero = simu_aero.time; end
    
    deltat_aero = 0.1; 
    if isfield(simu_aero,'deltat'); deltat_aero = simu_aero.deltat; end
    
    cutoff = 1.0/100; 
    if isfield(simu_aero,'cutoff'); cutoff = simu_aero.cutoff; end
    
    str = strrep(sprintf('%10.5e\t',[time_aero,deltat_aero,cutoff,lf_duration]),'e','d');
    fprintf(fid,'%s\t%s',str,lf_type); fprintf(fid,'\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! wind data: sort(1), fluid density(2), intensity(3), duration(4), direction vector (5-7)\n');

    density = 1.0; i_vinf  = 1.0; time = 1.0; d_vinf = [1,0,0]; str_sort = 'constant';
    if isfield(simu_aero,'density'); density = simu_aero.density; end
    if isfield(simu_aero,'i_vinf'); i_vinf = simu_aero.i_vinf; end
    if isfield(simu_aero,'time'); time = simu_aero.time; end
    if isfield(simu_aero,'d_vinf'); d_vinf = simu_aero.d_vinf; end
    if isfield(simu_aero,'sort'); str_sort = simu_aero.sort; end
    
    str_settings = strrep(sprintf('%10.8e\t',[density i_vinf time d_vinf]),'e','d');
    fprintf(fid,'%s\t%s\n',str_sort,str_settings);
    
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! deltat (1), tolerance (2), iteration limit (3), gravity flag (4), row2: flag for writing matrices (1) - on/off = 1/0, aero tangent (2) - on/off = 1/0, flag for automatic aero time stepping (3) - on/off = 1/0 \n');
    
    struc_simType = 'dynamic';
    if isfield(simu,'simType')
        struc_simType = simu.simType;
    elseif isfield(simu_struc,'simType')
        struc_simType = simu_struc.simType;
    end
    
    str_fsi_type = 'strong';
    if isfield(simu,'fsi_type')
        str_fsi_type = simu.fsi_type;
    end
    
    if strncmp(struc_simType,'kinematic',length('kinematic'))
        str_fsi_type = '';
    end
        
    deltat_struc = deltat_aero;
    if isfield(simu_struc,'deltat'); deltat_struc = simu_struc.deltat; end;
    
    tol_struc = 1e-6;
    if isfield(simu_struc,'tol'); tol_struc = simu_struc.tol; end;
    
    niter_struc = 1e-6;
    if isfield(simu_struc,'niter'); niter_struc = simu_struc.niter; end;
    
    grav = 0;
    if isfield(simu_struc,'grav'); grav = simu_struc.grav; end;
    
    fprintf(fid,'%s %s\n',struc_simType,str_fsi_type);
    str = strrep(sprintf('%10.5e\t',[deltat_struc,tol_struc]),'e','d');
    fprintf(fid,'%s\t%i\t%i',str,niter_struc,grav); fprintf(fid,'\n');
    
    flag_linearization = 0;
    if isfield(simu,'flag_linearization'); flag_linearization = simu.flag_linearization; end
    
    flag_output_matrix = 0;
    if isfield(simu,'flag_linearization'); flag_output_matrix = simu.flag_output_matrix; end
    fprintf(fid,'%i\t%i\t%i',flag_output_matrix,flag_linearization);fprintf(fid,'\n');

    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! gravity vector (1, 2, 3)\n');
    
    grav_vec = [0,0,0];
    if isfield(simu_struc,'grav_vec'); grav_vec = simu_struc.grav_vec; end;

    str = strrep(sprintf('%10.5e\t',grav_vec),'e','d');
    fprintf(fid,'%s',str); fprintf(fid,'\n');
    
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! type adaptive time step (1), output flag (2), deviation parameter (3), increasing percentage (4), initial step to activate adaptive time step (5), relation initial wake length - max wake length (6), , vector of blade tip nodes (7,8,9)\n');
    
    
    type_adaptive_timestep = 'off';
    if isfield(simu_timestep,'type_adaptive_timestep')
        type_adaptive_timestep = simu_timestep.type_adaptive_timestep; 
    end
    fprintf(fid,'%s \n',type_adaptive_timestep);
    
    flag_auto_aero_stepping_output = 0;
    if isfield(simu_timestep,'flag_auto_aero_stepping_output'); flag_auto_aero_stepping_output = simu_timestep.flag_auto_aero_stepping_output; end
    
    adaptive_timestep_deviation = 0;
    if isfield(simu_timestep,'adaptive_timestep_deviation'); adaptive_timestep_deviation = simu_timestep.adaptive_timestep_deviation; end
    
    adaptive_timestep_initialstep = 0;
    if isfield(simu_timestep,'adaptive_timestep_initialstep'); adaptive_timestep_initialstep = simu_timestep.adaptive_timestep_initialstep; end
    
    adaptive_timestep_increase_percentage = 0;
    if isfield(simu_timestep,'adaptive_timestep_increase_percentage'); adaptive_timestep_increase_percentage = simu_timestep.adaptive_timestep_increase_percentage; end
    
    adaptive_timestep_rel_wake_length = 0;
    if isfield(simu_timestep,'adaptive_timestep_rel_wake_length'); adaptive_timestep_rel_wake_length = simu_timestep.adaptive_timestep_rel_wake_length; end
    
    fprintf(fid,'%i\t%i\t%i\t%i\t%i\t\n',flag_auto_aero_stepping_output,adaptive_timestep_deviation,adaptive_timestep_increase_percentage,adaptive_timestep_initialstep,adaptive_timestep_rel_wake_length);

    adaptive_timestep_blade_tip_nodes = [0,0,0];
    if isfield(simu_timestep,'adaptive_timestep_blade_tip_nodes'); adaptive_timestep_blade_tip_nodes = simu_timestep.adaptive_timestep_blade_tip_nodes; end
    
    integer_values = int64(adaptive_timestep_blade_tip_nodes);
    str = sprintf('%-4d\t%-4d\t%-4d\n', integer_values);
    fprintf(fid, '%s\n', str);
    
    if strncmp(str_sort,'file',4)
        if isfield(simu_aero,'wind_field_file')
            if exist([currDir '\' simu_aero.wind_field_file], 'file')
              copyfile([currDir '\' simu_aero.wind_field_file],[caseDir]);
            else
              warningMessage = sprintf('Warning: Wind data file does not exist:\n%s', [currDir '\' simu_aero.wind_field_file]);
              disp(warningMessage)
            end
            fprintf(fid,'!!\n');
            fprintf(fid,'!! if sort = file then consider these addtional lines for inflow settings\n');
            fprintf(fid,'!! filename for wind field: filename with file extension (1)\n');
            
            % extension of file:
            [inz_Ext] = strfind(simu_aero.wind_field_file,'.');
            wind_field_file_ext = simu_aero.wind_field_file(inz_Ext(end)+1:end);

            if strncmp(wind_field_file_ext,'wnd',4)
                wind_field_type = 1;
            elseif strncmp(wind_field_file_ext,'bts',4)
                wind_field_type = 2;
            else
                wind_field_type = 1;
                disp(['Warning: wrong file type for wind_field_file! Needed (.wnd or .bts)']);    
            end 

            fprintf(fid,'%s\n',[simu_aero.wind_field_file]);
            fprintf(fid,'!!\n');
            fprintf(fid,'!! inflow settings\n');
            fprintf(fid,'!! grids center (1-3)\n');
            grid_center = [0,0,0];
            if isfield(simu_aero,'grid_center'); grid_center = simu_aero.grid_center; end
            str = strrep(sprintf('%10.8e\t',grid_center),'e','d');
            fprintf(fid,'%s\n',str);
        else
            disp(['Warning: no field: wind_field_file found!']);
        end
    end
    
fclose(fid);

fid = fopen('fsi_input.txt','w');
    nfsi = 0;
    if isfield(simu,'data'); nfsi = size(simu.data,1); end;
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! number of fsi (1)\n');
    fprintf(fid,'%i\n',nfsi);
    
    for i = 1:nfsi
        fprintf(fid,'!! input for fluid-structure interaction:\n');
        fprintf(fid,'!! (row 1) shell or beam (1), surface from fluid (2), number of structural bodies, (3) local search radius\n');
        fprintf(fid,'!! (row 2) bodies of structural elements\n');
        if ischar(simu.data{i,4})
            fprintf(fid,'%s\t%i\t%i\t%s\n',simu.data{i,1},simu.data{i,2},simu.data{i,3},simu.data{i,4});
            fprintf(fid,'%i\n',simu.data{i,5});
        else
            fprintf(fid,'%s\t%i\t%i\t%20.10fd0\n',simu.data{i,1},simu.data{i,2},simu.data{i,3},simu.data{i,4});
            fprintf(fid,'%i\t',simu.data{i,5});
            fprintf(fid,'\n');
        end
    end
fclose(fid);

% write fsi_radius_input_files
if isfield(simu,'data_input')
    if not(isempty(simu.data_input))
        for i = 1:nfsi
            fid = fopen([simu.data{i,4}],'w');
                fprintf(fid,'!! Input search radius for weights using radial-based function.\n');
                fprintf(fid,'!! fsi %i: %s %i to surface %i\n',i,simu.data{i,1},simu.data{i,2},simu.data{i,3});
                fprintf(fid,'!! number of nodes for fsi search radius\n');
                fprintf(fid,'%i\n',size(simu.data_input(i).node_fsi_radius,2));
                fprintf(fid,'!!\n');
                fprintf(fid,'!!\n');
                fprintf(fid,'!! fsi search radius (1)\n');
                str = strrep(sprintf('%10.5e\n',simu.data_input(i).node_fsi_radius),'e','d');
                fprintf(fid,'%s',str); fprintf(fid,'\n');
            fclose(fid);
        end
    end
end
% =================================================================================================================
cd(simu.currDir);
% disp('creating DeSiO-FSI input files');
% =================================================================================================================
return