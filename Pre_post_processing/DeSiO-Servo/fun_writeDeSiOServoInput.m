% =================================================================================================================
% Function for writing DeSiO-Aero input files
% 
% Author: Christian Hente
% Date: 05.05.2022
% 
% input:
%   simu - struc variable that contains DeSiO model components, i.e. rb,
%   pointmass12,constraints,simulation settings
%   mesh - struc variable that containts grid coordinates and
%   connectivities of structural finite beam mesh
% =================================================================================================================
function fun_writeDeSiOServoInput(simu)
% =================================================================================================================
currDir = cd;
if ~isempty(simu)
    % creating structure input files and directories
    if isfield(simu,'jobname')
        caseDir = [simu.currDir '\' simu.strfilename '\DeSiO\' simu.jobname];
    else
        caseDir = [simu.currDir '\' simu.strfilename '\DeSiO\'];
    end
    mkdir(caseDir);

    if isfield(simu,'controller_inputfile')
        if exist([currDir '\' simu.controller_inputfile], 'file')
          copyfile([currDir '\' simu.controller_inputfile],[caseDir]);
        else
          warningMessage = sprintf('Warning: file does not exist:\n%s', [currDir '\' simu.controller_inputfile]);
        end
    end
    
    if isfield(simu,'controller_rotor_performance')
        if exist([currDir '\' simu.controller_rotor_performance '.txt'], 'file')
          copyfile([currDir '\' simu.controller_rotor_performance '.txt'],[caseDir]);
        else
          warningMessage = sprintf('Warning: file does not exist:\n%s', [currDir '\' simu.controller_rotor_performance '.txt']);
        end
    end
    
    % writing input files
    cd(caseDir);
    fid = fopen('controllerinput.txt','w');
    
        fprintf(fid,'!! \n'); 
        fprintf(fid,'!! \n'); 
        fprintf(fid,'!! operating system: windows or linux \n');
        str_os = 'windows'; 
        if isfield(simu,'OS') str_os = simu.OS; end
        fprintf(fid,'%s\n',str_os);

        fprintf(fid,'!! \n'); 
        fprintf(fid,'!! \n');
        fprintf(fid,'!! controller input file name - e.g. DISCON_DATA.IN (case sensitive)\n');
        controller_inputfile = 'DISCON_DATA.IN';
        if isfield(simu,'controller_inputfile'); controller_inputfile = simu.controller_inputfile; end
        fprintf(fid,'%s\n',simu.controller_inputfile);

        fprintf(fid,'!! \n'); 
        fprintf(fid,'!! \n');
        fprintf(fid,'!! initial pitch angle \n');
        init_pitch = 0.0;
        if isfield(simu,'init_pitch'); init_pitch = simu.init_pitch; end
        str = strrep(sprintf('%10.8e',init_pitch),'e','d'); fprintf(fid,'%s \n',str);

        fprintf(fid,'!! \n'); 
        fprintf(fid,'!! \n');
        fprintf(fid,'!! number of blades \n');
        n_blades = 3;
        if isfield(simu,'n_blades'); n_blades = simu.n_blades; end
        fprintf(fid,'%i \n',n_blades);

        fprintf(fid,'!! row(1): beam body numbers \n'); 
        fprintf(fid,'!! row(2): revolutejoint_2 constraint numbers of \n');
        fprintf(fid,'!! row(2): rotation_local constraint numbers to blade pitching \n');
        blade_beam_bodies = [0,0,0]; if isfield(simu,'blade_beam_bodies'); blade_beam_bodies = simu.blade_beam_bodies; end
        constraint_revj2  = [0,0,0]; if isfield(simu,'constraint_revj2');  constraint_revj2  = simu.constraint_revj2; end
        constraint_rotloc = [0,0,0]; if isfield(simu,'constraint_rotloc'); constraint_rotloc = simu.constraint_rotloc; end
        fprintf(fid,'%i\t',blade_beam_bodies);  fprintf(fid,'\n');
        fprintf(fid,'%i\t',constraint_revj2);  fprintf(fid,'\n');
        fprintf(fid,'%i\t',constraint_rotloc); fprintf(fid,'\n');
    fclose(fid);
end

% =================================================================================================================
% close all;
cd(simu.currDir);
disp('creating DeSiO-Servo input files');
% =================================================================================================================
return