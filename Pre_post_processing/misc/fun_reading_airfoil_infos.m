% =================================================================================================================
function airfoil_data = fun_reading_airfoil_infos(airfoil_dir,airfoil_data,airfoil_ID)
% =================================================================================================================
    currDir = cd;
    % loop over all airfoil files to find the corresponding coord file
    cd(airfoil_dir);
    fid = fopen([airfoil_data(airfoil_ID).filename],'r');
    if fid >= 0
        j = 0; arr_coord = [];
        while ~feof(fid)
            tline = fgetl(fid);
            TF = contains(tline,'NumCoords','IgnoreCase',true);
            if TF ~= 0
                strline = textscan(tline,'%s %s',2);
                inz_tag = strfind(strline{1}{1},'"');
                airfoil_data(airfoil_ID).filename_coord = strline{1}{1}(inz_tag(1)+1:inz_tag(2)-1);
                [airfoil_data] = fun_read_airfoil_coord(airfoil_data(airfoil_ID).filename_coord,airfoil_data,airfoil_ID);
                break;
            end
        end
    fclose(fid);
    cd(currDir);
    else
        disp(['Error: file: ' airfoil_data(airfoil_ID).filename ' not found!'])
        return;
    end
end

% =================================================================================================================
function [airfoil] = fun_read_airfoil_coord(filename_coord,airfoil,a)
% =================================================================================================================
% loop over all files to write a new .yaml file only with the airfoils:
fid = fopen([filename_coord],'r');
    j = 0; arr_coord = []; flag_numcoord = 0;
    while ~feof(fid)
        tline = fgetl(fid);
        TF = contains(tline,'coordinates of airfoil shape','IgnoreCase',true);
        if TF ~= 0
            tline = fgetl(fid); tline = fgetl(fid);
            for k = 1:numCoords
%             while ~feof(fid)
                tline = fgetl(fid);
                strline = textscan(tline,'%10.5f %10.5f',2);
                arr_coord(k,1:2) = str2num(tline);
            end
        end
        TF = contains(tline,'x-y coordinate of airfoil reference','IgnoreCase',true);
        if TF ~= 0
            tline = fgetl(fid); tline = fgetl(fid);
            aerodynamic_center = str2num(tline);
        end
        TF = contains(tline,'NumCoords','IgnoreCase',true);
        if TF ~= 0
            if flag_numcoord == 0
                str_cell_ne = textscan(tline,'%10.5f %s');
                numCoords = str_cell_ne{1};
                flag_numcoord = 1;
            end
        end
    end
%     inz_coord = strfind(upper(filename_coord),'_COORD');
%     airfoil(a).name = filename_coord(1:inz_coord-1);
    airfoil(a).coordinates.x = arr_coord(:,1)';
    airfoil(a).coordinates.y = arr_coord(:,2)';
    airfoil(a).relative_thickness = 1.0;
    airfoil(a).aerodynamic_center = aerodynamic_center(1);
fclose(fid);
end