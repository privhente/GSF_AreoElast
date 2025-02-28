% =================================================================================================================
% Function to extract geometry from WindIO for structual mesh in DeSiO-Format.
% 
% Author: Christian Hente
% Date: 05.05.2022
%
% input:
%   strName - type of surface of (airfoil) cross-section
%   model - component beam object specified in WindIO
%   materials - material data specified in WindIO
%   scale_opt - scaling factor for scaling in longitudinal direction (optional input)
% output:
%     beam - beam object containing coordinates and connectivity for creating structural mesh
% =================================================================================================================
function beam = fun_get_beam_model(strName,model,materials,scale_opt,mesh_type_opt,mesh_factor_opt)
% =================================================================================================================
    scale       = 1;
    mesh_type   = 0;
    mesh_factor = 0;

    if nargin >= 4; scale = scale_opt; end
    if nargin >= 5; mesh_type = mesh_type_opt; end
    if nargin >= 6; mesh_factor = mesh_factor_opt; end
    
    % extracting beam data from WindIO
    beam  = fun_extract_beam_data(strName,model,materials,mesh_type,mesh_factor);
    
    % twist angle around pitch axis
    arr_twist = ones(beam.M+1,1)*0;
    if isfield(model,'twist')
       arr_twist = -beam.arr_twist;
    end
        
    % reference axis
    arr_xre   = [beam.arr_xre_x,beam.arr_xre_y,beam.arr_xre_z*scale];
   
    % airfoil orientation
    arr_airfoil_orientation = [];
    if isfield(beam,'airfoil_orientation')
        arr_airfoil_orientation = beam.airfoil_orientation;
    end   

    % set global (local) blabe basis
    
    i1 = [1;0;0];
    i2 = [0;1;0];
    i3 = [0;0;1];

    % beam position and director coordinates of the reference axis
    beam.arr_coordinates = [];
    for i = 1:model.M_struc+1
        
        % save the old directors from the previous segment
        if i == 1
            d1_old = [1; 0; 0];
            d2_old = [0; 1; 0];
            d3_old = [0; 0; 1];
            alpha_old = 0;
        else
            d1_old = d1';   % transposed because they were transposed for the output
            d2_old = d2';
            d3_old = d3';
            alpha_old = arr_twist(i-1);
        end
        
        % compute the new 3-director
        if ~isempty(arr_airfoil_orientation)
            R_axis = cos(arr_airfoil_orientation(i)*pi/180)*eye(3) + sin(arr_airfoil_orientation(i)*pi/180)*skew(i2)+(1-cos(arr_airfoil_orientation(i)*pi/180))*i2*i2';
            n3 = R_axis*i3;
        else
            n3 = zeros(3, 1);   % the current connection vector
            if i == model.M_struc+1
                n3(:) = arr_xre(i,:) - arr_xre(i-1,:);  % exeption at last node: connection interpolated backwards
            else
                n3(:) = arr_xre(i+1,:) - arr_xre(i,:);  % usual case: connection interpolated forward
            end
        end
        d3 = n3/norm(n3);   % 3-director as normed connection vector
        
%         % compute the new 3-director
%         n3 = zeros(3, 1);   % the current connection vector
%         if i == model.M_struc+1
%             n3(:) = arr_xre(i,:) - arr_xre(i-1,:);  % exeption at last node: connection interpolated backwards
%         else
%             n3(:) = arr_xre(i+1,:) - arr_xre(i,:);  % usual case: connection interpolated forward
%         end
%         d3 = n3/norm(n3);   % 3-director as normed connection vector
        
        % compute the angle between old and new 3-director and the
        % corresponding rotation axis
        theta = acos(dot(d3_old, d3));   % rotation angle
        
        % compute the rotation axis for the computation with the new
        % 3-director if the rotation angle is not very small
        if theta <= 1e-6   % special case (for the last segment): the old and new director are very similar --> then just use one of the old directors as 
            rot_matrix = eye(3);
        else
            rot_axis = cross(d3_old, d3)/norm(cross(d3_old, d3));   % direction vector of the rotation axis
            u_x = rot_axis(1);  % 1-component of the rotation axis
            u_y = rot_axis(2);  % 2-component of the rotation axis
            u_z = rot_axis(3);  % 3-component of the rotation axis
            
        % compute the rotation matrix
            rot_matrix = [
            cos(theta) + u_x^2*(1-cos(theta)),          u_x*u_y*(1-cos(theta)) - u_z*sin(theta),	u_x*u_z*(1-cos(theta)) + u_y*sin(theta);
            u_y*u_x*(1-cos(theta)) + u_z*sin(theta),	cos(theta) + u_y^2*(1-cos(theta)),          u_y*u_z*(1-cos(theta)) - u_x*sin(theta);
            u_z*u_x*(1-cos(theta)) - u_y*sin(theta),    u_z*u_y*(1-cos(theta)) + u_x*sin(theta),    cos(theta) + u_z^2*(1-cos(theta))];
        end
        
        % rotate the directors to get the new COS (without twist)
        d3_check = mtimes(rot_matrix, d3_old);
        d3_diff = d3 - d3_check;
        if norm(d3_diff)>1e-10
            disp('ERROR: rotation does not work as intended');
            disp(d3_old);
            disp(d3);
            disp(d3_check);
            return
        end
        d2_noTw = mtimes(rot_matrix, d2_old);
        d1_noTw = mtimes(rot_matrix, d1_old);
        
        % include the twist
        alpha_tw = arr_twist(i) - alpha_old;   % current twist angle
        
        % rotation axis for the twist rotation: the current 3-director
        rot_axis = d3;
        u_x = rot_axis(1);  % 1-component of the rotation axis
        u_y = rot_axis(2);  % 2-component of the rotation axis
        u_z = rot_axis(3);  % 3-component of the rotation axis
        
        % compute the rotation matrix for the twist rotation
        rot_matrix_Tw = [
        cos(alpha_tw) + u_x^2*(1-cos(alpha_tw)),          u_x*u_y*(1-cos(alpha_tw)) - u_z*sin(alpha_tw),	u_x*u_z*(1-cos(alpha_tw)) + u_y*sin(alpha_tw);
        u_y*u_x*(1-cos(alpha_tw)) + u_z*sin(alpha_tw),	cos(alpha_tw) + u_y^2*(1-cos(alpha_tw)),          u_y*u_z*(1-cos(alpha_tw)) - u_x*sin(alpha_tw);
        u_z*u_x*(1-cos(alpha_tw)) - u_y*sin(alpha_tw),    u_z*u_y*(1-cos(alpha_tw)) + u_x*sin(alpha_tw),    cos(alpha_tw) + u_z^2*(1-cos(alpha_tw))];
        
        % compute the new directors
        d1 = mtimes(rot_matrix_Tw, d1_noTw);
        d2 = mtimes(rot_matrix_Tw, d2_noTw);
        
        % transpose the vectors for output
        d1 = d1';
        d2 = d2';
        d3 = d3';
        
        % write the computed current director and reference point to the output-array
        beam.arr_coordinates(i,1:12) = [arr_xre(i,1:3),d1,d2,d3];
    end

    % beam connectivities
    for i = 1:model.M_struc
        beam.connectivity(i,1:3) = [i,i+1,i];
    end
% =================================================================================================================
return