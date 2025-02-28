% =================================================================================================================
% Function to extract uvlm data from WindIO for aerodynamical grid in DeSiO-Format.
% 
% Author: Christian Hente
% Date: 05.05.2022
%
% input:
%   strName - type of surface of (airfoil) cross-section
%   model - component uvlm object specified in WindIO
%   airfoils - airfoil data specified in WindIO
% output:
%     uvlm_ob - uvlm object containing coordinates and connectivity for creating aerodynamic grid
% =================================================================================================================
function uvlm_obj = fun_extract_uvlm_data(strName,model_uvlm,airfoils,mesh_type_opt,mesh_factor_opt,blade_vl_model_opt)
% =================================================================================================================
    mesh_type   = 0;
    mesh_factor = 0;
    blade_vl_model = '2D';
    
    if nargin >= 4; mesh_type = mesh_type_opt; end
    if nargin >= 5; mesh_factor = mesh_factor_opt; end
    if nargin >= 6; blade_vl_model = blade_vl_model_opt; end
    
    % span-wise and chord-wise discretization
    uvlm_obj.M =  model_uvlm.M_aero; % span-wise
    uvlm_obj.N =  model_uvlm.N_aero; % chord-wise
    
    % natural coordinates of span-and chord-wise discretization
    xhi_x_li = [0:1/uvlm_obj.M:1];   % span-wise
    delta_xhi_c = 1/uvlm_obj.M;
    
    % cos distributed mesh size in span-wise direction:
    mc = delta_xhi_c;
    mf = mc*(1-mesh_factor+1e-6);
    if mesh_type == 1
        % type: coarse - fine:
        dxhi_x   = 0.5*(mc-mf)*(1+cos(pi*xhi_x_li)) + mf;
        xhi_x_nl = 0.5*(mc-mf)*(xhi_x_li+1/pi*sin(pi*xhi_x_li)) + mf*xhi_x_li;
    elseif mesh_type == 2
        % different type: fine - coarse - fine:
        dxhi_x   = mc + (mf-mc)*0.5*(1+cos(2*pi*xhi_x_li));
        xhi_x_nl = mc*xhi_x_li + 0.5*(mf-mc)*(xhi_x_li+1/(2*pi)*sin(2*pi*xhi_x_li)) ;%+ (1-0.5*(mc+mf));
    else
        dxhi_x   = ones(length(xhi_x_li)).*mc;
        xhi_x_nl = xhi_x_li;
    end
    xhi_x_nl = xhi_x_nl/max(xhi_x_nl);
%     figure; grid on; plot(xhi_x_li,dxhi_x);
%     figure; grid on; plot(xhi_x_li,xhi_x_nl);
    
    xhi_y_c  = [0:1/(uvlm_obj.N):1]; % chord-wise
    xhi_x    = xhi_x_nl;
    
    nbr             = 0;                    % number of blade root cross-section
    arr_inz_airfoil = [];                   % index array for locating airfoils in span-wise direction
    
    % interpolating coordinates of reference axis according to span-wise discretization
    uvlm_obj.arr_xre_x = interp1(model_uvlm.reference_axis.x.grid,model_uvlm.reference_axis.x.values,xhi_x,'linear')';
    uvlm_obj.arr_xre_y = interp1(model_uvlm.reference_axis.y.grid,model_uvlm.reference_axis.y.values,xhi_x,'linear')';
    uvlm_obj.arr_xre_z = interp1(model_uvlm.reference_axis.z.grid,model_uvlm.reference_axis.z.values,xhi_x,'linear')';
%     figure(); title('ref_x'); hold on; grid on; plot(xhi_x,uvlm_obj.arr_xre_x,'xr'); plot(model_uvlm.reference_axis.x.grid,model_uvlm.reference_axis.x.values,'-k');
%     figure(); title('ref_y'); hold on; grid on; plot(xhi_x,uvlm_obj.arr_xre_y,'xr'); plot(model_uvlm.reference_axis.y.grid,model_uvlm.reference_axis.y.values,'-k');
%     figure(); title('ref_z'); hold on; grid on; plot(xhi_x,uvlm_obj.arr_xre_z,'xr'); plot(model_uvlm.reference_axis.z.grid,model_uvlm.reference_axis.z.values,'-k');
    
    % airfoil orientation from AeroDyn-File
    if isfield(model_uvlm,'airfoil_orientation')
        uvlm_obj.airfoil_orientation = interp1(model_uvlm.airfoil_orientation.grid,model_uvlm.airfoil_orientation.values,xhi_x,'linear')';
%         figure(); title('airfoil_orientation'); hold on; grid on; plot(xhi_x,uvlm_obj.airfoil_orientation,'xr'); plot(model_uvlm.airfoil_orientation.grid,model_uvlm.airfoil_orientation.values,'-k');
    end
    
    % if-statement according to type of surface cross-section
    if strncmp(strName,'blade',length(strName))
        uvlm_obj.airfoil = model_uvlm.airfoil_position;
        nbr              = model_uvlm.blade_root_position;
        if nbr > size(uvlm_obj.airfoil.labels,1)
            nbr = size(uvlm_obj.airfoil.labels,1)
        end
            
        % interpolating chord length in span-wise direction
        uvlm_obj.arr_c   = interp1(model_uvlm.chord.grid,model_uvlm.chord.values,xhi_x,'linear')';
%         figure(); title('chord'); hold on; grid on; plot(xhi_x,uvlm_obj.arr_c,'xr'); plot(model_uvlm.chord.grid,model_uvlm.chord.values,'-k');
        
        % interpolating twist angle according to span-wise discretization
        uvlm_obj.arr_twist = ones(uvlm_obj.M+1,1)*0.0;
        if isfield(model_uvlm,'twist')
            uvlm_obj.arr_twist = interp1(model_uvlm.twist.grid,model_uvlm.twist.values,xhi_x,'linear')';
%             figure(); title('twist'); hold on; grid on; plot(xhi_x,uvlm_obj.arr_twist,'xr'); plot(model_uvlm.twist.grid,model_uvlm.twist.values,'-k');
        end
        % interpolating airfoil in span-wise direction. This is important
        % to identify the blade root and blade's lifting surfaces
        airfoil_grid    = uvlm_obj.airfoil.grid;
        airfoil_values  = [1:length(uvlm_obj.airfoil.grid)];
        arr_inz_airfoil = fix(interp1(airfoil_grid,airfoil_values,xhi_x,'linear'))';
%         figure(); title('airfoilID'); hold on; grid on; plot(xhi_x,arr_inz_airfoil,'xr'); plot(airfoil_grid,airfoil_values,'-k');

        % camber_line_factor
        camber_line_factor = 0.5;
        if isfield(model_uvlm,'camber_line_factor')
            if isfield(model_uvlm.camber_line_factor,'values')
                uvlm_obj.arr_camber_line_factor = interp1(model_uvlm.camber_line_factor.grid,model_uvlm.camber_line_factor.values,airfoil_grid,'linear')';
            else
                camber_line_factor = model_uvlm.camber_line_factor;
            end
        end
%       
        arr_pitch_ax = ones(uvlm_obj.M+1,1)*1.0;
        if isfield(model_uvlm,'pitch_axis')
            uvlm_obj.arr_pitch_ax = interp1(model_uvlm.pitch_axis.grid,model_uvlm.pitch_axis.values,xhi_x,'linear')';
%             figure(); title('pitch_ax'); hold on; grid on; plot(xhi_x,uvlm_obj.arr_pitch_ax,'xr'); plot(model_uvlm.pitch_axis.grid,model_uvlm.pitch_axis.values,'-k');            
        end
        % natural coordinates of chord-wise discretization for whole
        % surface of cross-section
%         if strncmp(blade_vl_model,'2D',2)
            xhi_y_w = xhi_y_c;
%         elseif strncmp(blade_vl_model,'3D',2)
%             xhi_y_w = 0.5*(1-cos([0:(pi)/(uvlm_obj.N):pi]));
%         end
        
    elseif strncmp(strName,'pipe',length(strName))
        % interpolating chord length in span-wise direction
        uvlm_obj.arr_c          = interp1(model_uvlm.outer_diameter.grid,model_uvlm.outer_diameter.values-model_uvlm.thickness.values,xhi_x,'linear');
        uvlm_obj.arr_twist      = ones(uvlm_obj.M+1,1)*0.0; % zero twist
        uvlm_obj.arr_pitch_ax   = ones(uvlm_obj.M+1,1)*0.5; % location of pitch axis
        uvlm_obj.airfoil.labels = {model_uvlm.cross_section;model_uvlm.cross_section};  % artificial "airfoil" sections for pipe
        uvlm_obj.airfoil.grid   = [0.0, 1.0];               % artificial "airfoil" sections grid for pipe
        
        % mesh discretization around circular cross-section
        xhi_y_w = 0.5*(1-cos([0:(pi)/(uvlm_obj.N):pi]));
        
        % camber_line_factor
        camber_line_factor = 0.5;
        if isfield(model_uvlm,'camber_line_factor')
            if isfield(model_uvlm.camber_line_factor,'values')
                uvlm_obj.arr_camber_line_factor = interp1(model_uvlm.camber_line_factor.grid,model_uvlm.camber_line_factor.values,airfoil_grid,'linear')';
            else
                camber_line_factor = model_uvlm.camber_line_factor;
            end
        end        
    end
    
    % get information of airfoil list:
    strinfo = whos('airfoils');
    strlegend = {};
    hpp = [];
    
    % loop over airfoils in the model to calculate coordinates of
    % cross-section surfaces
    for i = 1:size(uvlm_obj.airfoil.labels,1)
        airfoil_name = uvlm_obj.airfoil.labels(i);

        % searching current airfoil from airfoil-list
        if strcmp(strinfo.class,'cell')
            for j = 1:size(airfoils,1)
                if size(airfoils,1) == 1
                    if strcmp(airfoils.name,airfoil_name)
                        airfoilj = airfoils;
                        break
                    end
                else
                    if strcmp(airfoils{j}.name,airfoil_name)
                        airfoilj = airfoils{j};
                        break
                    end
                end
            end
        elseif strcmp(strinfo.class,'struct')
            for j = 1:size(airfoils,2)
                if size(airfoils,2) == 1
                    if strcmp(airfoils.name,airfoil_name)
                        airfoilj = airfoils;
                        break
                    end
                else
                    if strcmp(airfoils(j).name,airfoil_name)
                        airfoilj = airfoils(j);
                        break
                    end
                end
            end
        end
        % extracting coordinates for upper and lower airfoil 
        airfoil_coord = [airfoilj.coordinates.x',airfoilj.coordinates.y'];
        [val,inzmin] = min(airfoil_coord(2:end-1,1));
        inz0 = inzmin+1;
        airfoil_coord_u = airfoil_coord(1:inz0(1),:);   [val,inzsort] = sort(airfoil_coord_u(:,1)); airfoil_coord_u = airfoil_coord_u(inzsort,:);
        airfoil_coord_l = airfoil_coord(inz0(1):end,:); [val,inzsort] = sort(airfoil_coord_l(:,1)); airfoil_coord_l = airfoil_coord_l(inzsort,:);
        
        % interpolating values in airfoil thickness direction according to
        % discretization in chord-wise direction to calculate camber
        % surface coordinates
        xhi_airf_u = fun_lin_interpolation(airfoil_coord_u,xhi_y_c);
        xhi_airf_l = fun_lin_interpolation(airfoil_coord_l,xhi_y_c);
        
        % coordinates of camber surface
        if isfield(uvlm_obj,'arr_camber_line_factor')
            clf = uvlm_obj.arr_camber_line_factor(i);
        else
            clf = camber_line_factor;
        end
        xhi_airf_c = xhi_airf_u(:,1)*clf + xhi_airf_l(:,1)*(1-clf);
%         xhi_airf_c = (xhi_airf_u(:,1) + xhi_airf_l(:,1))/2;
%         xhi_airf_c = xhi_airf_u(:,1);
%         xhi_airf_c = [[0,xhi_airf_c(1,1)];[1,xhi_airf_c(end,1)]];
%         xhi_airf_c = fun_lin_interpolation(xhi_airf_c,xhi_y_c);
        
%         if i == 1
%             figure(); hold on; grid on; axis equal;
%             rgbcolor = RGB_Color;
%         end
        
        if i >= 2
            u = fun_lin_interpolation(airfoil_coord_u,[0:1/100:1]);
            l = fun_lin_interpolation(airfoil_coord_l,[0:1/100:1]);
            c = (u(:,1) + l(:,1))/2;
%             hpp(end+1) = plot(airfoil_coord_u(:,1),airfoil_coord_u(:,2),'-','color',rgbcolor(i,:));
            strlegend{end+1} = airfoilj.name;
%             plot(airfoil_coord_l(:,1),airfoil_coord_l(:,2),'-','color',rgbcolor(i,:));
%             plot(xhi_y_c,xhi_airf_c,'-.','color',rgbcolor(i,:));
        end
        
        % interpolating values in airfoil thickness direction according to
        % discretization in chord-wise direction to calculate whole
        % cross-section surface coordinates
        xhi_airf_u = fun_lin_interpolation(airfoil_coord_u,xhi_y_w);
        xhi_airf_l = fun_lin_interpolation(airfoil_coord_l,xhi_y_w);
        
        % coordinates of whole cross-section surface
        xhi_airf_w = [xhi_airf_u(1:end-1);(xhi_airf_u(end)+xhi_airf_l(end))/2;xhi_airf_l(end-1:-1:1)];

        % if-statement to detect, if blade root ends. In case blade root
        % ends, switching from whole surface to camber surface
        if strncmp(blade_vl_model,'2D',2)
            if i == nbr+1
                if nbr ~= 0
                    xhi_airf_w = [xhi_airf_c(1:end-1);xhi_airf_c(end:-1:1)];
                end
            end
        end
        arr_val_c(:,i)   = xhi_airf_c;      % camber surface coordinates 
        arr_val_w(:,i)   = xhi_airf_w(:,1); % whole surface coordinates
    end
    
    legend(hpp,strlegend,'Interpreter', 'none');
    
    % interpolating camber surface coordinates in span-wise direction 
    for i = 1:size(arr_val_c,1)
        arr_val = interp1(uvlm_obj.airfoil.grid,arr_val_c(i,:)',xhi_x,'linear');
        arr_xhi_airf_c(:,i) = arr_val; %(:,1);
    end
    
    % interpolating whole surface coordinates in span-wise direction
    for i = 1:size(arr_val_w,1)
        arr_val = interp1(uvlm_obj.airfoil.grid,arr_val_w(i,:),xhi_x,'linear');
        arr_xhi_airf_w(:,i) = arr_val; %(:,1);
    end
    
%     rgbcolor = RGB_Color;
%     figure(); hold on; grid on;
%     mm = length(xhi_y_w)
%     xhi_y_w_mm = [xhi_y_w(1:end),xhi_y_w(end-1:-1:1)];
%     for i = 1:size(arr_xhi_airf_w,1)
%        plot(xhi_y_w_mm,arr_xhi_airf_w(i,:),'-','color',rgbcolor(i,:));
% %        plot(xhi_y_w_mm,arr_xhi_airf_w(i,:),'-','color',rgbcolor(i,:));
%     end

    if strncmp(blade_vl_model,'2D',2)
        if nbr~=0
            inz  = find(arr_inz_airfoil<=nbr);
            if length(inz)<length(arr_inz_airfoil)
                temp_coo = arr_xhi_airf_c(inz(end)+1,:);
                arr_xhi_airf_w(inz(end)+1,:) = [temp_coo(1:end),temp_coo(end-1:-1:1)];
            end
        end
    end
    
    
%     figure(); hold on; grid on;
%     xhi_y_w_mm = [xhi_y_w(1:end),xhi_y_w(end-1:-1:1)];
%     for i = 1:size(arr_xhi_airf_w,1)
%         plot(xhi_y_w_mm,arr_xhi_airf_w(i,:));
%     end
    
    uvlm_obj.arr_xhi_airf_c  = arr_xhi_airf_c;
    uvlm_obj.arr_xhi_airf_w  = arr_xhi_airf_w;
    uvlm_obj.xhi_y_w         = xhi_y_w;
    uvlm_obj.xhi_y_c         = xhi_y_c;
    uvlm_obj.arr_xhi_x       = xhi_x;
    uvlm_obj.nbr             = nbr;
    uvlm_obj.arr_inz_airfoil = arr_inz_airfoil;
    
return
