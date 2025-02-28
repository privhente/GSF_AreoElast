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
function beam = fun_extract_beam_data(strName,model_beam,materials,mesh_type_opt,mesh_factor_opt)
% =================================================================================================================
    mesh_type   = 0;
    mesh_factor = 0;
    
    if nargin >= 4; mesh_type = mesh_type_opt; end
    if nargin >= 5; mesh_factor = mesh_factor_opt; end

    % span-wise discretization
    beam.M = model_beam.M_struc;

    % initializing stiffness and mass matrix
    beam.arr_stiff_matrix = zeros(beam.M,21);
    beam.arr_mass_matrix = zeros(beam.M,6);
    O = 0*ones(beam.M,1);

    % natural coordinates of span-and chord-wise discretization
    xhi_x_li  = [0:1/beam.M:1];
    delta_xhi_c = 1/beam.M;

    % cos distributed mesh size in span-wise direction:
    mc = delta_xhi_c;
    mf = mc/(1/(1-mesh_factor+1e-6));
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
    
    xhi_x      = xhi_x_nl;
    beam.xhi_x = xhi_x;
    
    % element natural coordinates in span-wise direction
    xhi_x_elem = (xhi_x(1:end-1)+xhi_x(2:end))/2;

    % interpolating coordinates of reference axis according to span-wise discretization
    beam.arr_xre_x = interp1(model_beam.reference_axis.x.grid,model_beam.reference_axis.x.values,xhi_x,'linear')';
    beam.arr_xre_y = interp1(model_beam.reference_axis.y.grid,model_beam.reference_axis.y.values,xhi_x,'linear')';
    beam.arr_xre_z = interp1(model_beam.reference_axis.z.grid,model_beam.reference_axis.z.values,xhi_x,'linear')';
%     figure(); title('ref_x beam'); hold on; grid on; plot(xhi_x,beam.arr_xre_x,'xr'); plot(model_beam.reference_axis.x.grid,model_beam.reference_axis.x.values,'-k');
%     figure(); title('ref_y beam'); hold on; grid on; plot(xhi_x,beam.arr_xre_y,'xr'); plot(model_beam.reference_axis.y.grid,model_beam.reference_axis.y.values,'-k');
%     figure(); title('ref_z beam'); hold on; grid on; plot(xhi_x,beam.arr_xre_z,'xr'); plot(model_beam.reference_axis.z.grid,model_beam.reference_axis.z.values,'-k');

    % airfoil orientation from AeroDyn-File
    if isfield(model_beam,'airfoil_orientation')
        beam.airfoil_orientation = interp1(model_beam.airfoil_orientation.grid,model_beam.airfoil_orientation.values,xhi_x,'linear')';
%         figure(); title('airfoil_orientation'); hold on; grid on; plot(xhi_x,beam.airfoil_orientation,'xr'); plot(model_beam.airfoil_orientation.grid,model_beam.airfoil_orientation.values,'-k');
    end

    % interpolating twist angle according to span-wise discretization
    if isfield(model_beam,'twist');
        beam.arr_twist = interp1(model_beam.twist.grid,model_beam.twist.values,xhi_x,'linear')';
%         figure(); title('twist beam'); hold on; grid on; plot(xhi_x,beam.arr_twist,'xr'); plot(model_beam.twist.grid,model_beam.twist.values,'-k');
    end
    
    % if-statement according to type of surface cross-section
    if strfind(strName,'blade') ~= 0
        % interpolating along span-wise direction
        if isfield(model_beam.elastic_properties_mb,'six_x_six')
            beam.arr_stiff_matrix = [];
            beam.arr_mass_matrix  = [];
            % stiffness and mass/interia terms should be already converted to DeSiO-Format, i.e. Voigt notation
%             figure(); title('stiff beam'); hold on; grid on;
            for i = 1:21
                beam.arr_stiff_matrix(:,i) = interp1(model_beam.elastic_properties_mb.six_x_six.stiff_matrix.grid,model_beam.elastic_properties_mb.six_x_six.stiff_matrix.values(:,i),xhi_x_elem,'linear');
%                 plot(xhi_x_elem,beam.arr_stiff_matrix(:,i),'xr'); plot(model_beam.elastic_properties_mb.six_x_six.stiff_matrix.grid,model_beam.elastic_properties_mb.six_x_six.stiff_matrix.values(:,i),'-k');
            end
%             figure(); title('mass beam'); hold on; grid on;
            for i = 1:6
                beam.arr_mass_matrix(:,i)  = interp1(model_beam.elastic_properties_mb.six_x_six.inertia_matrix.grid,model_beam.elastic_properties_mb.six_x_six.inertia_matrix.values(:,i),xhi_x_elem,'linear');
%                 plot(xhi_x_elem,beam.arr_mass_matrix(:,i),'xr'); plot(model_beam.elastic_properties_mb.six_x_six.inertia_matrix.grid,model_beam.elastic_properties_mb.six_x_six.inertia_matrix.values(:,i),'-k');
            end
            beam.dissipation(1:2) = [model_beam.dissipation.alpha_s, model_beam.dissipation.alpha_v];
        else
            
            beam.arr_EA  = zeros(beam.M,1); beam.arr_GA1 = zeros(beam.M,1); beam.arr_GA2 = zeros(beam.M,1);
            beam.arr_EI1 = zeros(beam.M,1); beam.arr_EI2 = zeros(beam.M,1); beam.arr_GI3 = zeros(beam.M,1);
            beam.arr_ES1 = zeros(beam.M,1); beam.arr_ES2 = zeros(beam.M,1); beam.arr_EI12 = zeros(beam.M,1);
            beam.arr_GS1 = zeros(beam.M,1); beam.arr_GS2 = zeros(beam.M,1);
            
            beam.arr_rhoA  = zeros(beam.M,1); beam.arr_rhoI1 = zeros(beam.M,1); beam.arr_rhoI2  = zeros(beam.M,1);
            beam.arr_rhoS1 = zeros(beam.M,1); beam.arr_rhoS2 = zeros(beam.M,1); beam.arr_rhoI12 = zeros(beam.M,1);
            
            if isfield(model_beam.elastic_properties_mb,'EA')
                beam.arr_EA   = interp1(model_beam.elastic_properties_mb.EA.grid,model_beam.elastic_properties_mb.EA.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'GA1')
                beam.arr_GA1  = interp1(model_beam.elastic_properties_mb.GA1.grid,model_beam.elastic_properties_mb.GA1.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'GA2')
                beam.arr_GA2  = interp1(model_beam.elastic_properties_mb.GA2.grid,model_beam.elastic_properties_mb.GA2.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'EI1')
                beam.arr_EI1  = interp1(model_beam.elastic_properties_mb.EI1.grid,model_beam.elastic_properties_mb.EI1.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'EI2')
                beam.arr_EI2  = interp1(model_beam.elastic_properties_mb.EI2.grid,model_beam.elastic_properties_mb.EI2.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'GI3')
                beam.arr_GI3  = interp1(model_beam.elastic_properties_mb.GI3.grid,model_beam.elastic_properties_mb.GI3.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'ES1')
                beam.arr_ES1  = interp1(model_beam.elastic_properties_mb.ES1.grid,model_beam.elastic_properties_mb.ES1.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'ES2')
                beam.arr_ES2  = interp1(model_beam.elastic_properties_mb.ES2.grid,model_beam.elastic_properties_mb.ES2.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'GS1')
                beam.arr_GS1  = interp1(model_beam.elastic_properties_mb.GS1.grid,model_beam.elastic_properties_mb.GS1.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'GS2')
                beam.arr_GS2  = interp1(model_beam.elastic_properties_mb.GS2.grid,model_beam.elastic_properties_mb.GS2.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'EI12')
                beam.arr_EI12 = interp1(model_beam.elastic_properties_mb.EI12.grid,model_beam.elastic_properties_mb.EI12.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoA')
                beam.arr_rhoA   = interp1(model_beam.elastic_properties_mb.rhoA.grid,model_beam.elastic_properties_mb.rhoA.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoI1')
                beam.arr_rhoI1  = interp1(model_beam.elastic_properties_mb.rhoI1.grid,model_beam.elastic_properties_mb.rhoI1.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoI2')
            beam.arr_rhoI2  = interp1(model_beam.elastic_properties_mb.rhoI2.grid,model_beam.elastic_properties_mb.rhoI2.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoS1')
                beam.arr_rhoS1  = interp1(model_beam.elastic_properties_mb.rhoS1.grid,model_beam.elastic_properties_mb.rhoS1.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoS2')
                beam.arr_rhoS2  = interp1(model_beam.elastic_properties_mb.rhoS2.grid,model_beam.elastic_properties_mb.rhoS2.values,xhi_x_elem,'linear');
            end
            
            if isfield(model_beam.elastic_properties_mb,'rhoI12')
                beam.arr_rhoI12 = interp1(model_beam.elastic_properties_mb.rhoI12.grid,model_beam.elastic_properties_mb.rhoI12.values,xhi_x_elem,'linear');
            end
            % new test:
            zero_array = interp1([0.0, 1.0], [0.0, 0.0], xhi_x_elem, 'linear');
            % GA1 GA2 EA EI1 EI2 GI3 0 0 0 GS2 -GS1 0 0 0 0 0 ES1 -EI12 -ES2 0 0
%             beam.arr_stiff_matrix(:,1:21) = [beam.arr_GA1, beam.arr_GA2, beam.arr_EA, beam.arr_EI1, beam.arr_EI2, beam.arr_GI3, ...
%                                             O, O, O, beam.arr_GS2, -beam.arr_GS1, O, O, O, O, O, beam.arr_ES1, -beam.arr_EI12, -beam.arr_ES2, O, O];
            beam.arr_stiff_matrix(:,1:21) = [beam.arr_GA1(:), beam.arr_GA2(:), beam.arr_EA(:), beam.arr_EI1(:), beam.arr_EI2(:), beam.arr_GI3(:), ...
                                            zero_array(:), zero_array(:), zero_array(:), beam.arr_GS2(:), -beam.arr_GS1(:), ...
                                            zero_array(:), zero_array(:), zero_array(:), zero_array(:), zero_array(:), beam.arr_ES1(:),...
                                            -beam.arr_EI12(:), -beam.arr_ES2(:), zero_array(:), zero_array(:)];  
            % rhoA, rhoI2, rhoI1, rhoI12, rhoS1, rhoS2
            beam.arr_mass_matrix(:,1:6)   = [beam.arr_rhoA(:), beam.arr_rhoI2(:), beam.arr_rhoI1(:), beam.arr_rhoI12(:), beam.arr_rhoS1(:), beam.arr_rhoS2(:)];
            beam.dissipation(1:2)         = [model_beam.dissipation.alpha_s, model_beam.dissipation.alpha_v];
        end
    end
    if strfind(strName,'pipe')~=0
        % interpolating cross-section properties along span-wise direction
        E = 0.0; G = 0.0; rho = 0.0;
        if isfield(model_beam,'material')
            strmaterial = model_beam.material;
            
            % Abfrage, ob cell oder structure variable ! 
            arr_struct = whos('materials');
            if strncmp(arr_struct.class,'struct',6)
                for i = 1:size(materials,2)
                    tab_materials{i} = struct2table(materials(i));
                end
            elseif strncmp(arr_struct.class,'cell',4)
                for i = 1:size(materials,1)
                    tab_materials{i} = cell2table(materials(i));
                end
            end
            materials = tab_materials;
            for i = 1:size(materials,2)
%15MW original                 
                if strncmp(materials{1, i}.name,strmaterial,length(strmaterial))
                    if length(materials{i}.name) == length(strmaterial)
                        E   = materials{i}.E;
                        nu  = materials{i}.nu;
                        rho = materials{i}.rho;
                        G   = E/(2.0*(1.0+nu));
                    end
                    break
                end
%  5/25MW modifie               
%                   if strncmp(materials{1, i}.Var1.name,strmaterial,length(strmaterial))
%                     if length(materials{1, i}.Var1.name) == length(strmaterial)
%                         E   = materials{1, i}.Var1.E;
%                         nu  = materials{1, i}.Var1.nu;
%                         rho = materials{1, i}.Var1.rho;
%                         G   = E/(2.0*(1.0+nu));
%                         disp(materials{1, i}.Var1.name);
%                     end
%                     break
%                 end
            end
        end
        if isfield(model_beam,'elastic_properties_mb')
            E   = model_beam.elastic_properties_mb.E;
            G   = E/(2.0*(1+model_beam.elastic_properties_mb.nu));
            rho = model_beam.elastic_properties_mb.rho;
        end

        k1  = 1.0;
        k2  = 1.0;
        if isfield(model_beam,'shear_factor')
            k1 = model_beam.shear_factor.k1;
            k2 = model_beam.shear_factor.k2;
        end

        beam.arr_outer_diameter = interp1(model_beam.outer_diameter.grid,model_beam.outer_diameter.values,xhi_x_elem,'linear')';
        beam.arr_thickness      = interp1(model_beam.thickness.grid,model_beam.thickness.values,xhi_x_elem,'linear')';
        beam.arr_EA             = E*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_GA1            = k1*G*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_GA2            = k2*G*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_EI1            = E*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_EI2            = E*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_GI3            = G*pi/32*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_ES1            = beam.arr_EA*0;
        beam.arr_ES2            = beam.arr_EA*0;
        beam.arr_GS2            = beam.arr_EA*0;
        beam.arr_GS1            = beam.arr_EA*0;
        beam.arr_EI12           = beam.arr_EA*0;

        beam.arr_rhoA           = rho*pi/4*(beam.arr_outer_diameter.^2 - (beam.arr_outer_diameter-2*beam.arr_thickness).^2);
        beam.arr_rhoI1          = rho*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_rhoI2          = rho*pi/64*(beam.arr_outer_diameter.^4 - (beam.arr_outer_diameter-2*beam.arr_thickness).^4);
        beam.arr_rhoI12         = beam.arr_rhoA*0;
        beam.arr_rhoS1          = beam.arr_rhoA*0;
        beam.arr_rhoS2          = beam.arr_rhoA*0;

        % GA1 GA2 EA EI1 EI2 GI3 0 0 0 GS2 -GS1 0 0 0 0 0 ES1 -EI12 -ES2 0 0
        beam.arr_stiff_matrix(:,1:21) = [beam.arr_GA1, beam.arr_GA2, beam.arr_EA, beam.arr_EI1, beam.arr_EI2, beam.arr_GI3, ...
                                        O, O, O, beam.arr_GS2, -beam.arr_GS1, O, O, O, O, O, beam.arr_ES1, -beam.arr_EI12, -beam.arr_ES2, O, O];

        beam.arr_mass_matrix(:,1:6) = [beam.arr_rhoA, beam.arr_rhoI2, beam.arr_rhoI1, beam.arr_rhoI12, beam.arr_rhoS1, beam.arr_rhoS2];
        beam.dissipation(1:2)       = [ model_beam.dissipation.alpha_s, model_beam.dissipation.alpha_v];
    end
    
return
