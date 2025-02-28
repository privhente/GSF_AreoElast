% =================================================================================================================
addpath('H:\10_DeSiO\DeSiO-Aero\pre_post_processing\versions\v02_062021');

close all;
clear all;
clc;

currDir = cd;
% =================================================================================================================
straddfilename = 'wakelines1';
% NACA-mpxx
NACA = [0 0 0 9];
% Aspect ratio
AR = 4; 
% chord, profil length, here in y direction
c  = 1;
% meshing parameter
mx = 13; my = 10;
% free field velocity, axis of wind direction
vinf = 10;
% angle of attack of vinf, here (90-alpha), axis of wind direction at alpha = 90°
alpha = 5; n_axis= [1,0,0]'; % reference axis orthogonal to wind 
% wind density
density = 1.0;
% cutoff
cutoff = 0.01;
% blade length, here in x direction
L  = AR*c;
% Total simulation time
time = 1.0;
% max. wake cutting distance
factor_wcd = [1];
% time at wake is cut from simulation
time_cut = time;
% Simulation step time
deltat = 0.0625; % sqrt((c*L)/(my*mx))/vinf;

for i_arr = 1:length(factor_wcd)
    
    wcd = factor_wcd(i_arr)*time*vinf;
    % Number of waking rows to consider - After that number, wake is cut
    tnrows = round(time_cut/deltat);
    nrows  = round(wcd/(vinf*deltat));
    
    % ============================================================================================================================================
    % WRITING DESIO-AERO INPUT FILES
    % ============================================================================================================================================
    % max. curvature fraction
    m = NACA(1)/100;
    % distance of max. curvature fraction from profil nose
    p = NACA(2)*10/100;
    % airfoil thickness fraction xx
    t = (NACA(3)*10 + NACA(4))/100;

    % DeSiO-Input parameter
    strfilename = ['NACA_' num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)) ['_' straddfilename '_'] 'L' num2str(L) '_c' num2str(c) '_mx' num2str(mx) '_my' num2str(my) '_alpha' num2str(alpha) '_vinf' num2str(vinf) '_nrows' num2str(nrows)];

    % creating wake lines: [surface, property, vecA,vecB]
    wake_line = [1, 1, L,  c,  0, 0,  c,  0];

    %         wake_line = [1, 1, 0,  c,  0, 0,  0,  0; ...
    %                      1, 1, L,  0,  0, L,  c,  0; ...
    %                      1, 1, L,  c,  0, 0,  c,  0];

    wake_prop = [tnrows, nrows];

    % Creating velocity vector
    R      = cos(alpha*pi/180)*eye(3) + sin(alpha*pi/180)*skew(n_axis)+(1-cos(alpha*pi/180))*n_axis*n_axis';
    d_vinf = R*[0;1;0];

    % creating node coordinates
    k = 1; a = 0; xcoord0 = 0; nodes = []; ycoord = 0; profil = [];
    nn = (mx+1)*(my+1);
    for i = 1:nn
        if i == (mx+1)*k+1 
            xcoord0 = 0;
            ycoord  = k*c/my;
            k = k + 1; a = 0;
        end
        [airfoil]     = fun_NACA_airfoil(m,p,c,t,ycoord);
        xcoord        = xcoord0 + L/mx*a;
        zcoord        = airfoil.yc;
        nodes(i,1:4)  = [i,xcoord,ycoord,zcoord];
        a = a + 1;

        profil(i).upperline = [xcoord,airfoil.xu,airfoil.yu];
        profil(i).lowerline = [xcoord,airfoil.xl,airfoil.yl];
        profil(i).camber    = [xcoord,ycoord,zcoord];
    end

    % creating surface connectivity
    surfaces = [1:mx*my]; connectivity = [];
    for j = 1:my
        n1 = surfaces((j-1)*mx+1:(j-1)*mx+mx) + mx + j;
        n2 = n1+1;
        n3 = n1-mx;
        n4 = n1-(mx+1);
        connectivity = [connectivity; [n2',n1',n4',n3']];
    end

    wakes = []; wake.n_wakes = size(wake_line,1);
    for i = 1: wake.n_wakes
        vec_e = wake_line(i,6:8)-wake_line(i,3:5);
        [inz_x] = find(nodes(:,2)     >= min(wake_line(i,3),wake_line(i,6)) & nodes(:,2)     <= max(wake_line(i,3),wake_line(i,6)) );
        [inz]   = find(nodes(inz_x,3) >= min(wake_line(i,4),wake_line(i,7)) & nodes(inz_x,3) <= max(wake_line(i,4),wake_line(i,7)) );
        % ascend ordering according to angle
        [temp,inz_alpha] = sort(vec_e*nodes(inz_x(inz),2:4)');
        inz = inz_x(inz(inz_alpha));
        wakes(i).nsurf = wake_line(i,1);
        wakes(i).nseg  = length(inz);
        wakes(i).prop  = wake_line(i,2);
        for j = 1:length(inz)-1
            [row1,col1,v] = find(connectivity == [inz(j)]);
            [row2,col2,v] = find(connectivity(row1,:) == [inz(j+1)]);
            wakes(i).inf(j,1:3) = [inz(j:j+1)', row1(row2)];
        end
    end

    wake_props = []; wake.nprop = size(wake_prop,1);
    for i = 1:size(wake_prop,1)
        wake_props(i).tnrows = wake_prop(i,1);
        wake_props(i).nrows  = wake_prop(i,2);
    end

    % assigning simulation, mesh and wake data to matlab structure variables
    simu.time = time; simu.deltat = deltat; simu.cutoff = cutoff;
    simu.density = density; simu.i_vinf  = vinf; simu.d_vinf  = d_vinf;
    simu.currDir = currDir; simu.strfilename = strfilename;
    mesh.nodes = nodes; mesh.connectivity = connectivity; 
    mesh.mx = mx; mesh.my = my; mesh.nn = nn;
    wake.wakes = wakes; wake.wake_prop = wake_props;

    % plotting of the system
%     figure(); hold on; grid on; axis equal; view(45,30); axis([0 L 0 c -t*c t*c]);
%     xlabel('x-direction length'); 
%     ylabel('y-direction chord');
%     zlabel('z-direction thickness'); 
%     for i = 1:size(connectivity,1)
%         n1 = connectivity(i,1);  n2 = connectivity(i,2);  n3 = connectivity(i,3);  n4 = connectivity(i,4);
%         camber_coord    = [profil(n1).camber   ; profil(n2).camber   ; profil(n3).camber   ; profil(n4).camber];
%         lowerline_coord = [profil(n1).lowerline; profil(n2).lowerline; profil(n3).lowerline; profil(n4).lowerline];
%         upperline_coord = [profil(n1).upperline; profil(n2).upperline; profil(n3).upperline; profil(n4).upperline];
% 
%         surf = fill3(camber_coord(:,1)   , camber_coord(:,2)   , camber_coord(:,3)   ,'k','facealpha',0.1,'edgealpha',0.0);
%         surf = fill3(lowerline_coord(:,1), lowerline_coord(:,2), lowerline_coord(:,3),'g','facealpha',0.1,'edgealpha',0.5);
%         surf = fill3(upperline_coord(:,1), upperline_coord(:,2), upperline_coord(:,3),'g','facealpha',0.1,'edgealpha',0.5);
% 
%         s_factor = 0.6;
%         coordm = (camber_coord(1,:) + camber_coord(3,:))/2;
%         camber_coord_quiv = s_factor*camber_coord + [coordm;coordm;coordm;coordm]*(1-s_factor);
%         a = 0; b = 0;
%         for j = 1:4
%             a = a + 1;
%             if a == 4
%                 b = 1;
%             else
%                 b = a + 1;
%             end
%             d_quiv = (camber_coord_quiv(b,:)-camber_coord_quiv(a,:))/norm((camber_coord_quiv(b,:)-camber_coord_quiv(a,:)))*s_factor;
% %             quiver3(camber_coord_quiv(a,1),camber_coord_quiv(a,2),camber_coord_quiv(a,3),d_quiv(1),d_quiv(2),d_quiv(3),'color','r','linewidth',1);
%             if a == 1
% %                 plot3([camber_coord_quiv(a,1)],[camber_coord_quiv(a,2)],[camber_coord_quiv(a,3)],'.k','markersize',10);
%             end
%         end
% %         text(coordm(1),coordm(2),coordm(3),['G' num2str(i)],'fontsize',6);
%     end

%     for i = 1:wake.n_wakes
%         wake_line(i,6:8)-wake_line(i,3:5);
%         plot3([wake_line(i,3),wake_line(i,6)] , [wake_line(i,4),wake_line(i,7)], [wake_line(i,5),wake_line(i,8)] , '-b','linewidth',2);
%     end
% 
%     for i = 1:size(nodes,1)
% %         text(nodes(i,2), nodes(i,3), nodes(i,4), num2str(nodes(i,1)))
%     end

    % writing DeSiO-Aero files
    bol = fun_writeDeSiOAeroinput(simu,mesh,wake);
    close all;

    simu = [];
    mesh = [];
    wake = [];
end
return