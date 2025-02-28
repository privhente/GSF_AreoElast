function fun_DeSiO_create_input_Wing(mx,my,vinf)

my   = [4];
mx   = [20];
vinf = [10];
flag_support = 1;
flag_low_inertia = 0;
%   == 0 - all free
%   == 1 - rigid on left side
%   == 2 - rigid+spherical
%   == 3 - rigid+rigid
%   == 4 - sphericalsupport+rotationlocal_z+rotation_y
%   == 5 - sphericalsupport+sphericalsupport+rotation_z

% =================================================================================================================
addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));
currDir = cd;
% =================================================================================================================
% NACA-mpxx
NACA = [0 0 1 2];
% Aspect ratio
AR = 10/1; 
% chord, profil length, here in y direction
c  = 1.00;
% inclination of plate
phi = 0.0*pi/180;
% meshing parameter
% angle of attack of vinf, here (90-alpha), axis of wind direction at alpha = 90°
alpha = 15.0; n_axis= [1,0,0]'; % reference rotation axis orthogonal to wind 
% wind density
density = 1.225;
% cutoff
cutoff = 0.01;
% blade length, here in x direction
L  = AR*c;
% time at wake is cut from simulation
% Simulation step time
deltaL = sqrt(c*L/(mx*my));
% deltat = deltaL/vinf;
deltat = (c/my)/vinf;
time   = 3;
time_cut = time;
nrows = 200; %round(time_cut/deltat);
tnrows = 200; %nrows;
% Number of waking rows to consider - After that number, wake is cut
% wake diffusion
wake_diff = 0.0;

if flag_low_inertia==1
    strfilename = ['ex_plate_' num2str(flag_support) '_red_inertia'];
else
    strfilename = ['ex_plate_' num2str(flag_support)];
end
% structure
deltat_struc = deltat;
ne_struc     = mx;
cmat = [
    5.60000000E+08	2.10526316E+08	2.10526316E+08	2.98666667E+03	4.66666667E+07	1.754498246E+07	0.00000000E+00	0.00000000E+00	0.00000000E+00	0.00000000E+00	0.00000000E+00
];
mmat = [
    2.16000000E+01	1.15200000E-04	1.80000000E+00	0.00000000E+00	0.00000000E+00	0.00000000E+00
];
if flag_low_inertia==1
    mmat = [2.16000000E-02	1.15200000E-04	1.80000000E+00	0.00000000E+00	0.00000000E+00	0.00000000E+00];
end
diss = [0.0, 0.0];
simType = 'dynamic';
simfsitype = 'strong';
lfType  = 'constant';
lf_duration = 100*deltat;
flag_linearization = 0;
flag_output_matrix = 0;

% =========================================================================
% WRITING DESIO-AERO INPUT FILES
% =========================================================================
% max. curvature fraction
m = NACA(1)/100;
% distance of max. curvature fraction from profil nose
p = NACA(2)*10/100;
% airfoil thickness fraction xx
t = (NACA(3)*10 + NACA(4))/100;

% creating wake lines: [surface, property, vecA,vecB]
wake_line = [1, 1, L,  c,  0, 0,  c,  0];

% wake_line = [1, 1, 0,  c,  0, 0,  0,  0; ...
%              1, 2, L,  0,  0, L,  c,  0; ...
%              1, 3, L,  c,  0, 0,  c,  0];

wake_prop = [tnrows, nrows, wake_diff;...
    tnrows, nrows, wake_diff;...
    tnrows, nrows, wake_diff];

% Creating velocity vector
R      = cos(alpha*pi/180)*eye(3) + sin(alpha*pi/180)*skew(n_axis)+(1-cos(alpha*pi/180))*n_axis*n_axis';
d_vinf = (R*[0;1;0])';

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
    nodes(i,1:4)  = [i,xcoord,ycoord*cos(phi),ycoord*sin(phi)];
    nodes0(i,1:4) = [i,xcoord,ycoord,0];
    a = a + 1;

    profil(i).upperline = [xcoord,airfoil.xu*cos(phi),airfoil.yu*sin(phi)];
    profil(i).lowerline = [xcoord,airfoil.xl*cos(phi),airfoil.yl*sin(phi)];
    profil(i).camber    = [xcoord,ycoord*cos(phi),ycoord*sin(phi)];
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
    [inz_x] = find(nodes0(:,2)     >= min(wake_line(i,3),wake_line(i,6)) & nodes0(:,2)     <= max(wake_line(i,3),wake_line(i,6)) );
     inz_x  = [1:size(nodes0,1)]';
    
    [inz]   = find(nodes0(inz_x,3) >= min(wake_line(i,4),wake_line(i,7)) & nodes0(inz_x,3) <= max(wake_line(i,4),wake_line(i,7)) );
    % ascend ordering according to angle
    [temp,inz_alpha] = sort(vec_e*nodes0(inz_x(inz),2:4)');
    inz = inz_x(inz(inz_alpha));
    wakes(i).nsurf      = wake_line(i,1);
    wakes(i).nsegments  = length(inz)-1;
    wakes(i).nproperty  = wake_line(i,2);
    for j = 1:length(inz)-1
        [row1,col1,v] = find(connectivity == [inz(j)]);
        [row2,col2,v] = find(connectivity(row1,:) == [inz(j+1)]);
        wakes(i).inf(j,1:3) = [inz(j:j+1)', row1(row2)];
    end
end

% plotting of the system
figure(); hold on; grid on; axis equal; view(45,30); axis([0 L 0 c -t*c t*c]);
xlabel('x-direction length'); 
ylabel('y-direction chord');
zlabel('z-direction thickness'); 
for i = 1:size(connectivity,1)
    n1 = connectivity(i,1);  n2 = connectivity(i,2);  n3 = connectivity(i,3);  n4 = connectivity(i,4);
    camber_coord    = [profil(n1).camber   ; profil(n2).camber   ; profil(n3).camber   ; profil(n4).camber];
    lowerline_coord = [profil(n1).lowerline; profil(n2).lowerline; profil(n3).lowerline; profil(n4).lowerline];
    upperline_coord = [profil(n1).upperline; profil(n2).upperline; profil(n3).upperline; profil(n4).upperline];

    surf = fill3(camber_coord(:,1)   , camber_coord(:,2)   , camber_coord(:,3)   ,'k','facealpha',0.1,'edgealpha',0.5);
    surf = fill3(lowerline_coord(:,1), lowerline_coord(:,2), lowerline_coord(:,3),'g','facealpha',0.1,'edgealpha',0.1);
    surf = fill3(upperline_coord(:,1), upperline_coord(:,2), upperline_coord(:,3),'g','facealpha',0.1,'edgealpha',0.1);

    s_factor = 0.6;
    coordm = (camber_coord(1,:) + camber_coord(3,:))/2;
    camber_coord_quiv = s_factor*camber_coord + [coordm;coordm;coordm;coordm]*(1-s_factor);
    a = 0; b = 0;
    for j = 1:4
        a = a + 1;
        if a == 4
            b = 1;
        else
            b = a + 1;
        end
        d_quiv = (camber_coord_quiv(b,:)-camber_coord_quiv(a,:))/norm((camber_coord_quiv(b,:)-camber_coord_quiv(a,:)))*s_factor;
        quiver3(camber_coord_quiv(a,1),camber_coord_quiv(a,2),camber_coord_quiv(a,3),d_quiv(1),d_quiv(2),d_quiv(3),'color','r','linewidth',1);
        if a == 1
            plot3([camber_coord_quiv(a,1)],[camber_coord_quiv(a,2)],[camber_coord_quiv(a,3)],'.k','markersize',10);
        end
    end
    text(coordm(1),coordm(2),coordm(3),['G' num2str(i)],'fontsize',6,'color','r','fontweight','bold');
end

for i = 1:wake.n_wakes
    wake_line(i,6:8)-wake_line(i,3:5);
    plot3([wake_line(i,3),wake_line(i,6)] , [wake_line(i,4),wake_line(i,7)], [wake_line(i,5),wake_line(i,8)] , '-b','linewidth',2);
end

for i = 1:size(nodes,1)
    text(nodes(i,2), nodes(i,3), nodes(i,4), num2str(nodes(i,1)))
end

% writing DeSiO-Aero files
simu_aero.strfilename = [strfilename];
simu_aero.sort = 'constant'; 
simu_aero.time = time; 
simu_aero.deltat = deltat;
simu_aero.cutoff = cutoff;
simu_aero.density = density;
simu_aero.i_vinf  = vinf; 
simu_aero.d_vinf = d_vinf;
simu_aero.currDir = currDir;
simu_aero.strfilename = strfilename;
mesh.nodes = nodes(:,2:4);
mesh.connectivity = connectivity; 
mesh.mx = mx; 
mesh.my = my;
mesh.nn = nn;
wake.wakes = wakes;
wake.property = wake_prop;
wake.nprop = size(wake_prop,1);
fun_writeDeSiOAeroInput(simu_aero,mesh,wake);

% setting structural model: beam elements
nnode_struc = ne_struc + 1;
mesh.strname = 'wing'
for i = 1:nnode_struc
    d1 = [0,1,0]'; d2 = [0,0,1]'; d3 = [1,0,0]';
    R = [1,0,0;0,cos(phi) -sin(phi); 0,sin(phi) cos(phi)];
    coord = [L/(nnode_struc-1)*(i-1),0.5*c,0]';
    mesh_struc.nodes(i,1:12) = [(R*coord)', (R*d1)', (R*d2)', (R*d3)'];
end
for i = 1:ne_struc
   mesh_struc.connectivity(i,:) = [i,i+1,1];
end
mesh_struc.nn = nnode_struc;
mesh_struc.mx = ne_struc;

simu_struct.flag_matbeam = 0;
simu_struct.matbeam.cmat = cmat;
simu_struct.matbeam.mmat = mmat; 
simu_struct.matbeam.diss = diss; 

% constraints
simu_struct.constraints = [];
remo_internal_const = [];
if flag_support == 0
    % nothing
elseif flag_support == 1
    simu_struct.constraints(1).type  = 'rigidsupport';
    simu_struct.constraints(1).nodes = [1 0];
    remo_internal_const = [1];
elseif flag_support == 2
    simu_struct.constraints(1).type  = 'rigidsupport';
    simu_struct.constraints(1).nodes = [1 0];
    simu_struct.constraints(2).type  = 'sphericalsupport';
    simu_struct.constraints(2).nodes = [nnode_struc 0];
    remo_internal_const = [1];
elseif flag_support == 3
    simu_struct.constraints(1).type  = 'rigidsupport';
    simu_struct.constraints(1).nodes = [1 0];
    simu_struct.constraints(2).type  = 'rigidsupport';
    simu_struct.constraints(2).nodes = [nnode_struc 0];
    remo_internal_const = [1,nnode_struc];
elseif flag_support == 4
    simu_struct.constraints(1).type  = 'sphericalsupport';
    simu_struct.constraints(1).nodes = [1 0];
    simu_struct.constraints(2).type  = 'rotation_local';
    simu_struct.constraints(2).nodes = [1 0];
    simu_struct.constraints(2).dir   = [0,0,1];
    simu_struct.constraints(3).type  = 'rotation_local';
    simu_struct.constraints(3).nodes = [1 0];
    simu_struct.constraints(3).dir   = [0,1,0];
elseif flag_support == 5
    simu_struct.constraints(1).type  = 'sphericalsupport';
    simu_struct.constraints(1).nodes = [1 0];
    simu_struct.constraints(2).type  = 'sphericalsupport';
    simu_struct.constraints(2).nodes = [nnode_struc 0];
    simu_struct.constraints(3).type  = 'rotation_local';
    simu_struct.constraints(3).nodes = [1 0];
    simu_struct.constraints(3).dir   = [0,0,1];
end 

nodes = 1:nnode_struc;
ai = size(simu_struct.constraints,2);
for i = 1:nnode_struc
    if all(remo_internal_const-nodes(i))
        ai = ai + 1;
        simu_struct.constraints(ai).type  = 'internal';
        simu_struct.constraints(ai).nodes = [nodes(i) 0];
    end
end

% simulation settings
simu_struct.currDir     = currDir;
simu_struct.strfilename = strfilename;
simu_struct.type        = 'modal'
simu_struct.time        = time;
simu_struct.deltat      = deltat_struc;
simu_struct.tol         = 1.0e-8;
simu_struct.niter       = 50;
simu_struct.grav        = 0;
simu_struct.grav_vec    = [0,0,-9.81];

% write structure files
fun_writeDeSiOStructureInput(simu_struct,mesh_struc);

% fsi settings
simu_fsi.lf_type = lfType;
simu_fsi.lf_duration = lf_duration;
simu_fsi.currDir     = currDir;
simu_fsi.strfilename = [strfilename];
simu_fsi.radius_rbf  = c/2+10^(floor(log10(c/2)))^2;
simu_fsi.data(1,:) = {'beam',1,size(mesh_struc,2), simu_fsi.radius_rbf,[1]};

if size(simu_fsi,2) ~= 0 
    if size(simu_aero,2) ~= 0 && size(simu_struct,2) ~= 0
        path_fsi   = [simu_fsi.currDir '\' simu_fsi.strfilename '\DeSiO'];
        path_aero  = [simu_aero.currDir '\' simu_aero.strfilename '\DeSiO-Aero'];
        path_struc = [simu_struct.currDir '\' simu_struct.strfilename '\DeSiO-Structure'];
        mkdir(path_fsi);
        filePattern = fullfile(path_aero,'*.txt'); dd = dir(filePattern);
        for j = 1:size(dd,1)
            if isempty(strfind(dd(j).name,'simulationinput'))
                copyfile([path_aero '\' dd(j).name], path_fsi);
            end
        end
        filePattern = fullfile([path_struc '\'],'*.txt'); dd = dir(filePattern);
        for j = 1:size(dd,1)
            if isempty(strfind(dd(j).name,'simulationinput'))
                copyfile([path_struc '\' dd(j).name], path_fsi);
            end
        end
    end

    % create fsi input files
    simu_fsi.simType = simType;
    simu_fsi.fsi_type = simfsitype;
    simu_fsi.flag_linearization = flag_linearization;
    simu_fsi.flag_output_matrix = flag_output_matrix;
    fun_writeDeSiOFSIInput(simu_fsi,simu_struct,simu_aero);

    copyfile([currDir '\' 'DeSiO.bat'],[currDir '\' strfilename '\DeSiO']);    
    copyfile([currDir '\' 'evaluate_disp_rot_fm.m'],[currDir '\' strfilename '\DeSiO']);    
end
close all;
return