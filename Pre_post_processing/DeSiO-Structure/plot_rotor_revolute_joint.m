%% DeSiO -- Plotting of a Revolute Joint %%
% =================================================================================================================
% Script to plot the new rotor revolute joint constraint and test its
% properties.
% 
% Author:   Leon Minne
% Date:     17.04.2024
% Updated:  -
% =================================================================================================================

%% Preparation %%
% Path and file names
clc,clear,close all
path = 'C:\Users\leonm\Desktop\Bearbeitung\Revolute_Joint_Test_Blade\Test_Zwei_Verbindungsachsen\Revolute_Joint_Simulation_New_Offset_Test'; % path of files
addpath(genpath(path));

modelname = 'solution';

% Loading of DeSiO result files
time  = load([modelname '_t.dres']); % time
q     = load([modelname '_q.dres']); % generalized coordinates

n_vec = min([size(time,1),size(q,1)]); % modification the set vector length equal

time  = time(1:n_vec,:);
t     = time(:,1);
q     = q(1:n_vec,:);

% Parameters
node_a = 1; % Node A (the node to which the offset vectors are attached)
node_b = 2; % Node B (the node to rotate)

%phi_1  = [+4.72000000d+00,+0.00000000d+00,+7.40000000d-02]'; % Offset 1 for p_1 in Basis d_1a, d_2a, d_3a %[0,0,0]';%
%phi_2  = [-7.31073137d+00,+0.00000000d+00,+1.33848082d+00]'; % Offset 2 for p_2 in Basis d_1a, d_2a, d_3a %[0,0,0]';%

phi_1 = [1,1,1]';
phi_2 = [1,1,3]';

n_plot = n_vec; % maximum number of time steps to plot
nper   = 1;     % plot every nper % of total amount of time steps

%% Calculation %%
% Generalized coordinates and points with offsets
q_a  = q(:,12*node_a-11:12*node_a- 9);
d_1a = q(:,12*node_a- 8:12*node_a -6);
d_2a = q(:,12*node_a- 5:12*node_a -3);
d_3a = q(:,12*node_a- 2:12*node_a -0);

q_b  = q(:,12*node_b-11:12*node_b- 9);
d_1b = q(:,12*node_b- 8:12*node_b -6);
d_2b = q(:,12*node_b- 5:12*node_b -3);
d_3b = q(:,12*node_b- 2:12*node_b -0);

p_1  = q_a+phi_1(1)*d_1a+phi_1(2)*d_2a+phi_1(3)*d_3a;
p_2  = q_a+phi_2(1)*d_1a+phi_2(2)*d_2a+phi_2(3)*d_3a;

% Calculation of the axes
a = p_1-q_a;
b = p_2-q_a;
c = b-a;
d = q_b-p_1;
e = q_b-p_2;

% Axis lengths
lena = sqrt(a(:,1).^2+a(:,2).^2+a(:,3).^2);
lenb = sqrt(b(:,1).^2+b(:,2).^2+b(:,3).^2);
lenc = sqrt(c(:,1).^2+c(:,2).^2+c(:,3).^2);
lend = sqrt(d(:,1).^2+d(:,2).^2+d(:,3).^2);
lene = sqrt(e(:,1).^2+e(:,2).^2+e(:,3).^2);

% Dot products between directors A and axes
dotaa = [dot(a,d_1a,2),dot(a,d_2a,2),dot(a,d_3a,2)];
dotba = [dot(b,d_1a,2),dot(b,d_2a,2),dot(b,d_3a,2)];
dotca = [dot(c,d_1a,2),dot(c,d_2a,2),dot(c,d_3a,2)];
dotda = [dot(d,d_1a,2),dot(d,d_2a,2),dot(d,d_3a,2)];
dotea = [dot(e,d_1a,2),dot(e,d_2a,2),dot(e,d_3a,2)];

% Dot products between directors B and axes
dotab = [dot(a,d_1b,2),dot(a,d_2b,2),dot(a,d_3b,2)];
dotbb = [dot(b,d_1b,2),dot(b,d_2b,2),dot(b,d_3b,2)];
dotcb = [dot(c,d_1b,2),dot(c,d_2b,2),dot(c,d_3b,2)];
dotdb = [dot(d,d_1b,2),dot(d,d_2b,2),dot(d,d_3b,2)];
doteb = [dot(e,d_1b,2),dot(e,d_2b,2),dot(e,d_3b,2)];

%% Plotting %%
% Number of steps and difference in step numbers
dn     = round(nper/100*n_plot,0);

n_plot = max(min(n_plot,n_vec),1);
dn     = max(min(dn,n_vec),1);

% Plotting procedure
fig1 = figure('Name','Visualization of a Revolute Joint in the Current Configuration',...
       'NumberTitle','Off','Color',[0.9 0.9 0.9],...
       'Units','Normalized','Outerposition',[0 0 1 1]);
hold on

for i = 1:dn:n_plot
    % Plot directors
    quiver3(q_a(i,1),q_a(i,2),q_a(i,3),d_1a(i,1),d_1a(i,2),d_1a(i,3))
    quiver3(q_a(i,1),q_a(i,2),q_a(i,3),d_2a(i,1),d_2a(i,2),d_2a(i,3))
    quiver3(q_a(i,1),q_a(i,2),q_a(i,3),d_3a(i,1),d_3a(i,2),d_3a(i,3))
    quiver3(q_b(i,1),q_b(i,2),q_b(i,3),d_1b(i,1),d_1b(i,2),d_1b(i,3))
    quiver3(q_b(i,1),q_b(i,2),q_b(i,3),d_2b(i,1),d_2b(i,2),d_2b(i,3))
    quiver3(q_b(i,1),q_b(i,2),q_b(i,3),d_3b(i,1),d_3b(i,2),d_3b(i,3))
    
    % Plot axis and lines
    plot3([q_a(i,1),p_1(i,1),p_2(i,1),q_a(i,1)],...
        [q_a(i,2),p_1(i,2),p_2(i,2),q_a(i,2)],[q_a(i,3),p_1(i,3),p_2(i,3),q_a(i,3)])
    
    plot3([q_b(i,1),p_2(i,1),p_1(i,1),q_b(i,1)],...
        [q_b(i,2),p_2(i,2),p_1(i,2),q_b(i,2)],[q_b(i,3),p_2(i,3),p_1(i,3),q_b(i,3)])
end

view(45,15)
axis equal
grid minor
set(gca,'XMinorTick','On','YMinorTick','On','ZMinorTick','On')
title('Visualization of a Revolute Joint in the Current Configuration','FontSize',18,'Interpreter','LaTex')
xlabel('$x_1$','FontSize',18,'FontName','Times','Interpreter','LaTex')
ylabel('$x_2$','FontSize',18,'FontName','Times','Interpreter','LaTex')
zlabel('$x_3$','FontSize',18,'FontName','Times','Interpreter','LaTex')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','LaTex','FontName','Times')

fig2 = figure('Name','Visualization of a Revolute Joint in the Reference Configuration',...
       'NumberTitle','Off','Color',[0.9 0.9 0.9],...
       'Units','Normalized','Outerposition',[0 0 1 1]);
hold on

% Here, the displacement and rotation of node A is set to 0
for i = 1:dn:n_plot
    Ra = d_1a(1,:)'*d_1a(i,:)+d_2a(1,:)'*d_2a(i,:)+d_3a(1,:)'*d_3a(i,:);
    
    d_1an = d_1a(i,:)*Ra';
    d_2an = d_2a(i,:)*Ra';
    d_3an = d_3a(i,:)*Ra';
    
    d_1bn = d_1b(i,:)*Ra';
    d_2bn = d_2b(i,:)*Ra';
    d_3bn = d_3b(i,:)*Ra';
    
    q_an = q_a(1,:);
    q_bn = q_an+(q_b(i,:)-q_a(i,:))*Ra';
    
    p_1n = q_an+phi_1(1)*d_1an+phi_1(2)*d_2an+phi_1(3)*d_3an;
    p_2n = q_an+phi_2(1)*d_1an+phi_2(2)*d_2an+phi_2(3)*d_3an;
    
    % Plot directors
    quiver3(q_an(1,1),q_an(1,2),q_an(1,3),d_1an(1,1),d_1an(1,2),d_1an(1,3))
    quiver3(q_an(1,1),q_an(1,2),q_an(1,3),d_2an(1,1),d_2an(1,2),d_2an(1,3))
    quiver3(q_an(1,1),q_an(1,2),q_an(1,3),d_3an(1,1),d_3an(1,2),d_3an(1,3))
    quiver3(q_bn(1,1),q_bn(1,2),q_bn(1,3),d_1bn(1,1),d_1bn(1,2),d_1bn(1,3))
    quiver3(q_bn(1,1),q_bn(1,2),q_bn(1,3),d_2bn(1,1),d_2bn(1,2),d_2bn(1,3))
    quiver3(q_bn(1,1),q_bn(1,2),q_bn(1,3),d_3bn(1,1),d_3bn(1,2),d_3bn(1,3))
    
    % Plot axis and lines
    plot3([q_an(1,1),p_1n(1,1),p_2n(1,1),q_an(1,1)],...
        [q_an(1,2),p_1n(1,2),p_2n(1,2),q_an(1,2)],[q_an(1,3),p_1n(1,3),p_2n(1,3),q_an(1,3)])
    
    plot3([q_bn(1,1),p_2n(1,1),p_1n(1,1),q_bn(1,1)],...
        [q_bn(1,2),p_2n(1,2),p_1n(1,2),q_bn(1,2)],[q_bn(1,3),p_2n(1,3),p_1n(1,3),q_bn(1,3)])
end

view(45,15)
axis equal
grid minor
set(gca,'XMinorTick','On','YMinorTick','On','ZMinorTick','On')
title('Visualization of a Revolute Joint in the Reference Configuration','FontSize',18,'Interpreter','LaTex')
xlabel('$x_1$','FontSize',18,'FontName','Times','Interpreter','LaTex')
ylabel('$x_2$','FontSize',18,'FontName','Times','Interpreter','LaTex')
zlabel('$x_3$','FontSize',18,'FontName','Times','Interpreter','LaTex')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','LaTex','FontName','Times')

fig3 = figure('Name','Axes Lengths over Time',...
       'NumberTitle','Off','Color',[0.9 0.9 0.9],...
       'Units','Normalized','Outerposition',[0 0 1 1]);
hold on
plot(t,lena,'-k','LineWidth',1.5,'DisplayName','Axis A')
plot(t,lenb,'-r','LineWidth',1.5,'DisplayName','Axis B')
plot(t,lenc,'-m','LineWidth',1.5,'DisplayName','Axis C')
plot(t,lend,'-b','LineWidth',1.5,'DisplayName','Axis D')
plot(t,lene,'-g','LineWidth',1.5,'DisplayName','Axis E')
xlim([0,t(end)])
ylim([0,inf])
title('Axes Lengths over Time','FontSize',18,'Interpreter','LaTex')
xlabel('$t$ [s]','FontSize',18,'FontName','Times','Interpreter','LaTex')
ylabel('$\|\mathbf{v}(t)\|_2$ [m]','FontSize',18,'FontName','Times','Interpreter','LaTex')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','LaTex','FontName','Times')
legend('FontSize',16,'FontName','Times','Interpreter','LaTex','Location','EastOutside')

fig4 = figure('Name','Dot Products Node A over Time',...
       'NumberTitle','Off','Color',[0.9 0.9 0.9],...
       'Units','Normalized','Outerposition',[0 0 1 1]);
hold on
plot(t,dotaa,'-k','LineWidth',1.5,'DisplayName','Axis A')
plot(t,dotba,'-r','LineWidth',1.5,'DisplayName','Axis B')
plot(t,dotca,'-m','LineWidth',1.5,'DisplayName','Axis C')
%plot(t,dotda,'-b','LineWidth',1.5,'DisplayName','Axis D') % not constrained
%plot(t,dotea,'-g','LineWidth',1.5,'DisplayName','Axis E') % not constrained
xlim([0,t(end)])
ylim([0,inf])
title('Dot Products Node A over Time','FontSize',18,'Interpreter','LaTex')
xlabel('$t$ [s]','FontSize',18,'FontName','Times','Interpreter','LaTex')
ylabel('$\mathbf{d}_i^A(t)\cdot\mathbf{v}(t)$ [m]','FontSize',18,'FontName','Times','Interpreter','LaTex')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','LaTex','FontName','Times')
legend('FontSize',16,'FontName','Times','Interpreter','LaTex','Location','EastOutside')

fig5 = figure('Name','Dot Products Node B over Time',...
       'NumberTitle','Off','Color',[0.9 0.9 0.9],...
       'Units','Normalized','Outerposition',[0 0 1 1]);
hold on
%plot(t,dotab,'-k','LineWidth',1.5,'DisplayName','Axis A') % not constrained
%plot(t,dotbb,'-r','LineWidth',1.5,'DisplayName','Axis B') % not constrained
plot(t,dotcb,'-m','LineWidth',1.5,'DisplayName','Axis C')
plot(t,dotdb,'-b','LineWidth',1.5,'DisplayName','Axis D')
plot(t,doteb,'-g','LineWidth',1.5,'DisplayName','Axis E')
xlim([0,t(end)])
ylim([0,inf])
title('Dot Products Node B over Time','FontSize',18,'Interpreter','LaTex')
xlabel('$t$ [s]','FontSize',18,'FontName','Times','Interpreter','LaTex')
ylabel('$\mathbf{d}_i^B(t)\cdot\mathbf{v}(t)$ [m]','FontSize',18,'FontName','Times','Interpreter','LaTex')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','LaTex','FontName','Times')
legend('FontSize',16,'FontName','Times','Interpreter','LaTex','Location','EastOutside')
