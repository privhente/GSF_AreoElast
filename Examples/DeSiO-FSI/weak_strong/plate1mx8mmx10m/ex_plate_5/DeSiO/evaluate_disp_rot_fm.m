function evaluate_disp_rot_fm()
% =========================================================================
% Plotting position vector
% =========================================================================
clc;
clear all;
close all

addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));
currDir = cd;

% loading of DeSiO result files
model  = fsi_readmodel;
q      = load([model.strSimName '_q.dres']);
time   = load([model.strSimName '_t.dres']);
lambda = load([model.strSimName '_lambda.dres']);

inz = min(size(q,1),size(time,1));
time = time(1:inz,:);

node     = 25;
deltat   = time(2,1)-time(1,1);
strCase  = 'ex02';

markersize = 4;
linewidth  = 2;
fontsize   = 20;
fontname   = 'times';
fontweight = 'normal';

% determine rotations of beam nodes
cd(currDir);
for i = 1:length(node)
    fig_u = figure(); hold on; grid off;
    [u] = get_DeSiO_dof_solu(node(i),1,q);
    plot(time(:,1),u(:,1),'linestyle','-','lineWidth',linewidth,'color','r');
    plot(time(:,1),u(:,2),'linestyle','-','lineWidth',linewidth,'color','g');
    plot(time(:,1),u(:,3),'linestyle','-','lineWidth',linewidth,'color','b');
    str_leg = {'u_1','u_2','u_3'};
    ylabel('u in ft'); 
    xlabel('step');
    figure(fig_u); 
    set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
    legend(str_leg,'location','best','orientation','horizontal','fontsize',fontsize,'box','off');
    fig_u.Position = [100 100 1200 400]; pbaspect auto;
    savefig(fig_u,[strCase '_disp' '.fig']);
    print([strCase '_disp'],'-dpng', '-r500');    
    
%     fig_phi = figure(); hold on; grid off;
%     d1_0 = q(1,12*(node(i)-1)+4  : 12*(node(i)-1)+6);
%     d2_0 = q(1,12*(node(i)-1)+7  : 12*(node(i)-1)+9);
%     d3_0 = q(1,12*(node(i)-1)+10 : 12*(node(i)-1)+12);
%     phi = zeros(length(time),3);
%     for j = 2:length(time)
%         d1_n = q(j,12*(node(i)-1)+4  : 12*(node(i)-1)+6);
%         d2_n = q(j,12*(node(i)-1)+7  : 12*(node(i)-1)+9);
%         d3_n = q(j,12*(node(i)-1)+10 : 12*(node(i)-1)+12);
%         delta_d1 = d1_n - d1_0; delta_d2 = d2_n - d2_0; delta_d3 = d3_n - d3_0;
%         % global rotation increments
%         phi(j,1:3) = phi(j-1,1:3) + 0.5*( cross(d1_n,delta_d1) + cross(d2_n,delta_d2) + cross(d3_n,delta_d3) );
%         d1_0 = d1_n; d2_0 = d2_n; d3_0 = d3_n;
%     end
%     plot(time(:,1),phi(:,1).*180/pi,'linestyle','-','lineWidth',linewidth,'color','r');
%     plot(time(:,1),phi(:,2).*180/pi,'linestyle','-','lineWidth',linewidth,'color','g');
%     plot(time(:,1),phi(:,3).*180/pi,'linestyle','-','lineWidth',linewidth,'color','b');
%     str_leg = {'\phi_1','\phi_2','\phi_3'};
%     ylabel('\phi in ^\circ'); 
%     xlabel('step')
%     figure(fig_phi); 
%     set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
%     legend(str_leg,'location','best','orientation','horizontal','fontsize',fontsize,'box','off');
%     fig_phi.Position = [100 100 1200 400]; pbaspect auto;
%     savefig(fig_phi,[strCase '_rot' '.fig']);
%     print([strCase '_rot'],'-dpng', '-r500');    
%     
%     % gradient of rotation
%     fig_gradphi = figure(); hold on; grid off;
%     grad_phi = zeros(length(time),3);
%     for j = 2:size(time,1)-1
%         grad_phi(j,:) = (phi(j+1,:)-phi(j,:))/deltat;
%     end
%     plot(time(:,1),grad_phi(:,1),'linestyle','-','lineWidth',linewidth,'color','r');
%     plot(time(:,1),grad_phi(:,2),'linestyle','-','lineWidth',linewidth,'color','g');
%     plot(time(:,1),grad_phi(:,3),'linestyle','-','lineWidth',linewidth,'color','b');
%     str_leg = {'\nabla phi_1','\nabla phi_2','\nabla phi_3'};
%     ylabel('\nabla phi in ^\circ'); 
%     xlabel('step')
%     set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
%     legend(str_leg,'location','best','orientation','horizontal','fontsize',fontsize,'box','off');
%     fig_gradphi.Position = [100 100 1200 400]; pbaspect auto;
%     savefig(fig_gradphi,[strCase '_gradphi' '.fig']);
%     print([strCase '_gradphi'],'-dpng', '-r500');    
%     
%     % torsion moment at fixed side to calculate torsional stiffness
%     fig_moment = figure(); hold on; grid off; 
%     for j = 1:size(time,1)
%         q_node = q(j,12*(node(i)-1)+1:12*(node(i)-1)+12);
%         f_node = lambda(j,[1:12]);
%         m_node(j,:) = -( cross(f_node(4:6),q_node(4:6)) + cross(f_node(7:9),q_node(7:9)) + cross(f_node(10:12),q_node(10:12)) );
%     end
%     plot(time(:,1),m_node(:,1),'linestyle','-','lineWidth',linewidth,'color','r');
%     plot(time(:,1),m_node(:,2),'linestyle','-','lineWidth',linewidth,'color','g');
%     plot(time(:,1),m_node(:,3),'linestyle','-','lineWidth',linewidth,'color','b');
%     str_leg = {'M_{1,r}','M_{2,r}','M_{3,r}'};
%     ylabel('M in lbf ft');
%     xlabel('step')
%     figure(fig_moment); 
%     set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
%     legend(str_leg,'location','best','orientation','horizontal','fontsize',fontsize,'box','off');
%     fig_moment.Position = [100 100 1200 400]; pbaspect auto;
%     savefig(fig_moment,[strCase '_moments' '.fig']);
%     print([strCase '_moments'],'-dpng', '-r500');    
%     
%     % torsion moment at fixed side to calculate torsional stiffness
%     fig_force = figure(); hold on; grid off; 
%     f_node = [];
%     for j = 1:size(time,1)
%         f_node(j,1:3) = lambda(j,1:3);
%     end
%     plot(time(:,1),f_node(:,1),'linestyle','-','lineWidth',linewidth,'color','r');
%     plot(time(:,1),f_node(:,2),'linestyle','-','lineWidth',linewidth,'color','g');
%     plot(time(:,1),f_node(:,3),'linestyle','-','lineWidth',linewidth,'color','b');
%     str_leg = {'F_{1,r}','F_{2,r}','F_{3,r}'};
%     ylabel('F in lbf'); 
%     xlabel('step')
%     figure(fig_force); 
%     set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
%     legend(str_leg,'location','best','orientation','horizontal','fontsize',fontsize,'box','off');
%     fig_force.Position = [100 100 1200 400]; pbaspect auto;
%     savefig(fig_force,[strCase '_forces' '.fig']);
%     print([strCase '_forces'],'-dpng', '-r500');     
end
return