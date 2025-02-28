function uvlm_plot_CL_time()
% =================================================================================================================
    clc; clear all; close all;
    addpath('..\..\..\Pre_post_processing\DeSiO-Aero\');
    
    currDir = cd;
    strfilename = 'ex01d'
    
    % initializing figures
    f1 = figure(); hold on; grid on; xlabel('time s'); ylabel('C_L');
    % calculate and plot lift coeficient
    model = uvlm_readmodel();
    model.windrotax = [1;0;0];
    t   = fun_load_file([model.strSimName '_uvlm_t.dres']); t = t(:,1);
    dp  = fun_load_file([model.strSimName '_uvlm_dps_cp.dres']);
    qs  = fun_load_file([model.strSimName '_uvlm_qs_nodal.dres']);
    [CL] = uvlm_CL_surface_i(model,t,qs,dp);
    plot(t,CL(:,1),'.-','linewidth',2,'markersize',20,'color','b');
    plot(t,CL(:,2),'.-','linewidth',2,'markersize',20,'color','r');
    set(gca,'fontweight','bold','fontsize',12)
    print([strfilename '_CL_vs_time'],'-dpng', '-r500');
    close all;
% =================================================================================================================
return