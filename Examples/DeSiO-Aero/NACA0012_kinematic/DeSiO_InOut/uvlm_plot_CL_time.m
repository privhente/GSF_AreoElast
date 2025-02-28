function uvlm_plot_CL_time()
% =================================================================================================================
    clc; clear all; close all;
    addpath('..\..\..\Pre_post_processing\DeSiO-Aero\');
    
    currDir = cd;
    strfilename = 'ex01b'
    
    % initializing figures
    f1 = figure(); hold on; grid on; xlabel('time s'); ylabel('C_L');
    % calculate and plot lift coeficient
    model = uvlm_readmodel();
    model.windrotax = [1;0;0];
    t   = fun_load_file([model.strSimName '_uvlm_t.dres']); t = t(:,1);
    dp  = fun_load_file([model.strSimName '_uvlm_dps_cp.dres']);
    qs  = fun_load_file([model.strSimName '_uvlm_qs_nodal.dres']);
    vref = 1*[0.00000, 9.848077530122080e-01, 1.736481776669303e-01];
    [CL] = uvlm_CL_surface(model,t,qs,dp,vref);
    plot(t,CL,'.-','linewidth',2,'markersize',20);
    set(gca,'fontweight','bold','fontsize',12)
    print([strfilename '_CL_vs_time'],'-dpng', '-r500');
    close all;
% =================================================================================================================
return