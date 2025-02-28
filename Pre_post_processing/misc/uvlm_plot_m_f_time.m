function uvlm_plot_m_f_time()
% omega in rpm, theta in degree, X0 as vector
% =================================================================================================================
    str_software_path = 'C:\Users\hente\Desktop';

    strpath = fileread([str_software_path '\DeSiO_path.yaml']);
    strpath = textscan(strpath,'%s');
    addpath(genpath(strpath{1}{8}));

    theta   = 0*pi/180;
    omega   = 7.25*2*pi/60;
    X0      = [0;0;0];
    n_rotor = [cos(theta);0;-sin(theta)];
    
    % calculate and plot lift coeficient
    model = uvlm_readmodel();
    model.windrotax = [0;1;0];
    t   = fun_load_file([model.strSimName '_uvlm_t.dres']); t = t(:,1);
    dp  = fun_load_file([model.strSimName '_uvlm_dps_cp.dres']);
    qs  = fun_load_file([model.strSimName '_uvlm_qs_nodal.dres']);

    M = uvlm_m_surface(model,t,qs,dp,X0);
    F = uvlm_f_surface(model,t,qs,dp);
    
    m_aero = n_rotor'*M';
    f_aero = n_rotor'*F';
    
    % initializing figures
    fig = figure(); hold on; grid on; 
    xlabel('time s'); 
    yyaxis left; ylabel('M_{rotor} in MNm');
    plot(t,m_aero/1e6,'-','linewidth',2,'markersize',20,'color','b');
    yyaxis right; ylabel('P_{rotor} in MW');
    plot(t,m_aero*omega/1e6,'-','linewidth',2,'markersize',20,'color','r');
    legend('M_{rotor}','P_{rotor}');
    set(gca,'fontweight','bold','fontsize',12)
    savefig(fig,'moment.fig');

    fig = figure(); hold on; grid on; 
    xlabel('time s'); 
    ylabel('F_{thrust} in MN');
    plot(t,f_aero/1e6,'-','linewidth',2,'markersize',20,'color','b');
    legend('F_{thrust}');
    set(gca,'fontweight','bold','fontsize',12)
    savefig(fig,'f_thrust.fig');
    
% =================================================================================================================return

return