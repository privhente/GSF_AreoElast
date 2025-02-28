function res = uvlm_get_m_f_time(x0)
% =================================================================================================================
    addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));
    
    % calculate and plot lift coeficient
    model = uvlm_readmodel();
    model.windrotax = [0;1;0];
    t   = fun_load_file([model.strSimName '_uvlm_t.dres']); t = t(:,1);
    dp  = fun_load_file([model.strSimName '_uvlm_dps_cp.dres']);
    qs  = fun_load_file([model.strSimName '_uvlm_qs_nodal.dres']);

    M = uvlm_m_surface(model,t,qs,dp,x0);
    F = uvlm_f_surface(model,t,qs,dp);
    
    res.M = M;
    res.F = F;
    save('res','res');
% =================================================================================================================return
return