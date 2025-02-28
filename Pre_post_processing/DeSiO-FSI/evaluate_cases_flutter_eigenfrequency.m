function evaluate_cases_flutter_eigenfrequency()
% =========================================================================
% Plotting position vector
% =========================================================================
clc;
clear all;
close all;

addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));
currDir = cd;

flag_delete   = 1;
flag_plot_ev  = 1;
nev           = [1:5];
flag_legend   = 1;
flag_restart  = 1;

strcase  = 'flutter_eigenproblem';
fontname = 'times';
fontsize = 14;
fontweight = 'normal';

% Reading solution files
file = dir(currDir);
m    = []; arrfolder = {};
for i = 1:size(file,1)
    if file(i).isdir == 1 && length(file(i).name)>2
        folder   = file(i).name;
        i_1 = findstr(folder,'v');
        i_1 = i_1(1);
        i_2 = findstr(folder,'m');
        if i_1 ~= 0
            var1 =  str2num(folder(i_1+1:i_2(1)-1));
            m(end+1,1) = [var1];
            arrfolder{end+1} = folder;
        end
    end
end
[val, inz] = sort(m(:,1)); 
m = m(inz(1:end),:); 
arrlength = m;
%
for i = 1:length(m)
   folder_name{i,:} = arrfolder{inz(i)};
end

rgb_Color = RGB_Color;

str_leg   = {};
res       = [];
arr_ev    = [];
arr_ev_rs = [];
arr_ev_is = [];
arr_ev_st = [];
arr_ev_st_d0 = [];
arr_vinf  = [];
if flag_restart == 1
    for i = 1:length(m)
        
        flag_evaluate = 0;
        cd([currDir '\' folder_name{i} '\DeSiO']);
        model_fsi = fsi_readmodel;
    
        if exist(['res' num2str(m(i)) '.mat'], 'file') && flag_delete == 0
            res = load(['res' num2str(m(i))],'-mat');
            res = res.res;
            t = load([model_fsi.strSimName '_t.dres']); t = t(:,1);
        else
            flag_evaluate = 0;
            if exist(['check.log'], 'file')
                copyfile([currDir '/' 'evaluate_for_flutter_eigenvalue_analysis.m'], ['evaluate_for_flutter_eigenvalue_analysis.m']);
                [e_rs, e_is, e_s, e_s_d0, Xs] = evaluate_for_flutter_eigenvalue_analysis();
                flag_evaluate = 1;
            end
        end
        
        if flag_evaluate == 1
            % plot time history diagrams
            if flag_plot_ev == 1
                % eigenvalues
    %             figure(fig_ev); plot([1:nev],e_is(1:nev),'color',rgb_Color(i,:),'lineWidth',2,'linestyle','--','marker','.','markersize',18);
            end
    
            str_leg{end+1}  = ['vinf' num2str(m(i))];
            arr_vinf(end+1) = m(i);
            arr_ev_rs(end+1,:) = e_rs(1:50);%./(2*pi);
            arr_ev_is(end+1,:) = e_is(1:50);%./(2*pi);
            arr_ev_st(end+1,:) = e_s(1:50);%./(2*pi);
            arr_ev_st_d0(end+1,:) = e_s_d0(1:50);%./(2*pi);
        end
    end
    cd([currDir])
    save('arr_vinf','arr_vinf');
    save('arr_ev_rs','arr_ev_rs');
    save('arr_ev_is','arr_ev_is');
    save('arr_ev_st','arr_ev_st');
    save('arr_ev_st_d0','arr_ev_st_d0');
else
    arr_vinf = load('arr_vinf');   arr_vinf  = arr_vinf.arr_vinf;
    arr_ev_rs = load('arr_ev_rs'); arr_ev_rs = arr_ev_rs.arr_ev_rs;
    arr_ev_is = load('arr_ev_is'); arr_ev_is = arr_ev_is.arr_ev_is;
    arr_ev_st = load('arr_ev_st'); arr_ev_st = arr_ev_st.arr_ev_st;
    arr_ev_st_d0 = load('arr_ev_st_d0'); arr_ev_st_d0 = arr_ev_st_d0.arr_ev_st_d0;    
end
arr_ev_rs = arr_ev_rs(:,nev);
arr_ev_is = arr_ev_is(:,nev);
arr_ev_st = arr_ev_st(:,nev);
arr_ev_st_d0 = arr_ev_st_d0(:,nev);

% arr_vinf(4)  = [];
% arr_ev_rs(4,:) = [];
% arr_ev_is(4,:) = [];
% arr_ev_st(4,:) = [];

cd(currDir);

if flag_plot_ev == 1
    % find the transition from negativ to positive reaql eigenvalues
    arr_p_crit = [];
    for j = 1:size(arr_ev_rs,2)
        if arr_ev_rs(1,j)>=0
            arr_ev_rs(:,j) = arr_ev_rs(:,j)*(-1);
        end
        inz1_sign_changes = find(diff(sign(arr_ev_rs(:,j))) ~= 0);
        if ~isempty(inz1_sign_changes)
            for k = 1:length(inz1_sign_changes)
                [xi,yi]=polyxpoly([arr_vinf(inz1_sign_changes(k)) arr_vinf(inz1_sign_changes(k)+1)],[arr_ev_rs(inz1_sign_changes(k),j), arr_ev_rs(inz1_sign_changes(k)+1,j)], [arr_vinf],[arr_vinf*0]);
                arr_p_crit(end+1:end+length(xi),1:2) = [xi,yi];
            end
        end
    end
    
    % find intersection between beding and torsion eigenvalue
    omega = arr_ev_rs + sqrt(-1)*arr_ev_is;
    omega = [omega];
    
    fig_ev_vel = figure(); 
    subplot(2,3,[1,2,3]);
    hold on; grid on; 
    ylabel('\omega_r in rad/s'); %xlabel('V_{\infty} ft/s');
    set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname,'box','on','XTick',[120:10:170],'YTick',[-0.6:0.1:0.2]);
    pbaspect auto;
    strleg = {};
    for i = 1:length(nev)
        plot(arr_vinf,arr_ev_rs(:,i),'lineWidth',2,'marker','.','markersize',16);
        strleg(end+1) = {['ev' num2str(i)]};
    end
    if ~isempty(arr_p_crit)
        plot(arr_p_crit(:,1),arr_p_crit(:,2),'x','markersize',16,'color','k');
    end
    if flag_legend==1
%         leg = legend(strleg);
%         set(leg,'location','best','orientation','vertical','fontsize',fontsize,'box','off');
    end
    
    subplot(2,3,[4,5,6]);
    hold on; grid on; box on;
    ylabel('\omega_i in rad/s'); xlabel('v_{\infty} ft/s');
    set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname,'box','on','XTick',[120:10:170],'YTick',[0:1:6]);
    strleg = {};
    for i = 1:length(nev)
        plot(arr_vinf,arr_ev_is(:,i),'lineWidth',2,'marker','.','markersize',20);
        strleg(end+1) = {['ev' num2str(i)]};
    end
    if flag_legend==1
        strleg={'1st IPB', '1st OPB','1st torsional','2nd torsional','2nd IPB','2nd OPB'};
        leg = legend(strleg);
%         set(leg,'location','best','orientation','vertical','fontsize',fontsize,'box','off');
        set(leg,'location','northoutside','orientation','horizontal','fontsize',fontsize,'box','off');
    end

    figure();
    hold on; grid on; 
    ylabel('\omega_s in rad/s'); xlabel('V_{\infty} ft/s');
    set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
    strleg = {};
    for i = 1:length(nev)
        plot(arr_vinf,arr_ev_st(:,i),'lineWidth',2,'marker','.','markersize',20);
        strleg(end+1) = {['ev' num2str(i)]};
    end
    if flag_legend==1
        leg = legend(strleg);
        set(leg,'location','best','orientation','vertical','fontsize',fontsize,'box','off');
    end

    figure();
    hold on; grid on; 
    ylabel('\omega_s d0 in rad/s'); xlabel('V_{\infty} ft/s');
    set(gca,'fontWeight',fontweight,'fontsize',fontsize,'fontname',fontname);
    strleg = {};
    for i = 1:length(nev)
        plot(arr_vinf,arr_ev_st_d0(:,i),'lineWidth',2,'marker','.','markersize',20);
        strleg(end+1) = {['ev' num2str(i)]};
    end
    if flag_legend==1
        leg = legend(strleg);
        set(leg,'location','best','orientation','vertical','fontsize',fontsize,'box','off');
    end

    fig_ev_vel.Position = [100 100 1200 800]; pbaspect auto;
    savefig(fig_ev_vel,[strcase '_ev_vs_vinf' '.fig']);
    print([strcase '_ev_vs_vinf'],'-dpng', '-r500');
end
return