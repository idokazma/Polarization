%fields only
clear Iz Iez Hx Hy Ez Im_x Im_y Ez_sol Hx_sol Hy_sol field_current temp_field_current_1

clear E_final_Fil E_final_inc_Fil E_final_Fil_center E_final_inc_Fil_center E_final_MoM E_final_inc_MoM
for rho = 1 : size(field_current_mat_output,1)
    for OMEGA_r = 1 : size(field_current_mat_output,2)
%         for sce = 7:9
        for sce = 1 :  size(field_current_mat_output,3)
%         for sce = 1 :  3
        field_current{rho,OMEGA_r,sce} = field_current_mat_output{rho,OMEGA_r,sce};
        temp_field_current_1 = cell2mat(field_current{rho,OMEGA_r,sce});

        
        Iz(sce) = temp_field_current_1(1);
        Iez(sce) = temp_field_current_1(2);
        Hx(sce) = temp_field_current_1(3);
        Hy(sce) = temp_field_current_1(4);
        Ez(sce) = temp_field_current_1(5);
        Im_x(sce) = temp_field_current_1(6);
        Im_y(sce) = temp_field_current_1(7);
        Ez_sol(sce) = temp_field_current_1(8);
        Hx_sol(sce) = temp_field_current_1(9);
        Hy_sol(sce) = temp_field_current_1(10);
        Ez_sol_center(sce) = temp_field_current_1(11);
        Hx_sol_center(sce) = temp_field_current_1(12);
        Hy_sol_center(sce) = temp_field_current_1(13);
        Ez_inc_center(sce) = temp_field_current_1(14);
        Hx_inc_center(sce) = temp_field_current_1(15);
        Hy_inc_center(sce) = temp_field_current_1(16);

        Ez_inc_mom(sce) = E_SOL_st_MoM_MEAN_PLANE_E0{rho,OMEGA_r,sce};
        Ez_sol_mom(sce) = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,sce};
%         Ez_sol_mom_center(sce) = E_SOL_st_alpha_MoM_center{rho,OMEGA_r,sce};
        Ez_sol_mom_center(sce) = 0;

        end
        
        E_final_Fil(rho,OMEGA_r,:) = (Ez_sol);
        E_final_inc_Fil(rho,OMEGA_r,:) = (Ez);
        E_final_Fil_center(rho,OMEGA_r,:) = (Ez_sol_center);
        E_final_inc_Fil_center(rho,OMEGA_r,:) = (Ez_inc_center);
        E_final_MoM(rho,OMEGA_r,:) = (Ez_sol_mom);
        E_final_inc_MoM(rho,OMEGA_r,:) = (Ez_inc_mom);
%         Ez_sol = Ez_sol+E_final_Fil(rho,4)-E_final_MoM(rho,4)
%         lineq_for_x = transpose([Hx_inc_center;Hy_inc_center ;Ez]);
%         lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; Ez_inc_center]);

%         lineq_for_x = transpose([1*Hx;1*Hy; Ez]);
        type_calc = '_inc_center_scattered_diff';

        lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; 1*Ez_inc_center]);
%         lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; 1*Ez_inc_mom]);

%            lineq_for_x = transpose([1*Hx;1*Hy;1*Ez_inc_center]);
%            lineq_for_x = transpose([1*Hx;1*Hy;1*Ez]);

        a=linsolve((lineq_for_x),transpose([Ez_sol-Ez]));
%         a=linsolve((lineq_for_x),transpose([Ez_sol-Ez]));
%         a=linsolve((lineq_for_x),transpose([Ez]));
%         a = lineq_for_x\transpose([Ez_sol]);
%         a = lsqminnorm((lineq_for_x),transpose([Ez_sol]), 1e-10);
        alpha_ezhx_center(rho,OMEGA_r) = a(1);
        alpha_ezhy_center(rho,OMEGA_r) = a(2);
        alpha_ezez_center(rho,OMEGA_r) = a(3);
        
%        lineq_for_x = transpose([1*Hx;1*Hy;1*Ez_inc_mom]);
%         lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; Ez_inc_center]);
        a=linsolve((lineq_for_x),transpose([Ez_sol_mom-Ez_inc_mom]));
%         a=linsolve((lineq_for_x),transpose([Ez_sol_mom]));
%         a=linsolve((lineq_for_x),transpose([Ez_inc_mom]));
        alpha_ezhx_mom_center(rho,OMEGA_r) = a(1);
        alpha_ezhy_mom_center(rho,OMEGA_r) = a(2);
        alpha_ezez_mom_center(rho,OMEGA_r) = a(3);

%         lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; Ez_inc_center]);
        a=linsolve((lineq_for_x),transpose([Ez_sol-Ez]));
%         a=linsolve((lineq_for_x),transpose([Ez_sol]));
%         a=linsolve((lineq_for_x),transpose([Ez]));
%         a = lineq_for_x\transpose([Ez_sol]);
%         a = lsqminnorm((lineq_for_x),transpose([Ez_sol]), 1e-10);
        alpha_ezhx_center_inc(rho,OMEGA_r) = a(1);
        alpha_ezhy_center_inc(rho,OMEGA_r) = a(2);
        alpha_ezez_center_inc(rho,OMEGA_r) = a(3);
        
%        lineq_for_x = transpose([1*Hx;1*Hy;1*Ez_inc_mom]);
%         lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; Ez_inc_center]);

        a=linsolve((lineq_for_x),transpose([Ez_sol_mom-Ez_inc_mom]));
%         a=linsolve((lineq_for_x),transpose([Ez_sol_mom]));
%         a=linsolve((lineq_for_x),transpose([Ez_inc_mom]));
        alpha_ezhx_mom_center_inc(rho,OMEGA_r) = a(1);
        alpha_ezhy_mom_center_inc(rho,OMEGA_r) = a(2);
        alpha_ezez_mom_center_inc(rho,OMEGA_r) = a(3);
        
        
        lineq_for_x = transpose([1*Hx;1*Hy;1*Ez]);

        a=linsolve((lineq_for_x),transpose([Ez_sol]));
%         a = lineq_for_x\transpose([Ez_sol]);
%         a = lsqminnorm((lineq_for_x),transpose([Ez_sol]), 1e-10);
        alpha_ezhx(rho,OMEGA_r) = a(1);
        alpha_ezhy(rho,OMEGA_r) = a(2);
        alpha_ezez(rho,OMEGA_r) = a(3);
        
%        lineq_for_x = transpose([1*Hx;1*Hy;1*Ez_inc_mom]);
%         lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; Ez_inc_center]);

        lineq_for_x = transpose([1*Hx;1*Hy;1*Ez_inc_mom]);

        a=linsolve((lineq_for_x),transpose([Ez_sol_mom]));
        alpha_ezhx_mom(rho,OMEGA_r) = a(1);
        alpha_ezhy_mom(rho,OMEGA_r) = a(2);
        alpha_ezez_mom(rho,OMEGA_r) = a(3);
        
        
        lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; Ez_inc_center]);
        a=linsolve((lineq_for_x),transpose([Ez_sol_center]));
%         a = lineq_for_x\transpose([Ez_sol]);
%         a = lsqminnorm((lineq_for_x),transpose([Ez_sol]), 1e-10);
        alpha_ezhx_center_center(rho,OMEGA_r) = a(1);
        alpha_ezhy_center_center(rho,OMEGA_r) = a(2);
        alpha_ezez_center_center(rho,OMEGA_r) = a(3);
        
        
        lineq_for_x = transpose([1*Hx_inc_center;1*Hy_inc_center; Ez_inc_center]);
        a=linsolve((lineq_for_x),transpose([Ez_sol_mom_center]));
%         a = lineq_for_x\transpose([Ez_sol]);
%         a = lsqminnorm((lineq_for_x),transpose([Ez_sol]), 1e-10);
        alpha_ezhx_mom_center_center(rho,OMEGA_r) = a(1);
        alpha_ezhy_mom_center_center(rho,OMEGA_r) = a(2);
        alpha_ezez_mom_center_center(rho,OMEGA_r) = a(3);
        
        a=linsolve((lineq_for_x),transpose([Hx_sol]));
        alpha_hxhx(rho,OMEGA_r) = a(1);
        alpha_hxhy(rho,OMEGA_r) = a(2);
        alpha_hxez(rho,OMEGA_r) = a(3);
        
        a=linsolve((lineq_for_x),transpose([Hy_sol]));
        alpha_hyhx(rho,OMEGA_r) = a(1);
        alpha_hyhy(rho,OMEGA_r) = a(2);
        alpha_hyez(rho,OMEGA_r) = a(3);
       

        lineq_for_x = transpose([Ez_inc_mom]);
        a=linsolve((lineq_for_x),transpose([Ez_sol_mom]));
        alpha_TM_MoM(rho,OMEGA_r) = a;
        
        lineq_for_x = transpose([Ez]);
        a=linsolve((lineq_for_x),transpose([Ez_sol]));
        alpha_TM(rho,OMEGA_r) = a;
        
    end
end

%%
alpha_ezez_center = alpha_ezez_center+1;
alpha_ezez_mom_center=alpha_ezez_mom_center+1;

% alpha_ezhx_center = alpha_ezhx_center+alpha_ezhx_center_inc;
% alpha_ezhx_mom_center=alpha_ezhx_mom_center+alpha_ezhx_mom_center_inc;
% 
% alpha_ezhy_center = alpha_ezhy_center+alpha_ezhy_center_inc;
% alpha_ezhy_mom_center=alpha_ezhy_mom_center+alpha_ezhy_mom_center_inc;

%% abs fil
figure('Name', 'Fil_abs_sigma_zz');
plot(params.OMEGA_vec/params.omega,abs(alpha_ezez_center) ,'LineWidth',2);
ylabel('$$|\sigma^{ee}_{zz}|$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['Fil: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'Fil_abs_sigma_zrho');
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx_center) ,'LineWidth',2);
ylabel('$$|\sigma^{em}_{z\rho}|$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['Fil: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'Fil_abs_sigma_ztheta');
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy_center) ,'LineWidth',2);
ylabel('$$|\sigma^{em}_{z\theta}|$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['Fil: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

%% phase fil
figure('Name', 'Fil_phase_sigma_zz');
plot(params.OMEGA_vec/params.omega,angle(alpha_ezez_center) ,'LineWidth',2);
ylabel('$$\angle \sigma^{ee}_{zz}$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['Fil: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'Fil_phase_sigma_zrho');
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx_center) ,'LineWidth',2);
ylabel('$$\angle \sigma^{em}_{z\rho}$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['Fil: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'Fil_phase_sigma_ztheta');
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy_center) ,'LineWidth',2);
ylabel('$$\angle\sigma^{em}_{z\theta}$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['Fil: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

%% abs mom
figure('Name', 'MoM_abs_sigma_zz');
plot(params.OMEGA_vec/params.omega,abs(alpha_ezez_mom_center) ,'LineWidth',2);
ylabel('$$|\sigma^{ee}_{zz}|$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['MoM: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'MoM_abs_sigma_zrho');
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx_mom_center) ,'LineWidth',2);
ylabel('$$|\sigma^{em}_{z\rho}|$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['MoM: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'MoM_abs_sigma_ztheta');
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy_mom_center) ,'LineWidth',2);
ylabel('$$|\sigma^{em}_{z\theta}|$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['MoM: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

%% phase mom
figure('Name', 'MoM_phase_sigma_zz');
plot(params.OMEGA_vec/params.omega,angle(alpha_ezez_mom_center) ,'LineWidth',2);
ylabel('$$\angle \sigma^{ee}_{zz}$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['MoM: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'MoM_phase_sigma_zrho');
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx_mom_center) ,'LineWidth',2);
ylabel('$$\angle \sigma^{em}_{z\rho}$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['MoM: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)

figure('Name', 'MoM_phase_sigma_ztheta');
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy_mom_center) ,'LineWidth',2);
ylabel('$$\angle\sigma^{em}_{z\theta}$$','Interpreter','Latex'); 
xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
title(['MoM: n_b=',num2str(params.n_out),', n_s=',num2str(params.n_in),', r=',num2str(params.original_radius/params.lambda),'\lambda, \lambda=',num2str(params.lambda*1e6), '\mum'],'FontSize',16)


figure('Name', '1A_Phase');
subplot(2,3,1); plot(params.OMEGA_vec/params.omega,angle(alpha_ezez_center) ,'LineWidth',2); title ('Fil angle($$\sigma^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,2); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy_center) ,'LineWidth',2); title ('Fil angle($$\sigma^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,3); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx_center) ,'LineWidth',2); title ('Fil angle($$\sigma^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

        
subplot(2,3,4); plot(params.OMEGA_vec/params.omega,angle(alpha_ezez_mom_center) ,'LineWidth',2); title ('MoM angle($$\sigma^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,5); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy_mom_center) ,'LineWidth',2); title ('MoM angle($$\sigma^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,6); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx_mom_center) ,'LineWidth',2); title ('MoM angle($$\sigma^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
       
lgd1 = legend (num2str(params.shift_vec'/params.lambda));
title(lgd1,'\rho_c [\lambda]')



figure('Name', '1A_Abs');
subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez_center) ,'LineWidth',2); title ('Fil abs($$\sigma^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy_center) ,'LineWidth',2); title ('Fil abs($$\sigma^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx_center) ,'LineWidth',2); title ('Fil abs($$\sigma^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

        
subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez_mom_center) ,'LineWidth',2); title ('MoM abs($$\sigma^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy_mom_center) ,'LineWidth',2); title ('MoM abs($$\sigma^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx_mom_center) ,'LineWidth',2); title ('MoM abs($$\sigma^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
       
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
%% 
if false
    now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
    if (params.n_in>params.n_out)
        FolderName = (['C:\Users\idoka\OneDrive\Desktop\Ido\Thesis_code\Polarization\' , now_str, '_pos',type_calc]);   % using my directory
        mkdir(FolderName)

    else
        FolderName = (['C:\Users\idoka\OneDrive\Desktop\Ido\Thesis_code\Polarization\' , now_str, '_neg',type_calc]);   % using my directory
        mkdir([FolderName])    
    end
        
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = FigHandle.Name;
      set(0, 'CurrentFigure', FigHandle);
      savefig(fullfile(FolderName, [FigName '.fig']));
      saveas(FigHandle,fullfile(FolderName, [FigName '.png']) );

    end
end