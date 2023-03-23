%fields only
clear Iz Iez Hx Hy Ez Im_x Im_y Ez_sol Hx_sol Hy_sol field_current temp_field_current_1
for rho = 1 : size(field_current_mat_output,1)
    for OMEGA_r = 1 : size(field_current_mat_output,2)
        for sce = 1 :  size(field_current_mat_output,3)
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


        Ez_inc_mom(sce) = E_SOL_st_MoM_MEAN_PLANE_E0{rho,OMEGA_r,sce};
        Ez_sol_mom(sce) = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,sce};

        end
        
 
        
        lineq_for_x = transpose([Hx ; Hy ;Ez]);
        a=linsolve((lineq_for_x),transpose([Ez_sol]));
        alpha_ezhx(rho,OMEGA_r) = a(1);
        alpha_ezhy(rho,OMEGA_r) = a(2);
        alpha_ezez(rho,OMEGA_r) = a(3);
        
        lineq_for_x = transpose([Hx ; Hy ;Ez_inc_mom]);
        a=linsolve((lineq_for_x),transpose([Ez_sol_mom]));
        alpha_ezhx_mom(rho,OMEGA_r) = a(1);
        alpha_ezhy_mom(rho,OMEGA_r) = a(2);
        alpha_ezez_mom(rho,OMEGA_r) = a(3);
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

for i = 1:length(alpha_ezez)
    for k=1:length(alpha_ezez(:,1))
   
        if params.OMEGA_vec(i)==0
            alpha_ezez(k,i) = abs(alpha_ezez(k,i))*exp(1i*(angle(alpha_ezez(k,i-1) )+angle(alpha_ezez(k,i+1)))/2) ; 
            alpha_ezhy(k,i) = abs(alpha_ezhy(k,i))*exp(1i*(angle(alpha_ezhy(k,i-1) )+angle(alpha_ezhy(k,i+1)))/2) ;
            alpha_ezhx(k,i) = abs(alpha_ezhx(k,i))*exp(1i*(angle(alpha_ezhx(k,i-1) )+angle(alpha_ezhx(k,i+1)))/2) ;
            alpha_hxez(k,i) = abs(alpha_hxez(k,i))*exp(1i*(angle(alpha_hxez(k,i-1) )+angle(alpha_hxez(k,i+1)))/2) ;
            alpha_hyez(k,i) = abs(alpha_hyez(k,i))*exp(1i*(angle(alpha_hyez(k,i-1) )+angle(alpha_hyez(k,i+1)))/2) ;
            alpha_hxhx(k,i) = abs(alpha_hxhx(k,i))*exp(1i*(angle(alpha_hxhx(k,i-1) )+angle(alpha_hxhx(k,i+1)))/2) ;
            alpha_hxhy(k,i) = abs(alpha_hxhy(k,i))*exp(1i*(angle(alpha_hxhy(k,i-1) )+angle(alpha_hxhy(k,i+1)))/2) ;
            alpha_hyhx(k,i) = abs(alpha_hyhx(k,i))*exp(1i*(angle(alpha_hyhx(k,i-1) )+angle(alpha_hyhx(k,i+1)))/2) ;
            alpha_hyhy(k,i) = abs(alpha_hyhy(k,i))*exp(1i*(angle(alpha_hyhy(k,i-1) )+angle(alpha_hyhy(k,i+1)))/2) ;
        end
    end
end

figure;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,angle(alpha_ezez),'x-','LineWidth',2); title ('$$\angle\sigma^{ee}_{zz}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy),'x-','LineWidth',2); title ('$$\angle\sigma^{em}_{z\theta}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx),'x-','LineWidth',2); title ('$$\angle\sigma^{em}_{z\rho}$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,angle(alpha_hxez),'x-','LineWidth',2); title ('$$\angle\sigma^{me}_{\rho z}$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,angle(alpha_hyez),'x-','LineWidth',2); title ('$$\angle\sigma^{me}_{\theta z}$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,angle(alpha_hxhx),'x-','LineWidth',2); title ('$$\angle\sigma^{mm}_{\rho\rho}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,angle(alpha_hxhy),'x-','LineWidth',2); title ('$$\angle\sigma^{mm}_{\rho\theta}$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,angle(alpha_hyhx),'x-','LineWidth',2); title ('$$\angle\sigma^{mm}_{\theta\rho}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,angle(alpha_hyhy),'x-','LineWidth',2); title ('$$\angle\sigma^{mm}_{\theta\theta}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')



figure;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('$$|\sigma^{ee}_{zz}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('$$|\sigma^{em}_{z\theta}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('$$|\sigma^{em}_{z\rho}|$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hxez),'x-','LineWidth',2); title ('$$|\sigma^{me}_{\rho z}|$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_hyez),'x-','LineWidth',2); title ('$$|\sigma^{me}_{\theta z}|$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhx),'x-','LineWidth',2); title ('$$|\sigma^{mm}_{\rho\rho}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhy),'x-','LineWidth',2); title ('$$|\sigma^{mm}_{\rho\theta}|$$','Interpreter','Latex', 'FontSize', 16);xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhx),'x-','LineWidth',2); title ('$$|\sigma^{mm}_{\theta\rho}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhy),'x-','LineWidth',2); title ('$$|\sigma^{mm}_{\theta\theta}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('$$|\sigma^{ee}_{zz}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
figure;
plot(params.OMEGA_vec/params.omega,angle(alpha_ezez),'x-','LineWidth',2); title ('$$\angle\sigma^{ee}_{zz}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('$$|\sigma^{em}_{z\rho}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
figure;
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx),'x-','LineWidth',2); title ('$$\angle\sigma^{em}_{z\rho}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('$$|\sigma^{em}_{z\theta}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
figure;
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy),'x-','LineWidth',2); title ('$$\angle\sigma^{em}_{z\theta}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('abs($$\sigma^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('abs($$\sigma^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('abs($$\sigma^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hxez),'x-','LineWidth',2); title ('abs($$\sigma^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_hyez),'x-','LineWidth',2); title ('abs($$\sigma^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhx),'x-','LineWidth',2); title ('abs($$\sigma^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhy),'x-','LineWidth',2); title ('abs($$\sigma^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhx),'x-','LineWidth',2); title ('abs($$\sigma^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhy),'x-','LineWidth',2); title ('abs($$\sigma^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hxez),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_hyez),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhx),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhy),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhx),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhy),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_ezez),'x-','LineWidth',2); title ('real($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_ezhy),'x-','LineWidth',2); title ('real($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_ezhx),'x-','LineWidth',2); title ('real($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_hxez),'x-','LineWidth',2); title ('real($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,real(alpha_hyez),'x-','LineWidth',2); title ('real($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_hxhx),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hxhy),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,real(alpha_hyhx),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,real(alpha_hyhy),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,angle(alpha_ezez),'x-','LineWidth',2); title ('phase($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy),'x-','LineWidth',2); title ('phase($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx),'x-','LineWidth',2); title ('phase($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,angle(alpha_hxez),'x-','LineWidth',2); title ('phase($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,angle(alpha_hyez),'x-','LineWidth',2); title ('phase($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,angle(alpha_hxhx),'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,angle(alpha_hxhy),'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,angle(alpha_hyhx),'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,angle(alpha_hyhy),'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')


figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,real(alpha_ezez),'x-','LineWidth',2); title ('real($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,real(alpha_ezhy),'x-','LineWidth',2); title ('real($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,real(alpha_ezhx),'x-','LineWidth',2); title ('real($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_hxez),'x-','LineWidth',2); title ('real($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hyez),'x-','LineWidth',2); title ('real($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_hxhx),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_hxhy),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_hyhx),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_hyhy),'x-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,real(alpha_ezez),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,real(alpha_ezhy),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,real(alpha_ezhx),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_hxez),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hyez),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_hxhx),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_hxhy),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_hyhx),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_hyhy),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,imag(alpha_ezez),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,imag(alpha_ezhy),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,imag(alpha_ezhx),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_hxez),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,imag(alpha_hyez),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_hxhx),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_hxhy),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_hyhx),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_hyhy),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,angle(alpha_ezez)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,angle(alpha_hxez)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,angle(alpha_hyez)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,angle(alpha_hxhx)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,angle(alpha_hxhy)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,angle(alpha_hyhx)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,angle(alpha_hyhy)/pi*180,'x-','LineWidth',2); title ('phase($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')


figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez-alpha_ezez(1)),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_hxez),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hyez),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhx),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhy),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhx),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhy),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
plot(params.OMEGA_vec/params.omega, abs(alpha_TM),'x-'); hold on;
figure;
plot(params.OMEGA_vec/params.omega, abs(alpha_TM_MoM),'d-');
grid on;

figure;
subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez-alpha_ezez(1,:)),'x-','LineWidth',2); title ('Fil abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('Fil abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('Fil abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;


        
subplot(2,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez_mom-alpha_ezez_mom(1,:)),'x-','LineWidth',2); title ('MoM abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
       
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')


figure;
subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('Fil abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('Fil abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('Fil abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;


        
subplot(2,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
       
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')


figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_ezjz),'x-','LineWidth',2); title ('abs($$\alpha^{ez}_{jz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_hyjz),'x-','LineWidth',2); title ('abs($$\alpha^{h\theta}_{jz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_hxjz),'x-','LineWidth',2); title ('abs($$\alpha^{h\rho}_{jz}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_ezjy),'x-','LineWidth',2); title ('abs($$\alpha^{ez}_{j\theta}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_hxjx),'x-','LineWidth',2); title ('abs($$\alpha^{h\rho}_{j\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_hyjy),'x-','LineWidth',2); title ('abs($$\alpha^{h\theta}_{j\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_hxjy),'x-','LineWidth',2); title ('abs($$\alpha^{h\rho}_{j\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_hyjx),'x-','LineWidth',2); title ('abs($$\alpha^{h\theta}_{j\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_curr_ezjx),'x-','LineWidth',2); title ('abs($$\alpha^{ez}_{j\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

