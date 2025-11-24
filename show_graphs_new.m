mkdir temp_results

for i = 1:length(alpha_ezez)
    for k=1:length(alpha_ezez(:,1))
   
        if params.OMEGA_vec(i)==0
            alpha_ezez(k,i) = abs(alpha_ezez(k,i))*exp(1i*(angle(alpha_ezez(k,i-1) )+angle(alpha_ezez(k,i+1)))/2) 
            alpha_ezhy(k,i) = abs(alpha_ezhy(k,i))*exp(1i*(angle(alpha_ezhy(k,i-1) )+angle(alpha_ezhy(k,i+1)))/2) 
            alpha_ezhx(k,i) = abs(alpha_ezhx(k,i))*exp(1i*(angle(alpha_ezhx(k,i-1) )+angle(alpha_ezhx(k,i+1)))/2) 
            alpha_hxez(k,i) = abs(alpha_hxez(k,i))*exp(1i*(angle(alpha_hxez(k,i-1) )+angle(alpha_hxez(k,i+1)))/2) 
            alpha_hyez(k,i) = abs(alpha_hyez(k,i))*exp(1i*(angle(alpha_hyez(k,i-1) )+angle(alpha_hyez(k,i+1)))/2) 
            alpha_hxhx(k,i) = abs(alpha_hxhx(k,i))*exp(1i*(angle(alpha_hxhx(k,i-1) )+angle(alpha_hxhx(k,i+1)))/2) 
            alpha_hxhy(k,i) = abs(alpha_hxhy(k,i))*exp(1i*(angle(alpha_hxhy(k,i-1) )+angle(alpha_hxhy(k,i+1)))/2) 
            alpha_hyhx(k,i) = abs(alpha_hyhx(k,i))*exp(1i*(angle(alpha_hyhx(k,i-1) )+angle(alpha_hyhx(k,i+1)))/2) 
            alpha_hyhy(k,i) = abs(alpha_hyhy(k,i))*exp(1i*(angle(alpha_hyhy(k,i-1) )+angle(alpha_hyhy(k,i+1)))/2) 
        end
    end
end

hfig = figure;
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

saveas (hfig, ['temp_results', '\sigma_phase.fig'])


hfig = figure;
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
saveas (hfig, ['temp_results', '\sigma_abs.fig'])

hfig = figure;
plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('$$|\sigma^{ee}_{zz}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

saveas (hfig, ['temp_results', '\sigma_zz_abs.fig'])

hfig  = figure;
plot(params.OMEGA_vec/params.omega,angle(alpha_ezez),'x-','LineWidth',2); title ('$$\angle\sigma^{ee}_{zz}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
saveas (hfig, ['temp_results', '\sigma_zz_phase.fig'])

hfig = figure;
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('$$|\sigma^{em}_{z\rho}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
saveas (hfig, ['temp_results', '\sigma_zx_abs.fig'])

hfig = figure;
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx),'x-','LineWidth',2); title ('$$\angle\sigma^{em}_{z\rho}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
saveas (hfig, ['temp_results', '\sigma_zx_phase.fig'])

hfig = figure;
plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('$$|\sigma^{em}_{z\theta}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
saveas (hfig, ['temp_results', '\sigma_zy_abs.fig'])

hfig = figure;
plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy),'x-','LineWidth',2); title ('$$\angle\sigma^{em}_{z\theta}$$','Interpreter','Latex', 'FontSize', 16); xlabel ('$$\bar{\Omega}$$','Interpreter','Latex', 'FontSize', 16); grid on; grid minor;
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
saveas (hfig, ['temp_results', '\sigma_zy_abs.fig'])



hfig = figure;
subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('Fil abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('Fil abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('Fil abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

        
subplot(2,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx_mom),'x-','LineWidth',2); title ('MoM abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
       
lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')
saveas (hfig, ['temp_results', '\sigma_e_abs_mom_fil.fig'])



hfig = figure;
subplot(2,3,3); plot(params.OMEGA_vec/params.omega,angle(alpha_ezez),'x-','LineWidth',2); title ('Fil angle($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,2); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy),'x-','LineWidth',2); title ('Fil angle($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,1); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx),'x-','LineWidth',2); title ('Fil angle($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

        
subplot(2,3,6); plot(params.OMEGA_vec/params.omega,angle(alpha_ezez_mom),'x-','LineWidth',2); title ('MoM angle($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,5); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhy_mom),'x-','LineWidth',2); title ('MoM angle($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(2,3,4); plot(params.OMEGA_vec/params.omega,angle(alpha_ezhx_mom),'x-','LineWidth',2); title ('MoM angle($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
       
lgd1 = legend (num2str(params.shift_vec'/params.lambda));
title(lgd1,'\rho_c [\lambda]')
saveas (hfig, ['temp_results', '\sigma_e_phase_mom_fil.fig'])

