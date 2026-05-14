function TM_alpha_presentation_fields(results)
% TM_alpha_presentation_fields  Post-process results from eval_TM_results into
% the cell-array layout this presentation logic was originally written for,
% then run the regression fits and produce the comparison plots.
%
% Input:
%   results : cell array { i, j, t } of structs returned by eval_TM_results.

[ni, nj, nt] = size(results);

field_current_mat_output   = cell(ni, nj, nt);
E_SOL_st_MoM_MEAN_PLANE_E0 = cell(ni, nj, nt);
E_SOL_st_MoM_MEAN_PLANE    = cell(ni, nj, nt);

for ii = 1:ni
    for jj = 1:nj
        for tt = 1:nt
            r = results{ii,jj,tt};
            field_current_mat_output{ii,jj,tt} = r.filaments.field_current_mat;
            if isfield(r.mom, 'mean_E0'), E_SOL_st_MoM_MEAN_PLANE_E0{ii,jj,tt} = r.mom.mean_E0; end
            if isfield(r.mom, 'mean_E'),  E_SOL_st_MoM_MEAN_PLANE{ii,jj,tt}    = r.mom.mean_E;  end
        end
    end
end

% The legacy body below expects a `params` variable in scope.
params = results{1,1,1}.params;

%fields only
stm = [];
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
        
        stm = [stm; (-params.omega^2)*params.mu0*(params.er_in-params.er_out)*params.e0*(Ez_sol-Ez)-2*(-params.omega^2)*params.mu0/(params.c^2)*params.OMEGA_vec(OMEGA_r)*params.shift_vec(rho)*(Hy_sol-Hy)]
 
        
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
%         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=inv(lineq_for_x)*([Iz_1;Iz_2;Iz_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
        
%         lineq_for_x = [0 , 0 ,Ez_1; 0 ,0, Ez_2; 0 ,0, Ez_3];
%         a=linsolve(lineq_for_x,[Iz_1;Iz_2;Iz_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
     
%         
%         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=linsolve(lineq_for_x,[Iez_1;Iez_2;Iez_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
%         
%             lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=linsolve(lineq_for_x,[Iz_1-Iez_1;Iz_2-Iez_2;Iz_3-Iez_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
        %
%         lineq_for_x = [Hx , Hy ,Ez];
        a=linsolve((lineq_for_x),transpose([Hx_sol]));
        alpha_hxhx(rho,OMEGA_r) = a(1);
        alpha_hxhy(rho,OMEGA_r) = a(2);
        alpha_hxez(rho,OMEGA_r) = a(3);
        
%         lineq_for_x = [Hx , Hy ,Ez];
        a=linsolve((lineq_for_x),transpose([Hy_sol]));
        alpha_hyhx(rho,OMEGA_r) = a(1);
        alpha_hyhy(rho,OMEGA_r) = a(2);
        alpha_hyez(rho,OMEGA_r) = a(3);
        

        
        fact_a = -1i*params.omega*(params.er_in-params.er_out)*params.e0 * params.radius^2*pi;
        fact_b = -1i*params.omega*(1/params.c^2)*params.OMEGA_vec(OMEGA_r)*params.shift_vec(rho)*params.radius^2*pi;
        
        alpha_curr_ezjz(rho,OMEGA_r) = fact_a*alpha_ezez(rho,OMEGA_r)+fact_b*alpha_hxez(rho,OMEGA_r);
        alpha_curr_hyjz(rho,OMEGA_r) = fact_a*alpha_ezhy(rho,OMEGA_r)+fact_b*alpha_hxhy(rho,OMEGA_r);
        alpha_curr_hxjz(rho,OMEGA_r) = fact_a*alpha_ezhx(rho,OMEGA_r)+fact_b*(alpha_hxhx(rho,OMEGA_r)-1);

        
        fact_am = -1i*params.omega*(1/params.c^2)*params.OMEGA_vec(OMEGA_r)*params.shift_vec(rho)* params.radius^2*pi;
        fact_bm = -1i*params.omega*params.radius^2*pi*params.mu0*(params.mr_in - params.mr_out);
        

        
        alpha_curr_ezjx(rho,OMEGA_r) = fact_am*(alpha_ezez(rho,OMEGA_r)-1)+fact_bm*alpha_hxez(rho,OMEGA_r);
        alpha_curr_hyjx(rho,OMEGA_r) = fact_am*(alpha_ezhy(rho,OMEGA_r))+fact_bm*alpha_hxhy(rho,OMEGA_r);
        alpha_curr_hxjx(rho,OMEGA_r) = fact_am*(alpha_ezhx(rho,OMEGA_r))+fact_bm*(alpha_hxhx(rho,OMEGA_r)-1);

        alpha_curr_ezjy(rho,OMEGA_r) = 0*fact_am*alpha_hyez(rho,OMEGA_r)+fact_bm*alpha_hxez(rho,OMEGA_r);
        alpha_curr_hyjy(rho,OMEGA_r) = 0*fact_am*(alpha_hyhy(rho,OMEGA_r)-1)+fact_bm*alpha_hxhy(rho,OMEGA_r);
        alpha_curr_hxjy(rho,OMEGA_r) = 0*fact_am*alpha_hyhx(rho,OMEGA_r)+fact_bm*(alpha_hxhx(rho,OMEGA_r)-1);

        
%         %% fields alpha
%         a=linsolve(lineq_for_x,[Hx_sol_1;Hx_sol_2;Hx_sol_3]);
%         alpha_hxhx(rho,OMEGA_r) = a(1);
%         alpha_hxhy(rho,OMEGA_r) = a(2);
%         alpha_hxez(rho,OMEGA_r) = a(3);
%         
%         %         lineq_for_x = [Ez_1; Ez_2;Ez_3];
%         a=linsolve(lineq_for_x,[Hy_sol_1;Hy_sol_2;Hy_sol_3]);
%         alpha_hyhx(rho,OMEGA_r) = a(1);
%         alpha_hyhy(rho,OMEGA_r) = a(2);
%         alpha_hyez(rho,OMEGA_r) = a(3);
%         
%         
%           lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=linsolve(lineq_for_x,[Ez_sol_1;Ez_sol_2;Ez_sol_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
%         
        
        
        
        
        
        
        
        
        
        
        
%         lineq_for_x = [Ez_1;  Ez_2;  Ez_3];
%         alpha_TM(rho,OMEGA_r)=linsolve(lineq_for_x,[Ez_sol_1;Ez_sol_2;Ez_sol_3]);
        
         lineq_for_x = transpose([Ez_inc_mom]);
        a=linsolve((lineq_for_x),transpose([Ez_sol_mom]));
        alpha_TM_MoM(rho,OMEGA_r) = a;
        
         lineq_for_x = transpose([Ez]);
        a=linsolve((lineq_for_x),transpose([Ez_sol]));
        alpha_TM(rho,OMEGA_r) = a;
        
        %                 lineq_for_x = [Ez_1;  Ez_2;  Ez_3];
        %         alpha_TM(rho,OMEGA_r)=linsolve(lineq_for_x,[Iz_1;Iz_2;Iz_3]);
        %
%         lineq_for_x = [Ez_inc_MoM_1;  Ez_inc_MoM_2;  Ez_inc_MoM_3];
%         alpha_TM_MoM(rho,OMEGA_r)=linsolve(lineq_for_x,[Ez_sol_MoM_1;Ez_sol_MoM_2;Ez_sol_MoM_3]);
        
        
%         %         alpha_TM_MoM(rho,OMEGA_r) = E_SOL_st_alpha_MoM{rho,OMEGA_r,1};
% 
%         %         alpha_hxez(rho,OMEGA_r) = a(1);
%         %         alpha_hyez(rho,OMEGA_r) = a(2);
%         %         alpha_ezez(rho,OMEGA_r) = a(3);
%         
        
        %          b(rho,OMEGA_r) = mean([Iez_1 / Ez_1,Iez_2 / Ez_2,Iez_3 / Ez_3]);
        %          lineq_for_x = [Hx_1 , Hy_1 ; Hx_2 ,Hy_2; Hx_3 ,Hy_3];
        %         a=linsolve(lineq_for_x,[Iz_1-Iez_1;Iz_2-Iez_2;Iz_3-Iez_3]);
        %         alpha_lineqhx(rho,OMEGA_r) = a(1);
        %         alpha_lineqhy(rho,OMEGA_r) = a(2);
        %         alpha_lineqez(rho,OMEGA_r) = b(rho,OMEGA_r);
        
        %         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
        %         a=inv(lineq_for_x)*[Iz_1;Iz_2;Iz_3];
        %         alpha_lineqhx(rho,OMEGA_r) = a(1);
        %         alpha_lineqhy(rho,OMEGA_r) = a(2);
        %         alpha_lineqez(rho,OMEGA_r) = a(3);
        %
        
        %         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
        %         a=linsolve(lineq_for_x,[Iy_1;Iy_2;Iy_3]);
        %         alpha_lineqyx(rho,OMEGA_r) = a(1);
        %         alpha_lineqyy(rho,OMEGA_r) = a(2);
        %         alpha_lineqyh(rho,OMEGA_r) = a(3);
        
        %         lineq_for_x = [Ex_1 , Ey_1; Ex_2 ,Ey_2; Ex_3 ,Ey_3];
        %         a=linsolve(lineq_for_x,[Iy_1;Iy_2;Iy_3]);
        %         alpha_lineqyx(rho,OMEGA_r) = a(1);
        %         alpha_lineqyy(rho,OMEGA_r) = a(2);
        %         alpha_lineqyh(rho,OMEGA_r) = 0;
        %
        %         lineq_for_x = [Ex_1 , Ey_1; Ex_2 ,Ey_2; Ex_3 ,Ey_3];
        %         a=linsolve(lineq_for_x,[Ix_1;Ix_2;Ix_3]);
        %         alpha_lineqxx(rho,OMEGA_r) = a(1);
        %         alpha_lineqxy(rho,OMEGA_r) = a(2);
        %         alpha_lineqxh(rho,OMEGA_r) = 0;
        
        %         lineq_for_y = [Ex_1 , Ey_1; Ex_2 ,Ey_2];
        %         a=linsolve(lineq_for_y,[Iy_1;Iy_2]);
        %         alpha_lineqyx(rho,OMEGA_r) = a(1);
        %         alpha_lineqyy(rho,OMEGA_r) = a(2);
        
        
    end
end

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
%
% % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_yy),'x-','LineWidth',2); title ('\alpha_{\theta\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_yx),'x-','LineWidth',2); title ('\alpha_{\rho\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,abs(alpha_xx),'x-','LineWidth',2); title ('\alpha_{\rho\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs(alpha_xy),'x-','LineWidth',2); title ('\alpha_{\theta\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % legend (num2str(params.shift_vec/params.lambda'),'LineWidth',2);
% %
% %
% % figure;
% % subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyy),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxy),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxx),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % title(lgd1,'\rho_c [\lambda]')
% % subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyx),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% %
% % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyy),'*-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqxy),'*-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,real(alpha_lineqxx),'*-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % title(lgd1,'\rho_c [\lambda]')
% % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyx),'*-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% %
% % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyy),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxy),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxx),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % title(lgd1,'\rho_c [\lambda]')
% % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyx),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% %
% %
% %
% %
% % figure; subplot(2,2,4); (surf(abs(alpha_lineqyy))); title ('E at $$\hat{\theta} , \alpha_{\theta\theta}$$','Interpreter','Latex'); grid on; grid minor;
% % subplot(2,2,2); (surf(abs(alpha_lineqxy))); title ('E at $$\hat{\theta} , \alpha_{\rho\theta}$$','Interpreter','Latex'); grid on; grid minor;
% % subplot(2,2,1);(surf(abs(alpha_lineqxx))); title ('E at $$\hat{\rho} , \alpha_{\rho\rho}$$','Interpreter','Latex');grid on; grid minor;
% % subplot(2,2,3); (surf(abs(alpha_lineqyx))); title ('E at $$\hat{\rho} , \alpha_{\theta\rho}$$','Interpreter','Latex');grid on; grid minor;
% %
%
% figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez),'x-','LineWidth',2); title ('abs($$\alpha_{Ez}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhy),'x-','LineWidth',2); title ('abs($$\alpha_{H\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhx),'x-','LineWidth',2); title ('abs($$\alpha_{H\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Ix),'x-','LineWidth',2); title ('abs($$\alpha_{mEz\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Iy),'x-','LineWidth',2); title ('abs($$\alpha_{mEz\theta}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez),'x-','LineWidth',2); title ('real($$\alpha_{Ez}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhy),'x-','LineWidth',2); title ('real($$\alpha_{H\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhx),'x-','LineWidth',2); title ('real($$\alpha_{H\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Ix),'x-','LineWidth',2); title ('real($$\alpha_{mEz\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Iy),'x-','LineWidth',2); title ('real($$\alpha_{mEz\theta}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez),'x-','LineWidth',2); title ('imag($$\alpha_{Ez}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhy),'x-','LineWidth',2); title ('imag($$\alpha_{H\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhx),'x-','LineWidth',2); title ('imag($$\alpha_{H\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Ix),'x-','LineWidth',2); title ('imag($$\alpha_{mEz\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Iy),'x-','LineWidth',2); title ('imag($$\alpha_{mEz\theta}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%
%
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhy),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhx),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Ix),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Iy),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhy),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhx),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Ix),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Iy),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhy),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhx),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Ix),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Iy),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot((alpha_lineqez),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,(alpha_lineqhy),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,(alpha_lineqhx),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,(alpha_lineqez_Ix),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,(alpha_lineqez_Iy),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
%  figure;
%  plot(alpha_lineqez,'o-'); grid on; grid minor;
%   figure;
%  plot(alpha_lineqhx,'.-'); grid on; grid minor;
%
%
%  figure;
% subplot(1,3,1); plot(params.OMEGA_vec/params.omega,abs(b),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
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

% 
% figure;
% subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez - alpha_ezez(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy - alpha_ezhy(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx - alpha_ezhx(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_hxez - alpha_hxez(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hyez - alpha_hyez(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhx - alpha_hxhx(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhy - alpha_hxhy(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhx - alpha_hyhx(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhy - alpha_hyhy(:,1)),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 


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

figure; plot(params.OMEGA_vec/params.omega,reshape(abs(stm(:,5)),5,8));xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
figure; plot(params.OMEGA_vec/params.omega,reshape(abs(stm(:,3)),5,8));xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

figure; plot(params.OMEGA_vec/params.omega,reshape(abs(stm(:,1)),5,8));xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;




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

end
