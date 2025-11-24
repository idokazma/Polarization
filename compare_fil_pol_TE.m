for rho = 1 : size(alpha_te_fil_output,1)
    for OMEGA_r = 1 : size(alpha_te_fil_output,2)
%         temp_alpha = alpha_te_fil_output{rho,OMEGA_r,1};
%         alpha_yy(rho,OMEGA_r) = temp_alpha(2,2);
%         alpha_yx(rho,OMEGA_r) = temp_alpha(1,2);
        field_current{rho,OMEGA_r,1} = field_current_mat_output{rho,OMEGA_r,1};
        pol_current{rho,OMEGA_r,1} = Pvec_output{rho,OMEGA_r,1};
        
        fil_x(OMEGA_r,:) = field_current{rho,OMEGA_r,1}{6};
        fil_y(OMEGA_r,:) = field_current{rho,OMEGA_r,1}{7};

        hit_x(OMEGA_r,:) = field_current{rho,OMEGA_r,1}{3};
        hit_y(OMEGA_r,:) = field_current{rho,OMEGA_r,1}{4};
        
        pol_x(OMEGA_r,:) = pol_current{rho,OMEGA_r,1}(1:end/2);
        pol_y(OMEGA_r,:) = pol_current{rho,OMEGA_r,1}(end/2+1:end);
 
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


figure; p1= plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,abs(hit_x),'x-'); hold on;  p2=plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,abs(hit_y),'d-'); grid on; grid minor; 
title(['Ex and Ey vs \Omega/\omega, \lambda =',num2str(params.lambda*1e6),'\mum, Radius = ',num2str(params.radius/params.lambda,4),'\lambda,', ' epsilon = ',num2str(params.er_in)])
xlabel('\Omega/\omega');
ylabel('|E^{inc}| [V/m]')
legend([p1(1),p2(1)], {'Ex', 'Ey'})

figure; p1 = plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,abs(fil_x),'x-'); hold on;  p2=plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,abs(pol_x),'d-'); grid on; grid minor; xlabel('\Omega/\omega'); ylabel('Amp');
title(['|P_{x}| vs \Omega/\omega, \lambda =',num2str(params.lambda*1e6),'\mum, Radius = ',num2str(params.radius/params.lambda,4),'\lambda,', ' epsilon = ',num2str(params.er_in)])
xlabel('\Omega/\omega');
ylabel('Amp')
legend([p1(1),p2(1)], {'fil', 'pol'})


figure; p1=plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,abs(fil_y),'x-'); hold on;  p2=plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,abs(pol_y),'d-'); grid on; grid minor;  title('ABS Py'); xlabel('\Omega/\omega');ylabel('Amp');
title(['|P_{y}| vs \Omega/\omega, \lambda =',num2str(params.lambda*1e6),'\mum, Radius = ',num2str(params.radius/params.lambda,4),'\lambda,', ' epsilon = ',num2str(params.er_in)])
xlabel('\Omega/\omega');
ylabel('Amp')
legend([p1(1),p2(1)], {'fil', 'pol'})

figure; p1 = plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,sqrt(abs(fil_x).^2 + abs(fil_y).^2),'x-'); hold on;  p2=plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,sqrt(abs(pol_x).^2 + abs(pol_y).^2),'d-');grid on; grid minor;  title('ABS P');xlabel('\Omega/\omega');ylabel('Amp');
title(['|P| vs \Omega/\omega, \lambda =',num2str(params.lambda*1e6),'\mum, Radius = ',num2str(params.radius/params.lambda,4),'\lambda,', ' epsilon = ',num2str(params.er_in)])
xlabel('\Omega/\omega');
ylabel('Amp')
legend([p1(1),p2(1)], {'fil', 'pol'})


figure; plot(params.OMEGA_vec(1:size(hit_x,1))/params.omega,sqrt(abs(hit_x).^2 + abs(hit_y).^2),'x-');

% figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_yy),'x-','LineWidth',2); title ('\alpha_{\theta\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_yx),'x-','LineWidth',2); title ('\alpha_{\rho\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1);plot(params.OMEGA_vec/params.omega,abs(alpha_xx),'x-','LineWidth',2); title ('\alpha_{\rho\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs(alpha_xy),'x-','LineWidth',2); title ('\alpha_{\theta\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% legend (num2str(params.shift_vec/params.lambda'),'LineWidth',2);
% 
% 
% figure;
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyy),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxy),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxx),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyx),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% 
% figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyy),'*-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqxy),'*-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1);plot(params.OMEGA_vec/params.omega,real(alpha_lineqxx),'*-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyx),'*-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% 
% figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyy),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxy),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1);plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxx),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyx),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% 
% 
% 
% 
% figure; subplot(2,2,4); (surf(abs(alpha_lineqyy))); title ('E at $$\hat{\theta} , \alpha_{\theta\theta}$$','Interpreter','Latex'); grid on; grid minor;
% subplot(2,2,2); (surf(abs(alpha_lineqxy))); title ('E at $$\hat{\theta} , \alpha_{\rho\theta}$$','Interpreter','Latex'); grid on; grid minor;
% subplot(2,2,1);(surf(abs(alpha_lineqxx))); title ('E at $$\hat{\rho} , \alpha_{\rho\rho}$$','Interpreter','Latex');grid on; grid minor;
% subplot(2,2,3); (surf(abs(alpha_lineqyx))); title ('E at $$\hat{\rho} , \alpha_{\theta\rho}$$','Interpreter','Latex');grid on; grid minor;
% 

% % % figure;
% % % subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_eyey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_eyex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_exex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % %  lgd1 = legend (num2str(params.shift_vec'/params.lambda));
% % %  title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,3,4); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_exey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_hzex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hzey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % figure;
% % % subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_eyey),'x-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_eyex),'x-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_exex),'x-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % %  lgd1 = legend (num2str(params.shift_vec'/params.lambda));
% % %  title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,3,4); p4 =plot(params.OMEGA_vec/params.omega,real(alpha_exey),'x-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_hzex),'x-','LineWidth',2); title ('real($$\alpha_{\rho H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hzey),'x-','LineWidth',2); title ('real($$\alpha_{\theta H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % 
% % % figure;
% % % subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_eyey),'x-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_eyex),'x-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_exex),'x-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % %  lgd1 = legend (num2str(params.shift_vec'/params.lambda));
% % %  title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,3,4); p4 =plot(params.OMEGA_vec/params.omega,imag(alpha_exey),'x-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_hzex),'x-','LineWidth',2); title ('imag($$\alpha_{\rho H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,6); plot(params.OMEGA_vec/params.omega,imag(alpha_hzey),'x-','LineWidth',2); title ('imag($$\alpha_{\theta H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % 
% % % 
% % % figure;
% % % subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_eyey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_eyex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,1); plot(params.OMEGA_vec/params.omega,abs(alpha_exex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % % % title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_exey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,real(alpha_eyey),'*-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,real(alpha_eyex),'*-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega'); grid on; grid minor;
% % % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,real(alpha_exex),'*-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % % % title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,real(alpha_exey),'*-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,imag(alpha_eyey),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,imag(alpha_eyex),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,imag(alpha_exex),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % % % title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag(alpha_exey),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

% 
% figure;
% subplot(2,2,1); plot(params.OMEGA_vec/params.omega,real((alpha_lineqxx + alpha_lineqyx)./alpha_lineqxx(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,atan2(real(alpha_lineqxx),imag(alpha_lineqyx)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag((alpha_lineqxx + alpha_lineqyx)./alpha_lineqxx(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,atan2(imag(alpha_lineqxx),imag(alpha_lineqyx)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% 
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs((alpha_lineqxx + alpha_lineqyx)./alpha_lineqxx(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,angle(alpha_lineqxx + alpha_lineqyx)-angle(alpha_lineqxx(3,6)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% 
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs((alpha_lineqyy + alpha_lineqxy)./alpha_lineqyy(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,angle(alpha_lineqyy + alpha_lineqxy)-angle(alpha_lineqyy(3,6)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyx),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% 
% 
% 

% 
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_hzhz),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_eyhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_exhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hzex),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hzey),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%  
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  
%  title(lgd1,'\rho_c [\lambda]')
%  
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_hzhz),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_eyhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_exhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_hzex),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_hzey),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%  
% 
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%  
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_hzhz),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_eyhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_exhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_hzex),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_hzey),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%  
%  
 
 
 
 
 figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_hzhz),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_eyhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_exhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_hzex),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hzey),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_exex),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_exey),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_eyex),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_eyey),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));
 
 title(lgd1,'\rho_c [\lambda]')
  
 figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,real(alpha_hzhz),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,real(alpha_eyhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,real(alpha_exhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_hzex),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hzey),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_exex),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_exey),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_eyex),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_eyey),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));
 
 title(lgd1,'\rho_c [\lambda]')
 
 figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,imag(alpha_hzhz),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,imag(alpha_eyhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,imag(alpha_exhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_hzex),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,imag(alpha_hzey),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_exex),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_exey),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_eyex),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_eyey),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));
 
 title(lgd1,'\rho_c [\lambda]')
