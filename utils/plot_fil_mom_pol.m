E_SOL_st_MoM_final = [];
E_SOL_st_filaments_mul_final = [];
E_SOL_st_pol_final = [];
for i = 1:size(PARAMETERS_PLANE,1)
    for j = 1:size(PARAMETERS_PLANE,2)
        
        params_final = PARAMETERS_PLANE{i,j};
        E_SOL_st_MoM_final(:,j) = E_SOL_st_MoM_MEAN_PLANE{i,j};
        E_SOL_st_filaments_mul_final (:,j) = E_SOL_st_filaments_mul_E_PLANE{i,j}.';
        E_SOL_st_pol_final (:,j) = full(E_SOL_st_POL_PLANE{i,j}).';
        
        
        OMEGA_final(j) = params_final.OMEGA_vec(j)/params_final.omega;
        %         rho_c(mm,d) = params.shift_vec(d);

        
        
    end
end
figure;

for ll=1:size((E_SOL_st_MoM_final),1)
    
    p1(ll) = plot(OMEGA_final,abs(E_SOL_st_MoM_final(ll,:)),'x-'); hold on;
    p2(ll) = plot(OMEGA_final,abs(E_SOL_st_filaments_mul_final(ll,:) ),'s-');
    p3(ll) = plot(OMEGA_final,abs((E_SOL_st_pol_final(ll,:) )),'d-');
    
    legend([p1(1), p2(1), p3(1)], {'MoM' , 'Fil', 'Pol'});
    
    
end

title(['abs(E) vs \Omega/\omega, \lambda =',num2str(params_final.lambda*1e6),'\mum, Radius = ',num2str(params_final.radius/params_final.lambda,4),'\lambda,', ' epsilon = ',num2str(params_final.er_in)])
xlabel('\Omega/\omega');

grid on;
grid minor;



figure;

for ll=1:size((E_SOL_st_MoM_final),1)
    
    p1(ll) = plot(OMEGA_final,real(E_SOL_st_MoM_final(ll,:)),'x-'); hold on;
    p2(ll) = plot(OMEGA_final,real(E_SOL_st_filaments_mul_final(ll,:) ),'s-');
    p3(ll) = plot(OMEGA_final,real((E_SOL_st_pol_final(ll,:) )),'d-');
    
    legend([p1(1), p2(1), p3(1)], {'MoM' , 'Fil', 'Pol'});
    
    
end

title(['real(E) vs \Omega/\omega, \lambda =',num2str(params_final.lambda*1e6),'\mum, Radius = ',num2str(params_final.radius/params_final.lambda,4),'\lambda,', ' epsilon = ',num2str(params_final.er_in)])
xlabel('\Omega/\omega');

grid on;
grid minor;



figure;

for ll=1:size((E_SOL_st_MoM_final),1)
    
    p1(ll) = plot(OMEGA_final,imag(E_SOL_st_MoM_final(ll,:)),'x-'); hold on;
    p2(ll) = plot(OMEGA_final,imag(E_SOL_st_filaments_mul_final(ll,:) ),'s-');
    p3(ll) = plot(OMEGA_final,imag((E_SOL_st_pol_final(ll,:) )),'d-');
    
    legend([p1(1), p2(1), p3(1)], {'MoM' , 'Fil', 'Pol'});
    
    
end

title(['imag(E) vs \Omega/\omega, \lambda =',num2str(params_final.lambda*1e6),'\mum, Radius = ',num2str(params_final.radius/params_final.lambda,4),'\lambda,', ' epsilon = ',num2str(params_final.er_in)])
xlabel('\Omega/\omega');

grid on;
grid minor;

figure; scatter(params_final.sca_x*1e6, params_final.sca_y*1e6,500,abs(E_SOL_st_MoM_final(:,1)),'.');
colormap('jet');
grid on; grid minor; xlabel ('\lambda'); ylabel('\lambda');
title ('Scatterers Array');
