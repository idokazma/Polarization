len_n = params.len_n;
wid_n = params.wid_n;

len = params.len;
wid = params.wid;

y = params.y;
x = params.x;


X = params.X;
Y = params.Y;


[X_mom,Y_mom] = meshgrid(params.x_mom,params.y_mom);



delta_eps = params.er_in - params.er_out;
% scattereres (loc) = delta_eps+epsilon_bg;

% loc =find(sqrt((X_mom).^2 + (Y_mom).^2)<params.radius);

% loc =find(abs(X_mom*-0.25+Y_mom)<0.003*params.lambda & abs(X_mom*1+Y_mom)<0.005*params.lambda);
loc =find(abs(X_mom)<0.005*params.lambda);

rot_angle = 0*pi/6;
X_mom_new = real((X_mom+1i*(Y_mom)).*exp(1i*rot_angle));
Y_mom_new = imag((X_mom+1i*(Y_mom)).*exp(1i*rot_angle));
X_mom=X_mom_new;
Y_mom=Y_mom_new;
total_cells = len_n*wid_n;
epsilon_bg=params.er_out;
lambda = params.lambda;
k0 = 2*pi/lambda;

beta = 1;

V_i = abs(params.x_mom(2)-params.x_mom(1))*abs(params.y_mom(2)-params.y_mom(1));
norm_fact = (pi*params.radius^2)/(length(loc)*V_i);
effective_radius = sqrt(V_i*length(loc)/pi);


R_i_eff = sqrt(V_i/pi); %% R_i_eff = sqrt(V_i/pi);
gamma=R_i_eff/k0*besselh(1,1,R_i_eff*k0)+1i*2/(pi*k0^2);
V_i_bg =abs(params.x(2)-params.x(1))*abs(params.y(2)-params.y(1));
Mi = (1+1i*(2/pi*log(k0*params.n_out*sqrt(V_i)/(4*exp(1)))))*V_i;



sca_area = (effective_radius^2*pi);

sigma_res = zeros(length(params.shift_vec),length(params.OMEGA_vec));
sigma_res_sanity =zeros(length(params.shift_vec),length(params.OMEGA_vec));
sigma_res_hy = zeros(length(params.shift_vec),length(params.OMEGA_vec));
sigma_res_hx = zeros(length(params.shift_vec),length(params.OMEGA_vec));
sigma_res_de = zeros(length(params.shift_vec),length(params.OMEGA_vec));
H_t = -1i*params.omega*params.mr_out*params.mu0;

for shift_ind = 1:length(params.shift_vec)
            params.sca_x = 1*params.shift_vec(shift_ind);
            params.sca_y = 0*params.shift_vec(shift_ind);
for k = 1:length(params.OMEGA_vec)
    params.OMEGA = params.OMEGA_vec(k);
    disp(k/length(params.OMEGA_vec)*100)
    integ = 0;
    integ_2 = 0;
    integ_3 = 0;
    integ_hy = 0;
    integ_hx = 0;
    integ_de = 0;
    for j = 1:length(loc)
        x_temp = X_mom(loc(j))+params.sca_x;
        y_temp = Y_mom(loc(j))+params.sca_y;
        integ_2_temp = 0;
        H_t = -1i*params.omega*params.mr_out*params.mu0;
        for i=1:length(loc)
            xp_temp = X_mom(loc(i))+params.sca_x;
            yp_temp = Y_mom(loc(i))+params.sca_y;
            if i==j
                g = k0^2*delta_eps*Mi;
            else
                g = k0^2*delta_eps*V_i*scalar_green_mom([xp_temp;yp_temp],[x_temp;y_temp],params.n_out,k0,c,1*params.OMEGA);
            end
            integ = integ+V_i*g;
            integ_2_temp = integ_2_temp+g;
            
            
        end
        integ_3 =integ_3+V_i*(1/(1-integ_2_temp));
        integ_hy =integ_hy+V_i*(H_t*(x_temp-params.sca_x(1))/(1-integ_2_temp));
        integ_hx =integ_hx+V_i*(H_t*(y_temp)/(1-integ_2_temp));
        integ_de = integ_de+V_i*(+1i*params.OMEGA/params.omega*params.k0^2*y_temp*params.sca_x(1))/(1-integ_2_temp);
        %         integ_2_temp = integ_2_temp + 
    end
    sigma_res(shift_ind,k) = 1/(1-integ/sca_area);
    sigma_res_sanity(shift_ind,k) = integ_3/sca_area;
    sigma_res_hy(shift_ind,k) = integ_hy/sca_area;
    sigma_res_hx(shift_ind,k) = integ_hx/sca_area;
    sigma_res_de(shift_ind,k) = integ_de/sca_area;
    integ = 0;
end
end
figure;
% subplot(1,5,1);
% plot(abs(sigma_res.')); hold on;
subplot(1,3,1);

plot(abs(sigma_res_sanity.')); hold on;
subplot(1,3,2);

plot(abs(sigma_res_hy.')); hold on;
subplot(1,3,3);

plot(abs(sigma_res_hx.')); hold on;


if params.n_in>params.n_out
    figure('Name', 'ctm4_Abs_pos');

    subplot(1,3,1); plot(params.OMEGA_vec/params.omega,abs(sigma_res_sanity) ,'LineWidth',2); title ('$$|\sigma^{ee}_{zz}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
    subplot(1,3,2); plot(params.OMEGA_vec/params.omega,abs(sigma_res_hy) ,'LineWidth',2); title ('$$|\sigma^{em}_{z\theta}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
    subplot(1,3,3); plot(params.OMEGA_vec/params.omega,abs(sigma_res_hx) ,'LineWidth',2); title ('$$|\sigma^{em}_{z\rho}|$$','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

    lgd1 = legend (num2str(params.shift_vec'/params.lambda));

    title(lgd1,'\rho_c [\lambda]')

else
        figure('Name', 'ctm4_Abs_neg');

    subplot(1,3,1); plot(params.OMEGA_vec/params.omega,abs(sigma_res_sanity) ,'LineWidth',2); title ('$$|\sigma^{ee}_{zz}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
    subplot(1,3,2); plot(params.OMEGA_vec/params.omega,abs(sigma_res_hy) ,'LineWidth',2); title ('$$|\sigma^{em}_{z\theta}|$$','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
    subplot(1,3,3); plot(params.OMEGA_vec/params.omega,abs(sigma_res_hx) ,'LineWidth',2); title ('$$|\sigma^{em}_{z\rho}|$$','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

    lgd1 = legend (num2str(params.shift_vec'/params.lambda));

    title(lgd1,'\rho_c [\lambda]')
end
% subplot(1,5,5);
% 
% plot(abs(sigma_res_de)); hold on;
% 
% figure;
% plot(abs(sigma_res_de+sigma_res_sanity)); hold on;