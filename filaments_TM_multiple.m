%% set up the sources location and test points locations. test points are on the surface of the scatterer;

function [full_fields, mean_E_fil,mean_I_fil,alpha_Fil,output_mat] = filaments_TM_multiple(params,Esym)
filaments_total = [];
test_point_total = [];
for k = 1 : length(params.sca_x)
    R = params.radius;
    theta = linspace(0,2*pi,params.N_filaments+1);
    theta = theta(1:end-1);
    out_sources_loc = (params.R_out*R)*exp(1i*theta);
    in_sources_loc =(params.R_in*R)*exp(1i*theta);
    
    theta_test = linspace(0,2*pi,params.N_testpoints+1);
    theta_test = theta_test(1:end-1);
    test_loc = (1*R)*exp(1i*theta_test);
    
    %Transform to X,Y
    
    out_sources_loc = [real(out_sources_loc)+params.sca_x(k);imag(out_sources_loc)+params.sca_y(k)];
    in_sources_loc = [real(in_sources_loc)+params.sca_x(k);imag(in_sources_loc)+params.sca_y(k)];
    
    filaments_out = [out_sources_loc;ones(1,length(out_sources_loc))*k;zeros(1,length(out_sources_loc))];
    filaments_in = [in_sources_loc;ones(1,length(in_sources_loc))*k;ones(1,length(in_sources_loc))];
    filaments_total = [filaments_total , filaments_in , filaments_out];
    
    test_loc = [real(test_loc);imag(test_loc)];
    
    %calc the direction vector of the surface border (nhat)
    for i = 1:length(test_loc)
        direction(:,i) = test_loc(:,i)/norm(test_loc(:,i));
        nhat(i,:) = [direction(1,i),direction(2,i),0];
    end
    test_loc = test_loc+[params.sca_x(k);params.sca_y(k)];
    
    test_point = [test_loc ; ones(1,length(test_loc))*k; nhat'];
    
    test_point_total = [test_point_total, test_point];
    
    
end
filaments_total = transpose(filaments_total);
test_point_total = transpose(test_point_total);

%% set some EM parameters of the area
e0=params.e0;
mu0 = params.mu0;

mr_in = params.mr_in;
mr_out = params.mr_out;
er_in = params.er_in;
er_out = params.er_out;
lambda = params.lambda;


c = params.c;
f = params.f;
k0 = params.k0;
omega = params.omega;
OMEGA = params.OMEGA;

n_in = params.n_in;
n_out = params.n_out;

%% START THE ELECTROMAGNETIC CALCULATIONS

% set the Ez_incident and the Ht_incident

% now, for each Test point we are going to calculate the "impact" of a unit
% current located outside and inside the scatterer ON the surface (Test
% Points). as in Yehuda's article, we calculate first ZeI, then ZeII.
% ZeI (or ZeS) is the impact of the INSIDE currents, radiating with the
% OUTSIDE material parameters and Green function measured on the surface of
% the scatterer. ZeII is the OUTSIDE currents radiating towards the INSIDE
% with the INSIDE parameters.
syms xsym ysym zsym k0sym
% % % % Xsym = [xsym ysym zsym];
% % % %
% % % % if (~params.is_plane_wave)
% % % %     ez = 1i/4 * besselh (0 , 1 , k0sym *sqrt((xsym - params.source_loc_x).^2+(ysym- params.source_loc_y).^2));
% % % %     rot =  exp(1i*k0*OMEGA/c*(params.source_loc_y.*xsym-params.source_loc_x.*ysym));
% % % %     Esym = [0, 0, ez*rot];
% % % % end
% % % % Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym) - transpose(Xsym.*Esym(3)*params.OMEGA/(c^2*(mu0*mr_out))) ;


if (false)
    for i = 1:length(test_point_total)
        E_inc(i,:) = double(subs(params.Esym,[xsym,ysym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),k0]));
        H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),0,k0]));
        
        %     E_inc_green = scalar_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0,params.c,params.OMEGA);
        %     [H_inc_green_x , H_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0, c,params.OMEGA);
        %
    end
    disp ('End Eval Symbolic');
    
end
% % % % %
% % % % % H_inc_x = H_inc(:,1);
% % % % % H_inc_y = H_inc(:,2);
% % % % % H_inc_z = H_inc(:,3);
% % % % %
% % % % % E_inc_x = E_inc(:,1);
% % % % % E_inc_y = E_inc(:,2);
% % % % % E_inc_z = E_inc(:,3);

if (params.is_plane_wave == 1)
    
    
    for i = 1:length(test_point_total)
        E_inc(i,:) = double(subs(params.Esym,[xsym,ysym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),k0]));
        H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),0,k0]));
    end
    
    H_inc_x = H_inc(:,1);
    H_inc_y = H_inc(:,2);
    H_inc_z = H_inc(:,3);
    
    E_inc_x = E_inc(:,1);
    E_inc_y = E_inc(:,2);
    E_inc_z = E_inc(:,3);
    
else
    
    
    if (params.hit_plane == 0)
        
        E_inc_green = scalar_green([params.source_loc_x;params.source_loc_y],test_point_total',params, 0);
        [H_inc_green_x , H_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],test_point_total',params, 0);
        H_inc_green_x = H_inc_green_x;
        H_inc_green_y = H_inc_green_y;
        
        E_inc_z = E_inc_green*1 + 0*mean(E_inc_green);
        E_inc_x = 0*E_inc_z;
        E_inc_y = 0*E_inc_z;
        
        H_inc_x = H_inc_green_x*1 + 0*mean(H_inc_green_x);
        H_inc_y = H_inc_green_y*1 + 0*mean(H_inc_green_y);
        H_inc_z = 0*H_inc_green_x;
        
    else
        E_inc_green = inc_Ez_field([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0,params.c,params.OMEGA, params.direction);
        [H_inc_green_x , H_inc_green_y] = inc_Ht_field([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0, c,params.OMEGA, params.direction);
        H_inc_green_x = H_inc_green_x*1./(1i*params.omega*(params.mu0*params.mr_out));
        H_inc_green_y = H_inc_green_y*1./(1i*params.omega*(params.mu0*params.mr_out));
        
        E_inc_z = E_inc_green*1 + 0*mean(E_inc_green);
        E_inc_x = 0*E_inc_z;
        E_inc_y = 0*E_inc_z;
        
        H_inc_x = H_inc_green_x*1 + 0*mean(H_inc_green_x);
        H_inc_y = H_inc_green_y*1 + 0*mean(H_inc_green_y);
        H_inc_z = 0*H_inc_green_x;
        
    end
    
end


%%% just for test
%     E_inc_z = 1+0*E_inc_green;
%     E_inc_x = 0*E_inc_z;
%     E_inc_y = 0*E_inc_z;
%
%     H_inc_x = 0*H_inc_green_x;
%     H_inc_y = 1i/(params.omega*(params.mu0*params.mr_out))*(-1i*params.k0)+0*H_inc_green_y;
%     H_inc_z = 0*H_inc_green_x;
%% just for test2
%     E_inc_z = 0*E_inc_green + scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA);
%     E_inc_x = 0*E_inc_z;
%     E_inc_y = 0*E_inc_z;
%         [H_inc_green_x , H_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0, c,params.OMEGA);
%
%     H_inc_x = 0*E_inc_z + 1./(1i*params.omega*(params.mu0*params.mr_out)) * H_inc_green_x;
%     H_inc_y = 0*E_inc_z + 1./(1i*params.omega*(params.mu0*params.mr_out)) * H_inc_green_y;
%     H_inc_z = 0*H_inc_green_x;



for i=1:length(test_point_total)
    nhat = test_point_total(i,4:6);
    test_point_xy = [test_point_total(i,1:2)];
    
    for j=1:length(filaments_total)
        if (test_point_total(i,3) == filaments_total(j,3) && filaments_total(j,4) == 0) %% we are talking on the same scatterer and test point, and the filaments is outside the sca
            Mx(i,j) = -scalar_green(filaments_total(j,1:2)',test_point_xy',params,1);
            [Gx,Gy] = dyiadic_green(filaments_total(j,1:2)',test_point_xy',params,1);
            sol = cross(nhat,[Gx,Gy,0]);
            MHz(i,j) = -sol(3);%MINUS SIGN
            
        elseif (test_point_total(i,3) == filaments_total(j,3) && filaments_total(j,4) == 1)  %% we are talking on the same scatterer and test point, and the filaments is inside the sca
            Mx(i,j) = scalar_green(filaments_total(j,1:2)',test_point_xy',params, 0);
            [Gx,Gy] = dyiadic_green(filaments_total(j,1:2)',test_point_xy',params,0);
            sol = cross(nhat,[Gx,Gy,0]);
            MHz(i,j) = sol(3);
            
        elseif (test_point_total(i,3) ~= filaments_total(j,3) && filaments_total(j,4) == 1)  %% we are talking on other scatterer and test point, and the filaments is inside the sca
            Mx(i,j) = scalar_green(filaments_total(j,1:2)',test_point_xy',params, 0);
            [Gx,Gy] = dyiadic_green(filaments_total(j,1:2)',test_point_xy',params,0);
            sol = cross(nhat,[Gx,Gy,0]);
            MHz(i,j) = sol(3);
            
        else
            Mx(i,j) = 0;
            MHz(i,j) = 0;
            
            
        end
        
    end
    Vex(i) = -E_inc_z(i);
    sol = cross(nhat,[H_inc_x(i),H_inc_y(i),0]);
    Vhz(i) = -1*sol(3);
end

M_total = [Mx;MHz];
% figure; imagesc(real(M_total));
% P_inverse_M = ((M_total')*M_total)\(M_total');
% P_inverse_M = (transpose(M_total)*M_total)\(transpose(M_total));

% I_solution = P_inverse_M*(transpose([Vex,Vhz]));
% I_solution = (M_total)\(transpose([Vex,Vhz]));
I_solution = linsolve(M_total,transpose([Vex,Vhz]));
%  I_solution =P_inverse_M*transpose([Vex,Vhz]);




% figure; quiver(test_loc(1,:),test_loc(2,:),direction(1,:),direction(2,:));
% figure; scatter(test_loc(1,:),test_loc(2,:),20,abs(E_inc_z),'filled');
% hold on;
% scatter(in_sources_loc(1,:),in_sources_loc(2,:),20,abs(I_solution(1:length(in_sources_loc))),'filled');
% scatter(out_sources_loc(1,:),out_sources_loc(2,:),20,abs(I_solution((1+length(in_sources_loc)):2*length(in_sources_loc))),'filled');
X = params.X;
Y = params.Y;
loc_line = [X(:),Y(:)]';
E_sol = zeros(1,length(loc_line));

grid_inside_each_sca = ((X).^2 +(Y).^2)<=0.999999*R*R;
E_sol_final = zeros(length(params.sca_x),length(X(grid_inside_each_sca)));
Hx_sol_final = zeros(length(params.sca_x),length(X(grid_inside_each_sca)));
Hy_sol_final = zeros(length(params.sca_x),length(X(grid_inside_each_sca)));



for tt = 1:length(filaments_total)
    if (filaments_total(tt,4) == 0)
        in_loc_line = [X(grid_inside_each_sca)+params.sca_x(filaments_total(tt,3)),Y(grid_inside_each_sca)+params.sca_y(filaments_total(tt,3))]';
        E_sol_final(filaments_total(tt,3),:) = E_sol_final(filaments_total(tt,3),:) + I_solution(tt)*scalar_green(filaments_total(tt,1:2)',in_loc_line,params,1);
        [Hx,Hy] = dyiadic_green(filaments_total(tt,1:2)',in_loc_line,params,1);
        Hx_sol_final(filaments_total(tt,3),:) = Hx_sol_final(filaments_total(tt,3),:) + I_solution(tt)* Hx;
        Hy_sol_final(filaments_total(tt,3),:) = Hy_sol_final(filaments_total(tt,3),:) + I_solution(tt)* Hy;
        
        
    end
end



mean_E_fil = mean(E_sol_final,2);


mean_I_fil = mean_E_fil*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R);
%   alpha_Fil = mean_I_fil/scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA);
% alpha_Fil = mean_I_fil/1;

if ( params.hit_plane ==0)
hit_field_inside = scalar_green([params.source_loc_x;params.source_loc_y],in_loc_line,params,0);
%    alpha_Fil = mean((E_sol_final*(-1i*params.omega*params.e0*(params.er_in-params.er_out))))./mean(hit_field_inside)*pi*R*R;
else
hit_field_inside = inc_Ez_field([params.source_loc_x;params.source_loc_y],in_loc_line,n_out,k0,c,OMEGA, params.direction);
end
% alpha_Fil = mean_I_fil/mean(E_inc_green);
alpha_Fil = mean(mean_I_fil)./mean(hit_field_inside);



for mm = 1 : length(params.sca_x)
    
    in_loc_line = [X(grid_inside_each_sca)+params.sca_x(mm),Y(grid_inside_each_sca)+params.sca_y(mm)]';
    
    for i = 1:length(in_loc_line)
        %          E_inc_in_sca(i,:) = double(subs(params.Esym,[xsym,ysym,zsym,k0sym] ,[in_loc_line(1,i),in_loc_line(2,i),0,k0]));
        %             H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),0,k0]));
        if (params.hit_plane == 0)
            E_inc_green_inside(i,:) = scalar_green([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params,0);
            
            [H_inc_green_x_temp , H_inc_green_y_temp] = dyiadic_green([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params, 0);
            H_inc_green_x_inside(i,:) = H_inc_green_x_temp;
            H_inc_green_y_inside(i,:) = H_inc_green_y_temp;
        else
            E_inc_green_inside(i,:) = inc_Ez_field([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params.n_out,params.k0,params.c,params.OMEGA, params.direction);
            
            [H_inc_green_x_temp , H_inc_green_y_temp] = inc_Ht_field([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params.n_out,params.k0, c,params.OMEGA, params.direction);
            H_inc_green_x_inside(i,:) = 1./(1i*params.omega*(params.mu0*params.mr_out))*H_inc_green_x_temp;
            H_inc_green_y_inside(i,:) = 1./(1i*params.omega*(params.mu0*params.mr_out))*H_inc_green_y_temp;
        end
        
    end
    %      E_inc_x_inside_sca(mm,:) = E_inc_in_sca(:,1);
    %      E_inc_y_inside_sca(mm,:) = E_inc_in_sca(:,2);
    E_inc_z_inside_sca(mm,:) = E_inc_green_inside(:);
    H_inc_x_inside_sca(mm,:) = H_inc_green_x_inside(:);
    H_inc_y_inside_sca(mm,:) = H_inc_green_y_inside(:);
end

% E_inc_x_inside_sca = 0;
% E_inc_y_inside_sca = E_inc_y;


mean_Hx_fil = mean(Hx_sol_final,2);
mean_Hy_fil = mean(Hy_sol_final,2);

mean_Hx_inc = mean(H_inc_x_inside_sca,2);
mean_Hy_inc = mean(H_inc_y_inside_sca,2);

% if abs(mean_Hx_inc) < 1e-20
%     mean_Hx_inc = 0;
% end
% if abs(mean_Hy_inc) < 1e-20
%     mean_Hy_inc = 0;
% end
mean_Ez_inc = mean(E_inc_z_inside_sca,2);

% % % mean_Ez_inc = scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA);
% % % [mean_Hx_inc , mean_Hy_inc] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0, c,params.OMEGA);
% % % mean_Hx_inc = 0*mean_Ez_inc + 1./(1i*params.omega*(params.mu0*params.mr_out)) * mean_Hx_inc;
% % % mean_Hy_inc = 0*mean_Ez_inc + 1./(1i*params.omega*(params.mu0*params.mr_out)) * mean_Hy_inc;

%     E_inc_z = 0*E_inc_green + scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA);
%     E_inc_x = 0*E_inc_z;
%     E_inc_y = 0*E_inc_z;
%         [H_inc_green_x , H_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0, c,params.OMEGA);
%
%     H_inc_x = 0*E_inc_z + 1./(1i*params.omega*(params.mu0*params.mr_out)) * H_inc_green_x;
%     H_inc_y = 0*E_inc_z + 1./(1i*params.omega*(params.mu0*params.mr_out)) * H_inc_green_y;
%     H_inc_z = 0*H_inc_green_x;

for tt = 1 : length(Hy_sol_final)
    
    final_current_H(tt) = -Hx_sol_final(tt)*in_loc_line(1,tt) -  Hy_sol_final(tt)*in_loc_line(2,tt);
    rotation_current_H(tt) =  -H_inc_x_inside_sca(tt)*in_loc_line(1,tt) -  H_inc_y_inside_sca(tt)*in_loc_line(2,tt);
    Jm_fil_x(tt) = (E_sol_final(tt) - E_inc_z_inside_sca(tt))*in_loc_line(1,tt);
    Jm_fil_y(tt) = (E_sol_final(tt) - E_inc_z_inside_sca(tt))*in_loc_line(2,tt);
    
%     Jm_fil_x(tt) = (E_sol_final(tt))*in_loc_line(1,tt);
%     Jm_fil_y(tt) = (E_sol_final(tt))*in_loc_line(2,tt);
end
%
% final_current_H  = -mean_Hx_fil * params.sca_x - mean_Hy_fil * params.sca_y;
% rotation_current_H  = -mean_Hx_inc * params.sca_x - mean_Hy_inc * params.sca_y;

% H_current_add = 1/(params.c^2)*params.OMEGA*(final_current_H-rotation_current_H);
H_current_add = (final_current_H-rotation_current_H);
%      H_current_add = 1*(final_current_H);



% fprintf('\n ----------- H_diff is:');
% sum(H_current_add)
mean_Im_fil_x = -1i*params.omega*params.OMEGA/(params.c^2)* params.sca_x * (mean(E_sol_final) - mean_Ez_inc) *pi*R*R;
mean_Im_fil_y = -1i*params.omega*params.OMEGA/(params.c^2)* params.sca_y * (mean(E_sol_final)  - mean_Ez_inc) *pi*R*R;

% mean_Im_fil_x = -1i*params.omega*params.OMEGA/(params.c^2) * mean(Jm_fil_x) *pi*R*R;
% mean_Im_fil_y = -1i*params.omega*params.OMEGA/(params.c^2) * mean(Jm_fil_y) *pi*R*R;
% 


mean_I_fil = mean(E_sol_final)*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R) + 1i*params.omega/(params.c^2)*params.OMEGA*mean(H_current_add)*pi*R*R;
mean_I_fil_E_only = mean(E_sol_final)*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R);
%  mean_I_fil = mean_I_fil-mean_I_fil_E_only;
output_mat = {mean_I_fil,mean_I_fil_E_only,mean_Hx_inc,mean_Hy_inc,mean_Ez_inc , mean_Im_fil_x , mean_Im_fil_y,mean_E_fil,mean_Hx_fil,mean_Hy_fil};


%         if (params.hit_plane == 0)
% 
% 
%              E_inc_theo = scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x, params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA);
% %
%         [Hx_inc_theo , Hy_inc_theo] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.sca_x, params.sca_y]',params.n_out,params.k0, c,params.OMEGA);
%           Hx_inc_theo = 1./(1i*params.omega*(params.mu0*params.mr_out))*Hx_inc_theo;
%           Hy_inc_theo = 1./(1i*params.omega*(params.mu0*params.mr_out))*Hy_inc_theo;
% 
%            output_mat = {mean_I_fil,mean_I_fil_E_only,Hx_inc_theo,Hy_inc_theo,E_inc_theo , mean_Im_fil_x , mean_Im_fil_y,mean_E_fil,mean_Hx_fil,mean_Hy_fil};
% 
%         else
%               E_inc_theo = inc_Ez_field([params.source_loc_x;params.source_loc_y],[params.sca_x, params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA, params.direction);
% %
%         [Hx_inc_theo , Hy_inc_theo] = inc_Ht_field([params.source_loc_x;params.source_loc_y],[params.sca_x, params.sca_y]',params.n_out,params.k0, c,params.OMEGA, params.direction);
%           Hx_inc_theo = 1./(1i*params.omega*(params.mu0*params.mr_out))*Hx_inc_theo;
%           Hy_inc_theo = 1./(1i*params.omega*(params.mu0*params.mr_out))*Hy_inc_theo;
% 
%            output_mat = {mean_I_fil,mean_I_fil_E_only,Hx_inc_theo,Hy_inc_theo,E_inc_theo , mean_Im_fil_x , mean_Im_fil_y,mean_E_fil,mean_Hx_fil,mean_Hy_fil};
%         end
            
           %

% output_mat = {mean_I_fil,mean(E_sol_final)*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R),Hx_inc_theo,Hy_inc_theo,E_inc_theo};













E_SOL_st_fil_TM = -1;
Hx_SOL_st_fil_TM = -1;
Hy_SOL_st_fil_TM = -1;

    full_fields = {E_SOL_st_fil_TM,Hx_SOL_st_fil_TM,Hy_SOL_st_fil_TM};

Hx_sol_temp = zeros(size(E_sol));
Hy_sol_temp = zeros(size(E_sol));

if (params.calc_full_sol)
    
    X = X + params.sca_x;
    Y = Y + params.sca_y;
    loc_line = [X(:),Y(:)]';
    
    for tt = 1:length(filaments_total)
        if (filaments_total(tt,4) == 1)
          
%             E_sol = E_sol + 1i*omega*(mu0*mr_out)*I_solution(tt)*scalar_green(filaments_total(tt,1:2)',loc_line,n_out,k0,c,OMEGA);
            E_sol = E_sol +I_solution(tt)* scalar_green(filaments_total(tt,1:2)',loc_line,params,0);
            [Hx_temp, Hy_temp] = dyiadic_green(filaments_total(tt,1:2)',loc_line,params,0);
            Hx_sol_temp = Hx_sol_temp + I_solution(tt)* Hx_temp;
            Hy_sol_temp = Hy_sol_temp + I_solution(tt)* Hy_temp;

        end
    end
    
    out_loc = ones(size(X));
    for i = 1:length(params.sca_x)
        
        out_loc = and(out_loc,(X-params.sca_x(i)).^2 +(Y-params.sca_y(i)).^2>R*R);
        in_loc{i} = ((X-params.sca_x(i)).^2 +(Y-params.sca_y(i)).^2)<=R*R;
        
    end
    
    E_sol ( out_loc == 0) = 0;
    Hx_sol_temp ( out_loc == 0) = 0;
    Hy_sol_temp ( out_loc == 0) = 0;

    out_loc_line = [X(out_loc==1),Y(out_loc==1)]';
    
    for tt = 1:length(filaments_total)
        if (filaments_total(tt,4) == 0)
            in_loc_line = [X(in_loc{filaments_total(tt,3)}==1),Y(in_loc{filaments_total(tt,3)}==1)]';
%             E_sol(in_loc{filaments_total(tt,3)}) = E_sol(in_loc{filaments_total(tt,3)})  + 1i*omega*(mu0*mr_in)*I_solution(tt)*scalar_green(filaments_total(tt,1:2)',in_loc_line,n_in,k0,c,OMEGA);
            E_sol(in_loc{filaments_total(tt,3)}) = E_sol(in_loc{filaments_total(tt,3)})  +I_solution(tt)*scalar_green(filaments_total(tt,1:2)',in_loc_line,params,1);
            [Hx_temp, Hy_temp] =  dyiadic_green(filaments_total(tt,1:2)',in_loc_line,params,1);
            
             Hx_sol_temp(in_loc{filaments_total(tt,3)}) = Hx_sol_temp(in_loc{filaments_total(tt,3)}) + I_solution(tt)*Hx_temp;
             Hy_sol_temp(in_loc{filaments_total(tt,3)}) = Hy_sol_temp(in_loc{filaments_total(tt,3)}) + I_solution(tt)*Hy_temp;
            
        end
    end
    
    % % %     for i = 1:length(params.sca_x)
    % % %
    % % %         mean_E_fil(i) = mean(E_sol(in_loc{i}));
    % % %
    % % %     end
    
    % out_loc = (X.*X+Y.*Y)>R*R;
    % in_loc = (X.*X+Y.*Y)<=R*R;
    %
    
    % in_loc_line = [X(in_loc==1),Y(in_loc==1)]';
    %
    % E_sol = zeros(1,length(loc_line));
    %
    % for tt=1:length(in_sources_loc)
    %     E_sol(out_loc(:)) = E_sol(out_loc(:))+1i*omega*(mu0*mr_out)*I_solution(tt)*scalar_green(in_sources_loc(:,tt),out_loc_line,n_out,k0,c,OMEGA);
    % end
    %
    % for tt=1:length(out_sources_loc)
    %     E_sol(in_loc(:)) = E_sol(in_loc(:))+1i*omega*(mu0*mr_in)*I_solution(tt+length(in_sources_loc))*scalar_green(out_sources_loc(:,tt),in_loc_line,n_in,k0,c,OMEGA);
    % end
    
    % tic
    % inc_field = zeros(length(out_loc_line),3);
    % for i=1:length(out_loc_line)
    %     inc_field(i,:) =  double(subs(Esym,[xsym,ysym] ,[out_loc_line(1,i),out_loc_line(2,i)]));
    %
    %     if mod(i,100)==0
    %         est_t =toc;
    %         (length(out_loc_line)-i)/100*est_t
    %         tic
    %     end
    % end
    
    in_loc_temp = not(out_loc);
    in_loc_line = [X(in_loc_temp==1),Y(in_loc_temp==1)]';
    
    if (params.is_plane_wave)
        
        E_sol_total(out_loc(:)) = E_sol(out_loc(:)) +  (exp(-1i*k0*out_loc_line(1,:)));
        
        Hx_sol_total(out_loc(:)) = Hx_sol_temp(out_loc(:)) +  0;
        Hy_sol_total(out_loc(:)) = Hy_sol_temp(out_loc(:)) +  1i/(params.omega*(params.mu0*params.mr_out))*(-1i*params.k0*exp(-1i*params.k0*out_loc_line(1,:)));
    else
%         E_sol_total(out_loc(:)) = E_sol(out_loc(:)) +  scalar_green([params.source_loc_x;params.source_loc_y],[X(out_loc(:)),Y(out_loc(:))]',params.n_out,k0,c,params.OMEGA);
        E_sol_total(out_loc(:)) = E_sol(out_loc(:)) +   scalar_green([params.source_loc_x;params.source_loc_y],[X(out_loc(:)),Y(out_loc(:))]',params, 0);
       [Hx_temp, Hy_temp] =  dyiadic_green([params.source_loc_x;params.source_loc_y],[X(out_loc(:)),Y(out_loc(:))]',params,0);

        Hx_sol_total(out_loc(:)) = Hx_sol_temp(out_loc(:)) +  Hx_temp;
        Hy_sol_total(out_loc(:)) = Hy_sol_temp(out_loc(:)) +  Hy_temp;

    end
    
    E_sol_total(in_loc_temp(:)) = E_sol(in_loc_temp(:));
    
      Hx_sol_total(in_loc_temp(:)) = Hx_sol_temp(in_loc_temp(:));
      Hy_sol_total(in_loc_temp(:)) = Hy_sol_temp(in_loc_temp(:));

    
    E_sol_sca(out_loc(:)) = E_sol(out_loc(:));
    E_sol_sca(in_loc_temp(:)) = E_sol(in_loc_temp(:))-  exp(-1i*k0*in_loc_line(1,:));
    %
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol),'filled'); axis equal
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol),'filled');axis equal
    %
    %
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol_maybe),'filled'); axis equal
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol_maybe),'filled');axis equal
    
    % figure; imagesc(abs(reshape(E_sol_maybe,length(X),length(X)))); axis equal; title('filaments -- ABS')
    %
    % figure; imagesc(imag(reshape(E_sol_maybe,length(X),length(X)))); axis equal ; title('filaments --IMAG')
    % figure; imagesc(real(reshape(E_sol_maybe,length(X),length(X)))); axis equal ;  title('filaments --REAL')
    E_SOL_st_fil_TM = reshape(E_sol_total,length(X),length(X));
    
    Hx_SOL_st_fil_TM = reshape(Hx_sol_total,length(X),length(X));
    Hy_SOL_st_fil_TM = reshape(Hy_sol_total,length(X),length(X));

    full_fields = {E_SOL_st_fil_TM,Hx_SOL_st_fil_TM,Hy_SOL_st_fil_TM};
    
    E_sca = reshape(E_sol_sca,length(X),length(X));
end
end
