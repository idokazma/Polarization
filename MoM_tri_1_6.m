function [SOL, mean_E,mean_I, effective_radius , alpha_mom,mean_E0] = MoM_tri_1_6(params,E_inc_z,H_inc_x,H_inc_y,Esym)
tic
%
% params.radius = 1/50*1e-6;
% params.len = 2*params.radius;
% params.wid = 2*params.radius;
% params.len_n = 30; %number of cells
% params.wid_n = 30;
%
%
% % params.len = 5*1e-6;
% % params.wid = 5*1e-6;
% params.lambda = 1*1e-6;
% % params.center = 1*[1,0; 1*cos(2*pi/3),1*sin(2*pi/3); 1*cos(4*pi/3),1*sin(4*pi/3)]*1e-6;
% % phase_sca = 0:1/5:1;
% % phase_sca=phase_sca(1:end-1)*2*pi;
% % params.center = [cos(phase_sca);sin(phase_sca)]'*1e-6;
% params.center = pos'*1e-6;
% % params.radius = 1/50*1e-6;
% params.delta_eps = 10.4;
% params.omegaratio = Rot;

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

loc =find(sqrt((X_mom).^2 + (Y_mom).^2)<params.radius);


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


% scatterers_x = X(loc);
% scatterers_y = Y(loc);
scatterers_x = [];
scatterers_y = [];
for pp=1:length(params.sca_x)
    scatterers_x = [scatterers_x; X_mom(loc)+params.sca_x(pp)];
    scatterers_y = [scatterers_y; Y_mom(loc)+params.sca_y(pp)];
    
end


%% source gen
E_bg_dy = zeros(size(scatterers_x));
E_bg_st = zeros(size(scatterers_x));

% source_location = [5*ones(1,15);linspace(3,-3,15)];%;72.2,76.3];
source_location_x = params.source_loc_x;
source_location_y = params.source_loc_y;



c= params.c;
omega = params.omega;

Omega_B=params.OMEGA;
u0 = 4*pi*1e-7;

if (~params.is_plane_wave)
    if(params.hit_plane==0)
        E_bg_st =  scalar_green_mom([source_location_x;source_location_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA);
    else
       
%         E_bg_st =  scalar_green_mom([source_location_x;source_location_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA);
        E_bg_st =  inc_Ez_field([source_location_x;source_location_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA, params.direction);

    end
else
%     E_bg_st =  1+0*exp(-1i*k0*scatterers_x);
    E_bg_st =  1*exp(-1i*k0*scatterers_x);
end
I_amps=1;
% for(ll=1:size(source_location,2))
%
%
%     E_bg_dy = E_bg_dy + 1*1i/4*besselh(0,1,k0.*sqrt((scatterers_x-source_location(1,ll)).^2+(scatterers_y-source_location(2,ll)).^2)).*exp(1i*Omega_B/omega*k0^2.*(scatterers_x*source_location(2,ll)-scatterers_y*source_location(1,ll)));
%     E_bg_st = E_bg_st + 1*1i/4*besselh(0,1,k0.*sqrt((scatterers_x-source_location(1,ll)).^2+(scatterers_y-source_location(2,ll)).^2));
%
% end
syms xsym ysym zsym k0sym
Xsym = [xsym ysym zsym];


% tic
% for i = 1:length(scatterers_x)
%     E_inc(i,:) = double(subs(Esym,[xsym,ysym,k0sym] ,[scatterers_x(i),scatterers_y(i),k0]));
% end
% toc
%
% E_inc_x = E_inc(:,1);
% E_inc_y = E_inc(:,2);

% toc
% E_bg_st = E_bg_st/abs(max(max(E_bg_st)));
% E_bg_dy = E_bg_dy/abs(max(max(E_bg_dy)));
% figure; imagesc(real(E_bg_st)); set(gca,'YDir','normal'); axis equal;

%% scatter gen
% scattereres = (epsilon_bg+delta_eps)*ones(size(scatterers_y));

scattereres = (epsilon_bg+delta_eps)*ones(size(scatterers_y));

%
% center = [0,0];
% radius = 5.2;
% loc =find(and(((X-center(1)).^2 + (Y-center(2)).^2)<radius^2,((X-center(1)).^2 + (Y-center(2)).^2)>(radius-0.4)^2));

%  center = params.center;
%  radius = params.radius;

% E_middle_st = interp2(X,Y,E_bg_st,center(1),center(2));
% E_middle_dy = interp2(X,Y,E_bg_dy,center(1),center(2));

%  E_bg_st = E_bg_st/abs(E_middle_st);
% E_bg_dy = E_bg_dy/abs(E_middle_dy);




% scattereres (3,3) = 2;
% scattereres (3,500) = 2;

%  E_bg_dy_line = reshape(E_bg_dy,1,[]);
%  E_bg_st_line = reshape(E_bg_st,1,[]);

E_bg_dy_line = E_bg_dy(:);
E_bg_st_line = E_bg_st(:);

scattereres_line = reshape(scattereres,1,[]);

scatter_ind = find(scattereres_line~=epsilon_bg);
scatter_bg = find(scattereres_line==epsilon_bg);

%% solve for scatteres
% h = waitbar(0,'Good things happen for those who wait...');

% computations take place here

rho = zeros(length(scatter_ind),length(scatter_ind));
cross_rho = zeros(length(scatter_ind),length(scatter_ind));
GG = cross_rho;
rho_bg = zeros(length(scatter_bg),length(scatter_ind));
cross_rho_bg = zeros(length(scatter_bg),length(scatter_ind));

for (ii=1:length(scatter_ind))
    if (mod(ii,200)==0)
%         waitbar(ii / length(scatter_ind))
    end
    %     rho(:,ii) = sqrt((scatterers_x-scatterers_x(ii)).^2+(scatterers_y-scatterers_y(ii)).^2);
    %     cross_rho(:,ii) = scatterers_x*scatterers_y(ii)-scatterers_y*scatterers_x(ii);
    
    GG(:,ii) = scalar_green_mom([scatterers_x(ii);scatterers_y(ii)],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA);
    
    
    %      rho_bg(:,ii) = sqrt((X(scatter_bg)-X(scatter_ind(ii))).^2+(Y(scatter_bg)-Y(scatter_ind(ii))).^2);
    %      cross_rho_bg(:,ii) = X(scatter_bg)*Y(scatter_ind(ii))-Y(scatter_bg)*X(scatter_ind(ii));
    
end

% close(h)

% % % G_mat_st = -1i/4*besselh(0,1,k0.*rho)*V_i*k0^2*delta_eps;
% % % % toc
% % % % disp ('Green calc done 1');
% % % G_mat_dy = G_mat_st.*exp(1i*Omega_B/omega*k0^2.*cross_rho);

G_mat_dy = -GG*V_i*k0^2*delta_eps;

clear rho
clear cross_rho

for (ii=1:length(scatter_ind))
    
         G_mat_dy(ii,ii) = 1-(1i*pi/2*beta*gamma)*k0^2*delta_eps;
%     G_mat_dy(ii,ii) = 1/CylPolTM_IDO (params, R_i_eff,params.k0,params.er_out,params.er_in)*(pi*R_i_eff^2*params.omega*(params.er_in-params.er_out)*params.e0)/1i;
    %   G_mat_st(ii,ii) = 1-(1i*pi/2*beta*gamma)*k0^2*delta_eps;
end


% toc
% solution_mat_dy = (G_mat_dy\transpose(E_bg_dy_line(scatter_ind))).';
% solution_mat_st = (G_mat_st\transpose(E_bg_st_line(scatter_ind))).';
% %CHECKKKKKKKKKKKKKKKKKKKKKKK
% solution_mat_dy = (G_mat_dy\(E_bg_dy_line(scatter_ind))').';
% solution_mat_st = (G_mat_st\(E_bg_st_line(scatter_ind))').';
% solution_mat_dy = linsolve(G_mat_dy,(E_bg_dy_line(scatter_ind)));


solution_mat_st = linsolve(G_mat_dy,(E_bg_st_line(scatter_ind)));
start_ind = 1;
mean_E = zeros(1,length(params.sca_x));
mean_I = mean_E;
mean_alpha = mean_I;
for ll = 1:length(params.sca_x)
    
    mean_E(ll) = mean(solution_mat_st(start_ind:(start_ind+length(loc)-1)));
    mean_E0(ll) = mean(E_bg_st_line(start_ind:(start_ind+length(loc)-1)));
    mean_I(ll) = mean_E(ll)*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*V_i*length(loc));
    start_ind = (start_ind+length(loc)-1)+1;
    
    mean_alpha(ll) = mean_I(ll)/mean_E0(ll);
end
SOL = -1;
if(params.plane==1)
    alpha_mom = mean_I./(1+0*scalar_green_mom([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA));
else
    alpha_mom = mean_I./(0+1*scalar_green_mom([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA));
end
alpha_mom = mean_I/mean(E_bg_st_line);
fprintf('num of scatterers %d \n', length(scatterers_x))

if (params.calc_full_sol == 1)
    % disp ('Solved Scatterers 2');
    % toc
    
    % BG_mat_st = 0;
    BG_mat_st = 1i/4*besselh(0,1,k0.*rho_bg)*V_i*k0^2*delta_eps;
    % disp ('Solved Background 3');
    % toc
    % BG_mat_dy = 0;
    BG_mat_dy = BG_mat_st.*exp(1i*Omega_B/omega*k0^2.*cross_rho_bg);
    
    clear rho_bg
    clear cross_rho_bg
    %   full_sol_dy(scatter_bg) = E_bg_dy_line(scatter_bg) + BG_mat_dy*(solution_mat_dy);
    % % full_sol_dy(scatter_bg) = 0;
    %  full_sol_dy(scatter_ind) = solution_mat_dy;
    
    full_sol_st(scatter_bg) = E_bg_st_line(scatter_bg) + BG_mat_dy*(solution_mat_st);
    % full_sol_st(scatter_bg) = 0;
    full_sol_st(scatter_ind) = solution_mat_st;
%     figure; scatter(scatterers_x,scatterers_y,10,abs(full_sol_st),'filled');
    
    
    if (~params.is_plane_wave)
        SOL = scalar_green_mom([source_location_x;source_location_y],[X(:),Y(:)]',params.n_out,k0,c,params.OMEGA);
        SOL = reshape(SOL, size(X));
    else
        SOL = exp(-1i*k0*X);
    end
    % SOL = exp(-1i*k0*X);
    % SOL = scalar_green_mom([source_location_x;source_location_y],[X(:),Y(:)]',params.n_out,k0,c,params.OMEGA);
    % SOL = reshape(SOL, size(X));
    ind = [];
    ind_sca = [];
    for i = 1 : length(scatterers_x)
        if (mod(i,200)==0)
            waitbar(i / length(scatter_ind))
        end
        rho = sqrt((X-scatterers_x(i)).*(X-scatterers_x(i)) + (Y - scatterers_y(i)).*(Y - scatterers_y(i)));
        cross_rho =  X*scatterers_y(i)-Y*scatterers_x(i);
        if (find(rho==0))
            ind = [ind, find(rho==0)];
            ind_sca = [ind_sca,i];
        end
        %      SOL(rho~=0) = SOL(rho~=0) + solution_mat_st(i)*(1i/4*besselh(0,1,k0.*rho(rho~=0))*V_i*k0^2*delta_eps).*exp(1i*Omega_B/omega*k0^2.*cross_rho(rho~=0));
        add = solution_mat_st(i)*scalar_green_mom([scatterers_x(i);scatterers_y(i)],[X(rho~=0),Y(rho~=0)]',params.n_out,k0,c,params.OMEGA)*V_i*k0^2*delta_eps;
        SOL(rho~=0) = SOL(rho~=0) + transpose(add);
        SOL_2 = SOL;
%          SOL_2(rho==0) = SOL(rho==0) + solution_mat_st(i)*(-(1i*pi/2*beta*gamma)*k0^2*delta_eps);
%               SOL(rho==0) = SOL(rho==0) + solution_mat_st(i)*(-(1i*pi/2*beta*gamma)*k0^2*delta_eps);
    end
    
     SOL(ind) = solution_mat_st(ind_sca)/V_i*V_i_bg;
    toc
end


%
% epsilon_zero = 8.85418781762e-12;
% I_sol_dy = full_sol_dy*(-1i*omega*epsilon_zero*(delta_eps))*(V_i);
%
% for ii=1:length(pos)
%
%
% % delta_eps = params.delta_eps;
% % current_st(ii)= sum(I_sol_st(loc));
% % mean_Est (ii) = mean(SOL_st(loc));
% % current_dy(ii) = sum(I_sol_dy(loc));
% mean_Edy (ii) = mean(full_sol_dy(1+(ii-1)*length(X(loc)):ii*length(X(loc))));
%
% end
% % SOL_dy = reshape(full_sol_dy,len_n,[]);
% % SOL_st = reshape(full_sol_st,len_n,[]);
%  I_mom_dy = mean_Edy*(-1i*omega*epsilon_zero*(delta_eps))*(V_i*length(loc));%%(pi*(radius^2));
% % I_mom_dy = mean_Edy*(-1i*omega*epsilon_zero*(delta_eps))*(pi*(radius^2));



%
%
% toc
%
% disp ('Full SOL 4');
%
% % figure; imagesc(real(SOL_dy)); set(gca,'YDir','normal'); axis equal; title('Real Dy');
% % figure; imagesc(real(SOL_st)); set(gca,'YDir','normal'); axis equal;title('Real Sta');
% % figure; imagesc((SOL_st).*conj(SOL_st)); set(gca,'YDir','normal'); axis equal;title('ENG Sta');
% % figure; imagesc((SOL_dy).*conj(SOL_dy)); set(gca,'YDir','normal'); axis equal;title('ENG Dy');
% %
% % figure; imagesc(db(SOL_st)); set(gca,'YDir','normal'); axis equal;title('ENG Sta');
% % figure; imagesc(db(SOL_dy)); set(gca,'YDir','normal'); axis equal;title('ENG Dy');
%
% toc
%
% disp ('DONE 5');
% % figure; imagesc(db(SOL_st)); set(gca,'YDir','normal'); axis equal;title('ENG Sta');
% % figure; imagesc(db(SOL_dy)); set(gca,'YDir','normal'); axis equal;title('ENG Dy');
%
% figure; scatter(X(scatter_ind),Y(scatter_ind),[],(abs(SOL_dy(scatter_ind)./SOL_st(scatter_ind))),'filled'); set(gca,'YDir','normal'); axis equal;title('ENG Dy');
% figure; scatter(X(scatter_ind),Y(scatter_ind),[],(real(SOL_dy(scatter_ind))),'filled'); set(gca,'YDir','normal'); axis equal;title('ENG Dy');
% figure; scatter(X(scatter_ind),Y(scatter_ind),[],(abs(SOL_dy(scatter_ind))),'filled'); set(gca,'YDir','normal'); axis equal;title('|E|');
%
% epsilon_zero = 8.8e-12;
% E_in = E_bg_st;
% E_SOL_st = SOL_st;
% E_SOL_dy = SOL_dy;
% I_sol_dy = E_SOL_dy*(-1i*omega*epsilon_zero*(delta_eps))*(V_i);
% I_sol_st = E_SOL_dy*(-1i*omega*epsilon_zero*(delta_eps))*(V_i);
%
% figure; scatter(X(scatter_ind),Y(scatter_ind),[],(abs(I_sol_dy(scatter_ind))),'filled'); set(gca,'YDir','normal'); axis equal;title('ENG Dy');
%
% current_st = sum(SOL_st(loc)*V_i);
% current_dy = sum(SOL_dy(loc)*V_i);
%
% for ii=1:size(center,1)
% loc =find(sqrt((X-center(ii,1)).^2 + (Y-center(ii,2)).^2)<radius);
%
% delta_eps = params.delta_eps;
% current_st(ii)= sum(I_sol_st(loc));
% mean_Est (ii) = mean(SOL_st(loc));
% current_dy(ii) = sum(I_sol_dy(loc));
% mean_Edy (ii) = mean(SOL_dy(loc));
%
% end
% abs(current_dy./current_st-1);
% mean_Edy./mean_Est-1
%
% I_mom_dy = mean_Edy*(-1i*omega*epsilon_zero*(delta_eps))*(pi*(radius^2))
% end
% mean_Edy
% figure; subplot(1,2,1); imagesc(x,y,abs(E_mid).^2); set(gca,'YDir','normal'); axis equal; title ('ENG Mie Series');  axis image; xlabel('x(\lambda)'); ylabel('y(\lambda)'); colorbar
% subplot(1,2,2); imagesc(x,y,abs(SOL_st).^2); set(gca,'YDir','normal'); axis equal; title ('ENG MOM');  axis image; xlabel('x(\lambda)'); ylabel('y(\lambda)'); colorbar