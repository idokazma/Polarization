function [SOL, mean_E,mean_I, effective_radius , alpha_mom,mean_E0,E_center] = MoM(params,E_inc_z,H_inc_x,H_inc_y,Esym)

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
        E_bg_st =  1*scalar_green_mom([source_location_x;source_location_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA);
        %         E_bg_st =E_bg_st+params.I2*  scalar_green_mom([params.source_loc_x_2;params.source_loc_y_2],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA);
        if length(params.sca_x)==1
            if params.only_center==1
                E_bg_st = 1*E_bg_st - mean(E_bg_st);
            elseif params.only_center==0
                E_bg_st = 0*E_bg_st + scalar_green_mom([source_location_x;source_location_y],[params.sca_x(pp),params.sca_y(pp)]',params.n_out,k0,c,params.OMEGA);
            end
        end
        %         E_bg_st = 0*E_bg_st + mean(E_bg_st);
        %         E_bg_st = 0*E_bg_st + scalar_green_mom([source_location_x;source_location_y],[params.sca_x(pp),params.sca_y(pp)]',params.n_out,k0,c,params.OMEGA);
        %         E_bg_st = 1*E_bg_st - scalar_green_mom([source_location_x;source_location_y],[params.sca_x(pp),params.sca_y(pp)]',params.n_out,k0,c,params.OMEGA);
        %         E_bg_st =  scalar_green_mom([source_location_x;source_location_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA)- mean(scalar_green_mom([source_location_x;source_location_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA));
        
        if params.multisources
            E_bg_st =scalar_green_mom([params.sources_v2_x;params.sources_v2_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA);
            
        end
    else
        
        E_bg_st =  inc_Ez_field([source_location_x;source_location_y],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA, params.direction);
        
    end
else
    E_bg_st =  1*exp(-1i*k0*scatterers_x);
end
E_bg_dy_line = E_bg_dy(:);
E_bg_st_line = E_bg_st(:);

scattereres = (epsilon_bg+delta_eps)*ones(size(scatterers_y));

scattereres_line = reshape(scattereres,1,[]);

scatter_ind = find(scattereres_line~=epsilon_bg);
scatter_bg = find(scattereres_line==epsilon_bg);


%% solve for scatteres


rho = zeros(length(scatter_ind),length(scatter_ind));
cross_rho = zeros(length(scatter_ind),length(scatter_ind));
GG = cross_rho;
rho_bg = zeros(length(scatter_bg),length(scatter_ind));
cross_rho_bg = zeros(length(scatter_bg),length(scatter_ind));

for (ii=1:length(scatter_ind))
    
    GG(:,ii) = scalar_green_mom([scatterers_x(ii);scatterers_y(ii)],[scatterers_x,scatterers_y]',params.n_out,k0,c,params.OMEGA);
    
end


G_mat_dy = -GG*V_i*k0^2*delta_eps;

clear rho
clear cross_rho
Mi = (1i*pi/2*beta*gamma);
Mi = (1+1i*(2/pi*log(k0*params.n_out*sqrt(V_i)/(4*exp(1)))))*V_i;
for (ii=1:length(scatter_ind))
    
    G_mat_dy(ii,ii) = 1-Mi*k0^2*delta_eps;
end


solution_mat_st = linsolve(G_mat_dy,(E_bg_st_line(scatter_ind)));
start_ind = 1;
mean_E = zeros(1,length(params.sca_x));
mean_I = mean_E;
mean_alpha = mean_I;
E_center = solution_mat_st(find((scatterers_x==params.sca_x) & (scatterers_y==params.sca_y)));

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
    
    BG_mat_st = 1i/4*besselh(0,1,k0.*rho_bg)*V_i*k0^2*delta_eps;
    
    BG_mat_dy = BG_mat_st.*exp(1i*Omega_B/omega*k0^2.*cross_rho_bg);
    
    clear rho_bg
    clear cross_rho_bg
    
    full_sol_st(scatter_bg) = E_bg_st_line(scatter_bg) + BG_mat_dy*(solution_mat_st);
    full_sol_st(scatter_ind) = solution_mat_st;
    
    if (~params.is_plane_wave)
        SOL = scalar_green_mom([source_location_x;source_location_y],[X(:),Y(:)]',params.n_out,k0,c,params.OMEGA);
        SOL = reshape(SOL, size(X));
    else
        SOL = exp(-1i*k0*X);
    end
    
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
        add = solution_mat_st(i)*scalar_green_mom([scatterers_x(i);scatterers_y(i)],[X(rho~=0),Y(rho~=0)]',params.n_out,k0,c,params.OMEGA)*V_i*k0^2*delta_eps;
        SOL(rho~=0) = SOL(rho~=0) + transpose(add);
        SOL_2 = SOL;
    end
    
    SOL(ind) = solution_mat_st(ind_sca)/V_i*V_i_bg;
    toc
end
