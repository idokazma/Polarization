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


if (false)
    for i = 1:length(test_point_total)
        E_inc(i,:) = double(subs(params.Esym,[xsym,ysym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),k0]));
        H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),0,k0]));
        
    end
    disp ('End Eval Symbolic');
    
end

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
        %         OMEGA_orig = params.OMEGA;
        %         params.OMEGA = 0;
        if params.multisources
                E_inc_green = scalar_green([params.sources_v2_x;params.sources_v2_y],test_point_total',params, 0);
                [H_inc_green_x , H_inc_green_y] = dyiadic_green([params.sources_v2_x;params.sources_v2_y],test_point_total',params, 0);
                H_inc_green_x = H_inc_green_x;
                H_inc_green_y = H_inc_green_y;

                E_inc_z = E_inc_green;
                E_inc_x = 0*E_inc_z;
                E_inc_y = 0*E_inc_z;
                
                H_inc_x = H_inc_green_x;
                H_inc_y = H_inc_green_y;
                H_inc_z = 0*H_inc_green_x;
                
                [mean_Hx_inc_center,mean_Hy_inc_center] = dyiadic_green([params.sources_v2_x;params.sources_v2_y],[params.sca_x;params.sca_y],params, 0);
                mean_Ez_inc_center = scalar_green([params.sources_v2_x;params.sources_v2_y],[params.sca_x;params.sca_y],params, 0);
                
                E_inc_z = E_inc_z*0 + mean_Ez_inc_center;
                H_inc_x = H_inc_x*0 + mean_Hx_inc_center;
                H_inc_y = H_inc_y*0 + mean_Hy_inc_center;

        else
            
            E_inc_green = scalar_green([params.source_loc_x;params.source_loc_y],test_point_total',params, 0);
            [H_inc_green_x , H_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],test_point_total',params, 0);
            H_inc_green_x = H_inc_green_x;
            H_inc_green_y = H_inc_green_y;
%             
%             E_inc_green_2 = scalar_green([params.source_loc_x_2;params.source_loc_y_2],test_point_total',params, 0);
%             [H_inc_green_x_2 , H_inc_green_y_2] = dyiadic_green([params.source_loc_x_2;params.source_loc_y_2],test_point_total',params, 0);
%             H_inc_green_x = H_inc_green_x;
%             H_inc_green_y = H_inc_green_y;
            %         params.OMEGA = OMEGA_orig
            
            
            E_inc_z = E_inc_green;
            E_inc_x = 0*E_inc_z;
            E_inc_y = 0*E_inc_z;
            
            H_inc_x = H_inc_green_x;
            H_inc_y = H_inc_green_y;
            H_inc_z = 0*H_inc_green_x;
            
            [mean_Hx_inc_center,mean_Hy_inc_center] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.sca_x;params.sca_y],params, 0);
            mean_Ez_inc_center = scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x;params.sca_y],params, 0);
            
%                             [mean_Hx_inc_center,mean_Hy_inc_center] = dyiadic_green([params.sources_v2_x;params.sources_v2_y],[params.sca_x;params.sca_y],params, 0);
%                 mean_Ez_inc_center = scalar_green([params.sources_v2_x;params.sources_v2_y],[params.sca_x;params.sca_y],params, 0);
%                 
%                 E_inc_z = E_inc_z*1 -1* mean(E_inc_z);
%                 H_inc_x = H_inc_x*1 + 0*mean_Hx_inc_center;
%                 H_inc_y = H_inc_y*1 + 0*mean_Hy_inc_center;

        
        end
        %
        %         E_inc_z = E_inc_z*0 +1*mean_Ez_inc_center;
        %         H_inc_x = H_inc_x*0 +1*mean_Hx_inc_center;
        %         H_inc_y = H_inc_y*0 +1*mean_Hy_inc_center;
        %
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


Mx = zeros(length(test_point_total),length(filaments_total));
MHz = zeros(length(test_point_total),length(filaments_total));

for i=1:length(test_point_total)
    nhat = test_point_total(i,4:6);
    test_point_xy = [test_point_total(i,1:2)];
    
    activations = zeros(length(test_point_total),length(filaments_total));
    activations(i,:) =  (filaments_total(:,3)==test_point_total(i,3)) .*  (filaments_total(:,4)==0);
    Mx(i,activations(i,:)==1) = -scalar_green(filaments_total(activations(i,:)==1,1:2)',test_point_xy',params,1);
    [Gx,Gy] = dyiadic_green(filaments_total(activations(i,:)==1,1:2)',test_point_xy',params,1);
    sol = cross(kron(nhat,ones(length(Gx),1)).',[Gx;Gy;0*Gx]);
    MHz(i,activations(i,:)==1) = -sol(3,:);
    
    activations(i,:) =  2*(filaments_total(:,3)==test_point_total(i,3)) .*  (filaments_total(:,4)==1);
    Mx(i,activations(i,:)==2) = scalar_green(filaments_total(activations(i,:)==2,1:2)',test_point_xy',params,0);
    [Gx,Gy] = dyiadic_green(filaments_total(activations(i,:)==2,1:2)',test_point_xy',params,0);
    sol = cross(kron(nhat,ones(length(Gx),1)).',[Gx;Gy;0*Gx]);
    MHz(i,activations(i,:)==2) = sol(3,:);
    
    activations(i,:) =  3*(filaments_total(:,3)~=test_point_total(i,3)) .*  (filaments_total(:,4)==1);
    Mx(i,activations(i,:)==3) = scalar_green(filaments_total(activations(i,:)==3,1:2)',test_point_xy',params,0);
    [Gx,Gy] = dyiadic_green(filaments_total(activations(i,:)==3,1:2)',test_point_xy',params,0);
    sol = cross(kron(nhat,ones(length(Gx),1)).',[Gx;Gy;0*Gx]);
    MHz(i,activations(i,:)==3) = sol(3,:);
    
    
    Vex(i) = -E_inc_z(i);
    sol = cross(nhat,[H_inc_x(i),H_inc_y(i),0]);
    Vhz(i) = -1*sol(3);
    
    %     Vex(i) = -scalar_green([params.source_loc_x;params.source_loc_y],test_point_xy',params,0);
    %     [Gx,Gy] = dyiadic_green([params.source_loc_x;params.source_loc_y],test_point_xy',params,0);
    %     sol = cross(nhat,[Gx;Gy;0]);
    %     Vhz(i) = -1*sol(3);
    
end


M_total = [Mx;MHz];
I_solution = linsolve(M_total,transpose([Vex,Vhz]));




X = params.X;
Y = params.Y;
loc_line = [X(:),Y(:)]';
E_sol = zeros(1,length(loc_line));

grid_inside_each_sca = ((X).^2 +(Y).^2)<=0.999999*R*R;
grid_inside_each_sca = ((X).^2 +(Y).^2)<=1*R*R;
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

if ( params.hit_plane ==0)
    hit_field_inside = scalar_green([params.sources_v2_x;params.sources_v2_y],in_loc_line,params,0);
else
    hit_field_inside = inc_Ez_field([params.sources_v2_x;params.sources_v2_y],in_loc_line,n_out,k0,c,OMEGA, params.direction);
end
alpha_Fil = mean(mean_I_fil)./mean(hit_field_inside);



for mm = 1 : length(params.sca_x)
    
    in_loc_line = [X(grid_inside_each_sca)+params.sca_x(mm),Y(grid_inside_each_sca)+params.sca_y(mm)]';
    
    for i = 1:length(in_loc_line)
        if (params.hit_plane == 0)
            E_inc_green_inside(i,:) = scalar_green([params.sources_v2_x;params.sources_v2_y],[in_loc_line(1,i),in_loc_line(2,i)]',params,0);
            E_inc_green_inside(i,:) =E_inc_green_inside(i,:);
            
            [H_inc_green_x_temp , H_inc_green_y_temp] = dyiadic_green([params.sources_v2_x;params.sources_v2_y],[in_loc_line(1,i),in_loc_line(2,i)]',params, 0);
            H_inc_green_x_inside(i,:) = H_inc_green_x_temp;
            H_inc_green_y_inside(i,:) = H_inc_green_y_temp;
            
%             
%             [H_inc_green_x_temp_2 , H_inc_green_y_temp_2] = dyiadic_green([params.source_loc_x_2;params.source_loc_y_2],[in_loc_line(1,i),in_loc_line(2,i)]',params, 0);
%             H_inc_green_x_inside(i,:) = H_inc_green_x_temp+params.I2*H_inc_green_x_temp_2;
%             H_inc_green_y_inside(i,:) = H_inc_green_y_temp+params.I2*H_inc_green_y_temp_2;
            
            %             E_inc_green_2 = scalar_green([params.source_loc_x_2;params.source_loc_y_2],test_point_total',params, 0);
            %             [H_inc_green_x_2 , H_inc_green_y_2] = dyiadic_green([params.source_loc_x_2;params.source_loc_y_2],test_point_total',params, 0);
            
            
        else
            E_inc_green_inside(i,:) = inc_Ez_field([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params.n_out,params.k0,params.c,params.OMEGA, params.direction);
            
            [H_inc_green_x_temp , H_inc_green_y_temp] = inc_Ht_field([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params.n_out,params.k0, c,params.OMEGA, params.direction);
            H_inc_green_x_inside(i,:) = 1./(1i*params.omega*(params.mu0*params.mr_out))*H_inc_green_x_temp;
            H_inc_green_y_inside(i,:) = 1./(1i*params.omega*(params.mu0*params.mr_out))*H_inc_green_y_temp;
        end
        
    end
    E_inc_z_inside_sca(mm,:) = E_inc_green_inside(:);
    H_inc_x_inside_sca(mm,:) = H_inc_green_x_inside(:);
    H_inc_y_inside_sca(mm,:) = H_inc_green_y_inside(:);
end

mean_Hx_fil = mean(Hx_sol_final,2);
mean_Hy_fil = mean(Hy_sol_final,2);

mean_Hx_inc = mean(H_inc_x_inside_sca,2);
mean_Hy_inc = mean(H_inc_y_inside_sca,2);

mean_Ez_inc = mean(E_inc_z_inside_sca,2);

[mean_Hx_inc_center,mean_Hy_inc_center] = dyiadic_green([params.sources_v2_x;params.sources_v2_y],[params.sca_x;params.sca_y],params, 0);
mean_Ez_inc_center = scalar_green([params.sources_v2_x;params.sources_v2_y],[params.sca_x;params.sca_y],params, 0);

E_sol_final_point = zeros(length(params.sca_x),1);
Hx_sol_final_point = zeros(length(params.sca_x),1);
Hy_sol_final_point = zeros(length(params.sca_x),1);


for tt = 1:length(filaments_total)
    if (filaments_total(tt,4) == 0)
        E_sol_final_point(filaments_total(tt,3),:) = E_sol_final_point(filaments_total(tt,3),:) + sum(I_solution(tt)*scalar_green(filaments_total(tt,1:2)',[params.sca_x;params.sca_y],params,1));
        [Hx,Hy] = dyiadic_green(filaments_total(tt,1:2)',[params.sca_x;params.sca_y],params,1);
        Hx_sol_final_point(filaments_total(tt,3),:) = Hx_sol_final_point(filaments_total(tt,3),:) + sum(I_solution(tt)* Hx);
        Hy_sol_final_point(filaments_total(tt,3),:) = Hy_sol_final_point(filaments_total(tt,3),:) + sum(I_solution(tt)* Hy);
        
        
    end
end


if length(params.sca_x)==1
    
    for tt = 1 : length(Hy_sol_final)
        
        final_current_H(tt) = -Hx_sol_final(tt)*in_loc_line(1,tt) -  Hy_sol_final(tt)*in_loc_line(2,tt);
        rotation_current_H(tt) =  -H_inc_x_inside_sca(tt)*in_loc_line(1,tt) -  H_inc_y_inside_sca(tt)*in_loc_line(2,tt);
        Jm_fil_x(tt) = (E_sol_final(tt) - E_inc_z_inside_sca(tt))*in_loc_line(1,tt);
        Jm_fil_y(tt) = (E_sol_final(tt) - E_inc_z_inside_sca(tt))*in_loc_line(2,tt);
        
    end
else
    final_current_H = 0;
    rotation_current_H = 0;
    Jm_fil_x = 0;
    Jm_fil_y = 0;
end

H_current_add = (final_current_H-rotation_current_H);

mean_Im_fil_x = -1i*params.omega*params.OMEGA/(params.c^2)* params.sca_x * (mean(E_sol_final) - mean_Ez_inc) *pi*R*R;
mean_Im_fil_y = -1i*params.omega*params.OMEGA/(params.c^2)* params.sca_y * (mean(E_sol_final)  - mean_Ez_inc) *pi*R*R;



mean_I_fil = mean(E_sol_final)*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R) + 1i*params.omega/(params.c^2)*params.OMEGA*mean(H_current_add)*pi*R*R;
mean_I_fil_E_only = mean(E_sol_final)*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R);
output_mat = {mean_I_fil,mean_I_fil_E_only,mean_Hx_inc,mean_Hy_inc,mean_Ez_inc , mean_Im_fil_x , mean_Im_fil_y,mean_E_fil,mean_Hx_fil,mean_Hy_fil , E_sol_final_point, Hx_sol_final_point,Hy_sol_final_point, mean_Ez_inc_center, mean_Hx_inc_center, mean_Hy_inc_center};

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
            E_sol(in_loc{filaments_total(tt,3)}) = E_sol(in_loc{filaments_total(tt,3)})  +I_solution(tt)*scalar_green(filaments_total(tt,1:2)',in_loc_line,params,1);
            [Hx_temp, Hy_temp] =  dyiadic_green(filaments_total(tt,1:2)',in_loc_line,params,1);
            
            Hx_sol_temp(in_loc{filaments_total(tt,3)}) = Hx_sol_temp(in_loc{filaments_total(tt,3)}) + I_solution(tt)*Hx_temp;
            Hy_sol_temp(in_loc{filaments_total(tt,3)}) = Hy_sol_temp(in_loc{filaments_total(tt,3)}) + I_solution(tt)*Hy_temp;
            
        end
    end
    
    in_loc_temp = not(out_loc);
    in_loc_line = [X(in_loc_temp==1),Y(in_loc_temp==1)]';
    
    if (params.is_plane_wave)
        
        E_sol_total(out_loc(:)) = E_sol(out_loc(:)) +  (exp(-1i*k0*out_loc_line(1,:)));
        
        Hx_sol_total(out_loc(:)) = Hx_sol_temp(out_loc(:)) +  0;
        Hy_sol_total(out_loc(:)) = Hy_sol_temp(out_loc(:)) +  1i/(params.omega*(params.mu0*params.mr_out))*(-1i*params.k0*exp(-1i*params.k0*out_loc_line(1,:)));
    else
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
    
    E_SOL_st_fil_TM = reshape(E_sol_total,length(X),length(X));
    
    Hx_SOL_st_fil_TM = reshape(Hx_sol_total,length(X),length(X));
    Hy_SOL_st_fil_TM = reshape(Hy_sol_total,length(X),length(X));
    
    full_fields = {E_SOL_st_fil_TM,Hx_SOL_st_fil_TM,Hy_SOL_st_fil_TM};
    
    E_sca = reshape(E_sol_sca,length(X),length(X));
end
end
