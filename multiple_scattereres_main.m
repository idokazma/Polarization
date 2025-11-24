%% Compare All Methods -
addpath(genpath(pwd))


GA_generator;
params.te = 0;
params.tm = 1;

eps_vec = [11.4];
radius_vec = [2/100];
params.lambda = 1e-6; %[meters]

params.is_plane_wave = 0;
params.calc_full_sol = 0;
params.plane = 0;
params.num_of_scatterers = 15;
% params.len_n = 180;
% params.wid_n = 180;
% params.len = 0.5*1e-6;
% params.wid = 0.5*1e-6;

params.len_n = 18;
params.wid_n = 18;
params.len = radius_vec*params.lambda*2.1;
params.wid = radius_vec*params.lambda*2.1;

params.Iz = 1;
params.hit_plane = 0;
params.toplot=0;
params.plave_wave_direction = 1; %1 for x, 0 for y;

params.sca_x = 100*params.lambda*VogelArrayXY(1,1:params.num_of_scatterers);
params.sca_y = 100*params.lambda*VogelArrayXY(2,1:params.num_of_scatterers);



disp (max(sqrt(params.sca_x.^2+params.sca_y.^2)))
generate_parameters;

params.OMEGA_vec = 1*1e-3*params.omega*(-5e-4:1e-5:5e-4);
params.shift_vec = [10000] * params.lambda;


x_shift = 0*params.lambda;
y_shift = 500*params.lambda;
%% START SIMULATION
% Hitting Field
for t = 1:1
    for i = 1:1
        for j = 1:length(params.OMEGA_vec)
            
            disp(j/length(params.OMEGA_vec)*100)
            disp(-(j-length(params.OMEGA_vec))*toc/60);
            params.radius = radius_vec*lambda;
            params.OMEGA = params.OMEGA_vec(j);
            
            params.sca_x = 20*params.lambda*VogelArrayXY(1,1:params.num_of_scatterers) + x_shift;
            params.sca_y = 20*params.lambda*VogelArrayXY(2,1:params.num_of_scatterers) + y_shift;

            params.plave_wave_direction = 1;
            
            params.source_loc_x = x_shift;
            params.source_loc_y = y_shift;
            
            params.sources_v2_x = params.source_loc_x;          
            params.sources_v2_y = params.source_loc_y;  
            

            if (params.is_plane_wave ==1)
                generate_plane_wave;
            else
                generate_source_wave;
            end
            
            if (params.tm)
                tic
                EVAL_TM_RESULTS;
                toc
            end
            
            if (params.te)
                tic
                params.radius = radius_vec*lambda;
                EVAL_TE_RESULTS;
                toc
            end
        end
    end
end




now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');

if (params.tm)
    plot_fil_mom_pol;
end
if (params.te)
    compare_fil_pol_TE;
    TE_alpha_presentation;
end


