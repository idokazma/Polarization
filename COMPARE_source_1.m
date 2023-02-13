%% Compare All Methods -
% 3 scatterer
% radius [1/100, 1/50] wavelength
% located {[0,0], [0.2,0] ,[-0.1,0.1]}
% center located source
% mu = 1
% epsilon = 4 , 11
% TM Only
% No rotation
GA_generator;
params.te = 0;
params.tm = 1;

eps_vec = [11.4];
radius_vec = [2/100];
params.lambda = 1e-6; %[meters]

% x_shift = params.lambda*30;
% y_shift = params.lambda*0;
%
% params.source_loc_x = params.lambda*[0]+x_shift;
% params.source_loc_y = params.lambda*[0]+y_shift;
%
% params.sca_x = params.lambda*[ 0 ] + x_shift;
% params.sca_y = params.lambda*[ 1 ]+ y_shift;

%  params.sca_x = params.lambda*VogelArrayXY(1,1:10)+x_shift;
%  params.sca_y = params.lambda*VogelArrayXY(2,1:10)+y_shift;


params.is_plane_wave = 0;
params.calc_full_sol = 0;
params.plane = 0;

params.len_n = 600;
params.wid_n = 600;
params.len = 1*1e-6;
params.wid = 1*1e-6;
params.Iz = 1
params.hit_plane = 0
params.plave_wave_direction = 1; %1 for x, 0 for y;

             params.sca_x = params.lambda*VogelArrayXY(1,1:10);
             params.sca_y = params.lambda*VogelArrayXY(2,1:10);
disp (max(sqrt(params.sca_x.^2+params.sca_y.^2)))
generate_parameters;

% params.OMEGA = 10e-3*params.omega;
params.OMEGA_vec = 1*linspace(-5e-6*params.omega,5e-6*params.omega,15);
params.OMEGA_vec = 1*params.omega*(0:1e-6:4e-6);

params.OMEGA_vec = 1e-3*params.omega*(-2e-4:5e-5:2e-4);
%  params.OMEGA_vec = 0*5e-6*params.omega;

% params.OMEGA_vec = 1e-5*params.omega;

params.shift_vec = [0] * params.lambda;

fprintf('.\n.\nOMEGA*rho_c/c = %f \n.\n.\n', max(params.OMEGA_vec)*max(params.shift_vec)/3e8)
fprintf('.\n.\nOMEGA*max(rho_c)/c = %f \n.\n.\n', max((params.sca_x.^2+params.sca_y.^2))*max(params.OMEGA_vec)/3e8)

% OMEGA_vec = 1 * 20e-6*params.omega;

%% START SIMULATION
% Hitting Field

for t = 1:1
    
    for i = 1:length(params.shift_vec)
              params.er_in = eps_vec;
              params.n_in = sqrt(params.mr_in*params.er_in);
        %     params.n_out = sqrt(params.mr_out*params.er_out);
        %
        for j = 1:length(params.OMEGA_vec)
            
            params.radius = radius_vec*lambda;
            params.OMEGA = params.OMEGA_vec(j);
            params.shift_vec(i)
            
            x_shift = 0*params.shift_vec(i)*params.lambda;
            y_shift = 0*params.lambda;
            
            params.source_loc_x = 0*params.lambda*[0] + 0*params.shift_vec(i)+ x_shift;
            params.source_loc_y = 0*params.lambda*[-1000] + y_shift;
            
%             params.sca_x = params.lambda*[ 0 ] + 1*params.shift_vec(i) + x_shift;
%             params.sca_y = params.lambda*[ 0 ] + y_shift;
            
            if (t == 1)
                params.plave_wave_direction = 1;
            else
                params.plave_wave_direction = 0;
                if (params.tm)
                    params.plane = 0;
                end
            params.source_loc_x = 0*params.lambda*[0] + 0*params.shift_vec(i)+ x_shift;
            params.source_loc_y = 0*params.lambda*[1000] + y_shift;
            
%             params.sca_x = params.lambda*[ 0 ] + 1*params.shift_vec(i)+ x_shift;
%             params.sca_y = params.lambda*[ 0 ]   + y_shift;
            
            if (t==3)
            
                
            params.source_loc_x = 0*params.lambda*[-1000] + 0*params.shift_vec(i)+ x_shift;
            params.source_loc_y = 0*params.lambda*[0] + y_shift;
            
%             params.sca_x = params.lambda*[ 0 ] +  1* params.shift_vec(i)+ x_shift;
%             params.sca_y = params.lambda*[ 0 ]   + y_shift;
            end
%             params.source_loc_x = params.lambda*[10] + x_shift;
%             params.source_loc_y = params.lambda*[100] + y_shift;
                
                
                %             params.sca_x = params.lambda*[ 0 ] + x_shift;
                %             params.sca_y = params.lambda*[ 10 ] + y_shift;
            end
            
            
            if (params.is_plane_wave ==1)
                generate_plane_wave;
            else
                generate_source_wave;
            end
            
            
            
            
            %MIE SERIES
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
       ppp=1  ;  
 
    end
    
end


now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');

if (params.tm)
    TM_alpha_presentation;
end
if (params.te)
    TE_alpha_presentation;
end    

% save (['alpha_TE_filaments_Im_source_100_',now_str]);

% PRINT_RESULTS_WITHOUTMIE;
% PRINT_RESULTS;


