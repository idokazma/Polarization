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
radius_vec = [1/100];
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
params.hit_plane=0;   

params.Iz = 1;

params.len_n = 901;
params.wid_n = 901;
params.len = 1*0.25e-6;
params.wid = 1*0.25e-6;

params.plave_wave_direction = 1; %1 for x, 0 for y;


generate_parameters;

% params.OMEGA = 10e-3*params.omega;
% params.OMEGA_vec = 1*linspace(-5e-6*params.omega,5e-6*params.omega,15);
params.OMEGA_vec = 1*params.omega*(0:1e-5:3e-5);

% params.OMEGA_vec = (1e-8*params.omega*(0:1:20));
%   params.OMEGA_vec = 5e-6*params.omega;

%params.OMEGA_vec = 1e-5*params.omega;

params.shift_vec = 1*[5,10,20] * params.lambda;
% params.shift_vec = 0;

dis = 1000;
fprintf('.\n.\nOMEGA*rho_c/c = %f \n.\n.\n', max(params.OMEGA_vec)*max(params.shift_vec)/3e8)

sources = dis *lambda * exp(1i*linspace(0,2*pi, 5));
sources = sources(1:end-1);
% sources = [sources,sources/10];
now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');

Run_name = ['RUN_' , now_str];
% OMEGA_vec = 1 * 20e-6*params.omega;

%% START SIMULATION
% Hitting Field

for t = 1:length(sources)
    
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
            
            params.sca_x = params.lambda*[ 0 ] + 1*params.shift_vec(i) + x_shift;
            params.sca_y = params.lambda*[ 0 ] + 0*params.shift_vec(i) + y_shift;
            
            params.source_loc_x = real(sources(t));
            params.source_loc_y = imag(sources(t));
%             switch t
%                 case 1
%                 params.source_loc_x = params.lambda*[cos(2*pi/6)*dis] + 0*params.shift_vec(i)+ x_shift;
%                 params.source_loc_y = params.lambda*[sin(2*pi/6)*dis] + y_shift;
%                 
%                 params.source_loc_x = params.lambda*[0] + 0*params.shift_vec(i)+ x_shift;
%                 params.source_loc_y = params.lambda*[-dis] + y_shift;
%                 
%                 params.direction = -1i;
% 
%                 case 2
%                 params.source_loc_x = params.lambda*[cos(2*pi/6)*dis] + 0*params.shift_vec(i)+ x_shift;
%                 params.source_loc_y = params.lambda*[-sin(2*pi/6)*dis] + y_shift;
%                 
%                  params.source_loc_x = params.lambda*[0] + 0*params.shift_vec(i)+ x_shift;
%                  params.source_loc_y = params.lambda*[dis] + y_shift;
%                  params.direction = 1i;
% 
%                 case 3
%                 params.source_loc_x = -1*norm([params.source_loc_x,params.source_loc_y])+ 0*params.lambda*[-dis] + 0*params.shift_vec(i)+ x_shift;
%                 params.source_loc_y = params.lambda*[0] + y_shift;
%                 params.direction = 1;
%                 
%                 case 4 
%                 params.direction = -1;
%             end  
            
                % %             
%             params.source_loc_x = params.lambda*[0] + 0*params.shift_vec(i)+ x_shift;
%             params.source_loc_y = params.lambda*[-dis] + y_shift;
            
            
            

            
%              params.source_loc_x = params.lambda*[0] + 0*params.shift_vec(i)+ x_shift;
%             params.source_loc_y = params.lambda*[dis] + y_shift;
% %             
%             params.sca_x = params.lambda*[ 0 ] + 1*params.shift_vec(i)+ x_shift;
%             params.sca_y = params.lambda*[ 0 ] + 0*params.shift_vec(i)+ y_shift;
%           
%             params.source_loc_x = params.lambda*[10] + x_shift;
%             params.source_loc_y = params.lambda*[100] + y_shift;
                
                
                %             params.sca_x = params.lambda*[ 0 ] + x_shift;
                %             params.sca_y = params.lambda*[ 10 ] + y_shift;
            
            
            
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
    end
    
  save(Run_name); disp('~~~~~~~~~SAVED~~~~~~~~~~~~');  
end


now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');

if (params.tm)
 TM_alpha_presentation_fields_multiple;
%     TM_alpha_presentation;
end
if (params.te)
    TE_alpha_presentation;
end    

% save (['alpha_TE_filaments_Im_source_100_',now_str]);

% PRINT_RESULTS_WITHOUTMIE;
% PRINT_RESULTS;


