% %% Pol in Rot
% % This code was written by Ido Kazma as part of a master thesis
% % @ Tel aviv university under the supervision of Prof. Ben z. Steinberg.
% % 
% % The goal of this code is to simulate scattering problems in a rotating
% % frame of reference and predict the excited fields and how they behave in
% % a through rotation.
% 
% %% Generate scenario parameters
% 
% params.te = 0;
% params.tm = 1;
% 
% eps_vec = [11.4]; %epsilon of the scatterer
% radius_vec = [1/10];%radius of scatterer in terms of wavelengths
% params.lambda = 0.6e-6; %[meters] %wavelength of radiation
% params.Iz = 1; % single source radiating current
% 
% params.multiscatterer = 1;
% 
% params.is_plane_wave = 0;
% params.calc_full_sol = 0;
% params.plane = 0;
% params.hit_plane=0;   
% 
% % geometric parameters of the grid for calculation.
% 
% params.len_n = 3;
% params.wid_n = 3;
% params.len = radius_vec*1e-6*3.1;
% params.wid = radius_vec*1e-6*3.1;
% 
% params.plave_wave_direction = 1; %1 for x, 0 for y;
% 
% if (params.multiscatterer == 1)
%     GA_generator;
%     x_shift = params.lambda*30;
%     y_shift = params.lambda*0;
% 
%     params.sca_x = params.lambda*VogelArrayXY(1,1:10)+x_shift;
%     params.sca_y = params.lambda*VogelArrayXY(2,1:10)+y_shift;
%     
% 
% end
% 
% params.sca_x = ArrayPoints(:,1)*1e-6;
% params.sca_y = ArrayPoints(:,2)*1e-6;
% 
% generate_parameters;
% 
% params.OMEGA_vec = 1*params.omega*(0:1e-5:3e-5);
% params.shift_vec = 1*[1,2,4,10,20,50,100,200] * params.lambda;
% params.shift_vec = 1*[1,2,4,8,16,32,64,128] * params.lambda;
% 
% %% Generate sources locations
% 
% 
% 
% dis = 200.345; 
% sources = dis *lambda * exp(1i*2*pi*linspace(0,1, 5));
% sources = sources(1:end-1);
% 
% params.sca_x = ArrayPoints(:,1)*1e-6;
% params.sca_y = ArrayPoints(:,2)*1e-6;
% params.OMEGA_vec = 0;
% params.shift_vec = 0;
% 
% dis = 3; 
% sources = dis *1e-6 * exp(1i*2*pi*linspace(0,1, 5));
% sources = sources(end-1);
% now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
% Run_name = ['RUN_' , now_str];
% 
% %% START SIMULATION
% % Hitting Field
% counter = 0;
% total_runs = length(params.shift_vec)*length(params.OMEGA_vec)*length(sources);
% for t = 1:length(sources) % for every source
%     
%     for i = 1:length(params.shift_vec) % for every shift of scatterer location
%               params.er_in = eps_vec;
%               params.n_in = sqrt(params.mr_in*params.er_in);
%         %     params.n_out = sqrt(params.mr_out*params.er_out);
%         
%         for j = 1:length(params.OMEGA_vec) % for every Omega (rotation rate)
%             
%             params.radius = radius_vec*lambda; %reupdate the initial scatterer radius
%             params.OMEGA = params.OMEGA_vec(j);
%             
% %             params.sca_x = params.lambda*[ 0 ] + 1*params.shift_vec(i);
% %             params.sca_y = params.lambda*[ 0 ] + 0*params.shift_vec(i);
% %             
%             params.source_loc_x = real(sources(t))+params.sca_x;
%             params.source_loc_y = imag(sources(t))+params.sca_y;
% 
%             counter=counter+1;
%             formatSpec = 'Calculating...Overall: %2.1f%% Shift: %2.1f%%, Omega: %2.1f%%, Source %2.1f%%.\n';
%             fprintf(formatSpec,100*counter/total_runs,100*i/length(params.shift_vec),100*j/length(params.OMEGA_vec), 100*t/length(sources))
%             
%             if (params.is_plane_wave ==1)
%                 generate_plane_wave;
%             else
%                 generate_source_wave;
%             end
% 
%             if (params.tm)
%                 EVAL_TM_RESULTS;
%             end
%             
%             if (params.te)
%                 params.radius = radius_vec*lambda;
%                 EVAL_TE_RESULTS;
%             end
%         end
%     end
%     
%   save(Run_name); disp('~~~~~~~~~SAVED~~~~~~~~~~~~');  
% end
% 
% 
% now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
% 
% if (params.tm)
%  TM_alpha_presentation_fields_multiple;
% end
% if (params.te)
%     TE_alpha_presentation;
% end    
% 
% % save (['alpha_TE_filaments_Im_source_100_',now_str]);
% 
% % PRINT_RESULTS_WITHOUTMIE;
% % PRINT_RESULTS;
% 
% 
%% Pol in Rot
% This code was written by Ido Kazma as part of a master thesis
% @ Tel aviv university under the supervision of Prof. Ben z. Steinberg.
% 
% The goal of this code is to simulate scattering problems in a rotating
% frame of reference and predict the excited fields and how they behave in
% a through rotation.

%% Generate scenario parameters

addpath(genpath(pwd))


params.te = 0;
params.tm = 1;

eps_vec = [11.4]; %epsilon of the scatterer
radius_vec = [1/100];%radius of scatterer in terms of wavelengths
params.lambda = 1e-6; %[meters] %wavelength of radiation
params.Iz = 1; % single source radiating current

params.multiscatterer = 0;

params.is_plane_wave = 0;
params.calc_full_sol = 0;
params.plane = 0;
params.hit_plane=0;   

% geometric parameters of the grid for calculation.

params.len_n = 61;
params.wid_n = 61;
params.len = radius_vec*params.lambda*2.1;
params.wid = radius_vec*params.lambda*2.1;

params.plave_wave_direction = 1; %1 for x, 0 for y;

if (params.multiscatterer == 1)
    GA_generator;
    x_shift = params.lambda*30;
    y_shift = params.lambda*0;

    params.sca_x = params.lambda*VogelArrayXY(1,1:10)+x_shift;
    params.sca_y = params.lambda*VogelArrayXY(2,1:10)+y_shift;
end

generate_parameters;

params.OMEGA_vec = 1*params.omega*(0:1e-5:3e-5);
params.shift_vec = 1*[1,2,4,10,20,50,100,200] * params.lambda;
params.shift_vec = 1*[1,2,4,8,16,32,64,128] * params.lambda;

%% Generate sources locations
dis = 200.345; 
sources = dis *lambda * exp(1i*2*pi*linspace(0,1, 5));
sources = sources(1:end-1);



now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
Run_name = ['RUN_' , now_str];

%% START SIMULATION
% Hitting Field
counter = 0;
total_runs = length(params.shift_vec)*length(params.OMEGA_vec)*length(sources);
for t = 1:length(sources) % for every source
    
    for i = 1:length(params.shift_vec) % for every shift of scatterer location
              params.er_in = eps_vec;
              params.n_in = sqrt(params.mr_in*params.er_in);
        %     params.n_out = sqrt(params.mr_out*params.er_out);
        
        for j = 1:length(params.OMEGA_vec) % for every Omega (rotation rate)
            
            params.radius = radius_vec*lambda; %reupdate the initial scatterer radius
            params.OMEGA = params.OMEGA_vec(j);
            
            params.sca_x = params.lambda*[ 0 ] + 1*params.shift_vec(i);
            params.sca_y = params.lambda*[ 0 ] + 0*params.shift_vec(i);
            
            params.source_loc_x = real(sources(t))+params.sca_x;
            params.source_loc_y = imag(sources(t))+params.sca_y;

            counter=counter+1;
            formatSpec = 'Calculating...Overall: %2.1f%% Shift: %2.1f%%, Omega: %2.1f%%, Source %2.1f%%.\n';
            fprintf(formatSpec,100*counter/total_runs,100*i/length(params.shift_vec),100*j/length(params.OMEGA_vec), 100*t/length(sources))
            
            if (params.is_plane_wave ==1)
                generate_plane_wave;
            else
                generate_source_wave;
            end

            if (params.tm)
                EVAL_TM_RESULTS;
            end
            
            if (params.te)
                params.radius = radius_vec*lambda;
                EVAL_TE_RESULTS;
            end
        end
    end
    
  save(Run_name); disp('~~~~~~~~~SAVED~~~~~~~~~~~~');  
end


now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');

if (params.tm)
 TM_alpha_presentation_fields_multiple;
end
if (params.te)
    TE_alpha_presentation;
end    

% save (['alpha_TE_filaments_Im_source_100_',now_str]);

% PRINT_RESULTS_WITHOUTMIE;
% PRINT_RESULTS;