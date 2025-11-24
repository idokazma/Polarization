%% Pol in Rot
% This code was written by Ido Kazma as part of a master thesis
% @ Tel aviv university under the supervision of Prof. Ben z. Steinberg.
% 
% The goal of this code is to simulate scattering problems in a rotating
% frame of reference and predict the excited fields and how they behave in
% a through rotation.

%% Generate scenario parameters

addpath(genpath(pwd))
% clear all

params.te = 0;
params.tm = 1;

eps_vec = [11.4]; %epsilon of the scatterer
radius_vec = [2/100];%radius of scatterer in terms of wavelengths
params.lambda = 1e-6; %[meters] %wavelength of radiation
params.Iz = 1; % single source radiating current

params.multiscatterer = 0;

params.is_plane_wave = 0;
params.calc_full_sol = 0;
params.plane = 0;
params.hit_plane=0;   

% geometric parameters of the grid for calculation.

params.len_n = 51;
params.wid_n = 51;
params.len = radius_vec*params.lambda*2.1;
params.wid = radius_vec*params.lambda*2.1;

params.plave_wave_direction = 1; %1 for x, 0 for y;

if (params.multiscatterer == 1)
    GA_generator;
    x_shift = params.lambda*30;
    y_shift = params.lambda*0;

    params.sca_x = 100*params.lambda*VogelArrayXY(1,1:10)+x_shift;
    params.sca_y = 100*params.lambda*VogelArrayXY(2,1:10)+y_shift;
end

generate_parameters;
% params.only_center=1;
params.OMEGA_vec = 0+1e-0*params.omega*(-5e-5:(2*5e-6):5e-5);
% params.shift_vec = 1*[0,20, 40,60,80,100] * params.lambda;
params.shift_vec = 1*[0,25,50,75,100] * params.lambda;
% params.shift_vec = 1*[0,50,100] * params.lambda;

%% Generate sources locations
dis = 50;  %distance of source from scattereres
sources = dis * lambda * exp(1i*(2*pi*linspace(0,1,5))) ;
sources = sources(1:end-1);


rng('default') 
sources = [sources, sources/10];
sources_2 = sources ;
params.sources_v2 = sources + 0*0.01*1i* lambda;
% sources_2 = -sources ;
% number_of_sources = 1;
% params.sources_v2 = lambda*(20*rand(number_of_sources,length(sources))+dis).*exp(1i*(2*pi*rand(number_of_sources,length(sources))));
% params.sources_v2 = [params.sources_v2,-params.sources_v2];
% params.sources_v2 = [params.sources_v2,conj(params.sources_v2)];

% sources = [sources, conj(sources)];
% sources = [sources, -real(sources) + 1i*imag(sources)];
% sources = [sources, -(sources)];
% sources_2 = [sources_2, sources_2/10];
% sources_2 = [sources_2, conj(sources_2)];
% sources_2 = [sources_2, -real(sources_2) + 1i*imag(sources_2)];
% sources_2 = [sources_2, -sources_2];

% dis = 10.98472;  %distance of source from scattereres
% sources_ = dis *lambda * exp(1i*2*pi*linspace(0,1, 4)) ;
% sources = [sources_(1:end-1), sources];
now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
Run_name = ['RUN_' , now_str];

%% START SIMULATION
% Hitting Field
counter = 0;
total_runs = length(params.shift_vec)*length(params.OMEGA_vec)*length(sources);

for t = 1:length(sources) % for every source
                    params.I2 = 0;

    for i = 1:length(params.shift_vec) % for every shift of scatterer location
%               params.er_in = eps_vec;
%               params.n_in = sqrt(params.mr_in*params.er_in);
        %     params.n_out = sqrt(params.mr_out*params.er_out);
        
        for j = 1:length(params.OMEGA_vec) % for every Omega (rotation rate)
            
            params.radius = radius_vec*lambda; %reupdate the initial scatterer radius
            params.OMEGA = params.OMEGA_vec(j);
            
            params.sca_x = 1*params.shift_vec(i);
            params.sca_y = 0*params.shift_vec(i);
            
            params.source_loc_x = real(sources(t))+params.sca_x;
            params.source_loc_y = imag(sources(t))+params.sca_y;

            params.sources_v2_x = real(params.sources_v2(:,t))+params.sca_x;          
            params.sources_v2_y = imag(params.sources_v2(:,t))+params.sca_y;          
            
            
            params.source_loc_x_2 = real(sources_2(t))+params.sca_x;
            params.source_loc_y_2 = imag(sources_2(t))+params.sca_y;
            counter=counter+1;
            formatSpec = 'Calculating...Overall: %2.3f%% Shift: %2.1f%%, Omega: %2.1f%%, Source %2.1f%%.\n';
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
