%% Pol in Rot
% This code was written by Ido Kazma as part of a master thesis
% @ Tel aviv university under the supervision of Prof. Ben z. Steinberg.
%
% The goal of this code is to simulate scattering problems in a rotating
% frame of reference and predict the excited fields and how they behave in
% a through rotation.
clear all
%% Generate scenario parameters

params.te = 0;
params.tm = 1;

eps_vec = [11.4]; %epsilon of the scatterer
radius_vec = [1/40];%radius of scatterer in terms of wavelengths
params.Iz = 1; % single source radiating current

params.multiscatterer = 1;

params.is_plane_wave = 0;
params.calc_full_sol = 0;
params.plane = 0;
params.hit_plane=0;

% geometric parameters of the grid for calculation.

params.len_n = 81;
params.wid_n = 81;
params.len = radius_vec*1e-6*2.1;
params.wid = radius_vec*1e-6*2.1;

params.plave_wave_direction = 1; %1 for x, 0 for y;

load('arraypoints.mat')

params.sca_x = ArrayPoints(:,1)*1e-6;
params.sca_y = ArrayPoints(:,2)*1e-6;
params.lambda = 1e-6;
generate_parameters;


%% Generate sources locations



params.sca_x = ArrayPoints(:,1).'*1e-6*0.75;
params.sca_y = ArrayPoints(:,2).'*1e-6*0.75;
params.sca_x = (-5:(1/2*0.75):5)*1e-6;
params.sca_y = (0:(1/2*0.75):5)*1e-6;
[params.sca_x,params.sca_y] = meshgrid(params.sca_x,params.sca_y);
params.sca_x = params.sca_x(:).';
params.sca_y = params.sca_y(:).';
params.sca_y = params.sca_y - mean(params.sca_y);
params.sca_x = params.sca_x - mean(params.sca_x);
% figure;
% scatter(params.sca_x,params.sca_y); hold on;
% scatter( params.source_loc_x,  params.source_loc_y);
% scatter( obsereved_point(1,:), obsereved_point(2,:));
% grid minor
% params.sca_x = ArrayPoints(:,1).'*1e-6*1;
% params.sca_y = ArrayPoints(:,2).'*1e-6*1;

params.OMEGA_vec =  [-2,-1,0,1,2,3,4] * 1*1e-5 * params.omega;
params.shift_vec = 0;



dis = 3;
sources = dis  * exp(1i*2*pi*linspace(0,1, 5));
sources = sources(end-1);
now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
Run_name = ['RUN_' , now_str];

%% START SIMULATION
% Hitting Field
counter = 0;
total_runs = length(params.shift_vec)*length(params.OMEGA_vec)*length(sources);

params.wavelengths = (0.2:0.005:3.8)*1e-6; %[meters] %wavelength of radiation
params.wavelengths = (0.68:0.0003:1.2)*1e-6; %[meters] %wavelength of radiation
params.wavelengths = (0.68:0.002:1.2)*1e-6; %[meters] %wavelength of radiation
obsereved_point = [-3,0,3;5,5,5]*1e-6;
obsereved_point = [-3,0,3;5,5,5]* 1e-6;
obsereved_point = [-3,0,3,;5,5,5 ]* 1e-6;
% obsereved_point = [-3,0,3,;5,-5,5 ]* 1e-6;
params.shift_x_rec = 1*5000*1e-6;
params.shift_y_rec = -0*1e-6;
% obsereved_point = [-2,-2,-2,;4,4,4 ]* 1e-6;
% sources = 1*[1, -3]* 1e-6;
fact = 1
results = [];
for j = 1:length(params.OMEGA_vec) % for every Omega (rotation rate)
    obsereved_point = [0,0,0,;-5,-5,-5 ]* 1e-6;
    sources = 1*[1, 6]* 1e-6;
    for rec = 1:2
        E_sol = zeros(length(obsereved_point),length(params.wavelengths));
        E_hit = zeros(length(obsereved_point),length(params.wavelengths));
        dis = 5;
        
        if rec == 2
            obsereved_point =1*[1,1,1,;6,6,6 ]* 1e-6;
            sources =1*[0, -5]* 1e-6;
            %         obsereved_point = obsereved_point*-1
            %         sources = sources*-1;
            %     else
            %         obsereved_point = -1* [-0,0,0,;dis,dis,dis ]* 1e-6;
            %         sources = 1*[0, dis]* 1e-6;
            %
        end
        obsereved_point = obsereved_point + [params.shift_x_rec;params.shift_y_rec]
        sources = sources + [params.shift_x_rec,params.shift_y_rec];
        params.source_loc_x = (sources(1));
        params.source_loc_y = (sources(2));
        params.sca_x = (-5:(1/2*0.75):5)*1e-6;
        params.sca_y = (0:(1/2*0.75):4)*1e-6;
        [params.sca_x,params.sca_y] = meshgrid(params.sca_x,params.sca_y);
        params.sca_x = params.sca_x(:).';
        params.sca_y = params.sca_y(:).';
        params.sca_y = params.sca_y - mean(params.sca_y)+params.shift_y_rec;
        params.sca_x = params.sca_x - mean(params.sca_x)+params.shift_x_rec;
        %     params.sca_x = ArrayPoints(:,1).'*1e-6*0.75 + params.shift_x_rec;
        %     params.sca_y = ArrayPoints(:,2).'*1e-6*0.75;
        %     params.sca_y = params.sca_y - 1*mean(params.sca_y) +params.shift_y_rec ;
        
        for w = 1:length(params.wavelengths)
            total_runs = length(params.wavelengths);
            
            params.lambda = params.wavelengths(w);
            generate_parameters;
            
            
            for t = 1:1 % for every source
                
                
                for i = 1:length(params.shift_vec) % for every shift of scatterer location
                    params.er_in = eps_vec;
                    params.n_in = sqrt(params.mr_in*params.er_in);
                    %     params.n_out = sqrt(params.mr_out*params.er_out);
                    
                    
                    %                 params.radius = radius_vec*lambda; %reupdate the initial scatterer radius
                    params.radius = radius_vec*1e-6; %reupdate the initial scatterer radius
                    params.OMEGA = params.OMEGA_vec(j);
                    
                    %             params.sca_x = params.lambda*[ 0 ] + 1*params.shift_vec(i);
                    %             params.sca_y = params.lambda*[ 0 ] + 0*params.shift_vec(i);
                    %
                    
                    counter=counter+1;
                    formatSpec = 'Calculating...Overall: %2.1f%% Shift: %2.1f%%, Omega: %2.1f%%, Source %2.1f%%.\n';
                    fprintf(formatSpec,100*counter/total_runs,100*i/length(params.shift_vec),100*j/length(params.OMEGA_vec), 100*t/length(sources))
                    
                    if (params.is_plane_wave ==1)
                        generate_plane_wave;
                    else
                        generate_source_wave;
                    end
                    
                    if (params.tm)
                        [Pvec, E_SOL_st_POL, alpha_Pol] = RotatingArray_2D_TM(params, params.E_inc_z, 0);
                        disp(abs(alpha_Pol))
                        
                        
                        for pp = 1:length(obsereved_point)
                            E_hit(pp,w) = scalar_green([params.source_loc_x;params.source_loc_y],obsereved_point(:,pp),params, 0);
                            for tt = 1:length(Pvec)
                                E_sol(pp,w) = E_sol(pp,w) + -params.omega*params.mu0/4/(1i/4)*Pvec(tt) *scalar_green([params.sca_x(tt), params.sca_y(tt)]',obsereved_point(:,pp),params,0);
                            end
                        end
                        
                        
                        
                    end
                    
                    if (params.te)
                        params.radius = radius_vec*lambda;
                        EVAL_TE_RESULTS;
                    end
                end
            end
            
            %         save(Run_name); disp('~~~~SAVED~~~~~');
        end
        now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
        figure(101);
        results = [results; 20*log10(abs(E_hit(2,:)+E_sol(2,:))./abs(E_hit(2,:)))];
        plot(params.wavelengths,20*log10(abs(E_hit(2,:)+E_sol(2,:))./abs(E_hit(2,:))));
        hold on;
        grid minor;
        
        
    end
    
    
end

figure;
scatter(params.sca_x,params.sca_y); hold on;
scatter( params.source_loc_x,  params.source_loc_y,'x');
scatter( obsereved_point(1,:), obsereved_point(2,:));
grid minor


figure()
plot(params.wavelengths,results(4:2:end,:).','-db');
hold on
plot(params.wavelengths,results(3:2:end,:).','-xr');
grid on
grid minor
if (params.tm)
    TM_alpha_presentation_fields_multiple;
end
if (params.te)
    TE_alpha_presentation;
end




% save (['alpha_TE_filaments_Im_source_100_',now_str]);

% PRINT_RESULTS_WITHOUTMIE;
% PRINT_RESULTS;