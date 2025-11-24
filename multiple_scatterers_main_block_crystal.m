%% Compare All Methods -
addpath(genpath(pwd))


GA_generator;
params.te = 0;
params.tm = 1;

eps_vec = [11.4];
% radius_vec = [2/100];
radius_vec = [1/10];
params.lambda = 1e-6; %[meters]

params.is_plane_wave = 0;
params.calc_full_sol = 0;
params.plane = 0;
params.num_of_scatterers = 15;
% params.len_n = 180;
% params.wid_n = 180;
% params.len = 0.5*1e-6;
% params.wid = 0.5*1e-6;

params.len_n = 4;
params.wid_n = 4;
params.len = radius_vec*params.lambda*2.1;
params.wid = radius_vec*params.lambda*2.1;

params.Iz = 1;
params.hit_plane = 0
params.plave_wave_direction = 1; %1 for x, 0 for y;

params.sca_x = 1*params.lambda*VogelArrayXY(1,1:params.num_of_scatterers);
params.sca_y = 1*params.lambda*VogelArrayXY(2,1:params.num_of_scatterers);

load('arraypoints.mat')
% params.block_x = ArrayPoints(:,1)*1e-6*0.75;
% params.block_y = ArrayPoints(:,2)*1e-6*0.75;

% r=(max(sqrt(params.sca_x.^2 + params.sca_y.^2))+10*params.lambda):

x_all = -(max(sqrt(params.sca_x.^2 + params.sca_y.^2))+10*params.lambda):(params.lambda/2*0.75):(max(sqrt(params.sca_x.^2 + params.sca_y.^2))+10*params.lambda);
y_all = -(max(sqrt(params.sca_x.^2 + params.sca_y.^2))+10*params.lambda):(params.lambda/2*0.75):(max(sqrt(params.sca_x.^2 + params.sca_y.^2))+10*params.lambda);
[AX, AY] = meshgrid(x_all,y_all);

loc_filter = (sqrt(AX.^2 + AY.^2) < (max(sqrt(params.sca_x.^2 + params.sca_y.^2))+4*params.lambda))&(sqrt(AX.^2 + AY.^2) > (max(sqrt(params.sca_x.^2 + params.sca_y.^2))+2.5*params.lambda));

AX_block = AX(loc_filter);
AY_block = AY(loc_filter);

clear AX AY x_all y_all
temp_sca_x = [params.sca_x,AX_block.'];
temp_sca_y = [params.sca_y,AY_block.'];

array_r = sqrt(max(params.sca_x.^2 + params.sca_y.^2));
temp_r = sqrt(max(temp_sca_x.^2 + temp_sca_y.^2));


params.sca_x = [params.sca_x,AX_block.'];
params.sca_y = [params.sca_y,AY_block.'];

% params.sca_x = params.sca_x*(array_r/temp_r);
% params.sca_y = params.sca_y*(array_r/temp_r);

disp (max(sqrt(params.sca_x.^2+params.sca_y.^2)))
generate_parameters;

params.OMEGA_vec = 100*params.omega*(-5e-7:2e-8:5e-7);
params.shift_vec = [0000] * params.lambda;


x_shift = 000*params.lambda;
y_shift = 0*params.lambda;

params.sca_x = params.sca_x+x_shift;
%% START SIMULATION
% Hitting Field
for t = 1:1
    for i = 1:1
        for j = 1:length(params.OMEGA_vec)
            

            params.radius = radius_vec*lambda;
            params.OMEGA = params.OMEGA_vec(j);
            
%             params.sca_x = params.lambda*VogelArrayXY(1,1:params.num_of_scatterers) + x_shift;
%             params.sca_y = params.lambda*VogelArrayXY(2,1:params.num_of_scatterers) + y_shift;
            
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
                EVAL_TM_RESULTS_only_pol;
                toc
            end
            
            if (params.te)
                tic
                params.radius = radius_vec*lambda;
                EVAL_TE_RESULTS;
                toc
            end
            
            disp(j/length(params.OMEGA_vec)*100)
            disp(-(j-length(params.OMEGA_vec))*toc/60);
        end
    end
end




now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');

if (params.tm)
    plt_pol_block;
end
% if (params.te)
%     TE_alpha_presentation;
% end


