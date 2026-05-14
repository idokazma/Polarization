function multiple_scatterers()
% multiple_scatterers  Compare MoM, Filaments, and Polarizability theory for
% a Vogel-spiral array of dielectric cylinders in TM mode under rotation.

params.te = 0;
params.tm = 1;

eps_vec      = 11.4;
radius_ratio = 2/100;
params.lambda = 1e-6;
params.Iz     = 1;
params.OMEGA  = 0;

params.is_plane_wave = 0;
params.calc_full_sol = 0;
params.plane         = 0;
params.hit_plane     = 0;

params.len_n = 200;
params.wid_n = 200;
params.len   = 0.5e-6;
params.wid   = 0.5e-6;

params.plane_wave_direction = 1;

VogelArrayXY = GA_generator(20);
params.sca_x = params.lambda*VogelArrayXY(1, 1:15);
params.sca_y = params.lambda*VogelArrayXY(2, 1:15);

disp(max(sqrt(params.sca_x.^2 + params.sca_y.^2)));

params = generate_parameters(params);

params.OMEGA_vec = 1e-1*params.omega*(-5e-4:5e-5:5e-4);
params.shift_vec = 0 * params.lambda;

fprintf('OMEGA*rho_c/c = %f\n', max(params.OMEGA_vec)*max(params.shift_vec)/3e8);
fprintf('OMEGA*max(rho_c)/c = %f\n', max(params.sca_x.^2 + params.sca_y.^2)*max(params.OMEGA_vec)/3e8);

ni = length(params.shift_vec);
nj = length(params.OMEGA_vec);
results = cell(ni, nj);

for i = 1:ni
    params.er_in = eps_vec;
    params.n_in  = sqrt(params.mr_in*params.er_in);
    for j = 1:nj
        params.radius = radius_ratio*params.lambda;
        params.OMEGA  = params.OMEGA_vec(j);
        params.sca_x  = params.lambda*VogelArrayXY(1, 1:15);
        params.sca_y  = params.lambda*VogelArrayXY(2, 1:15);
        params.source_loc_x = 0;
        params.source_loc_y = 0;

        if params.is_plane_wave
            params = generate_plane_wave(params);
        else
            params = generate_source_wave(params);
        end

        tic;
        results{i,j} = eval_TM_results(params);
        toc;
    end
end

plot_fil_mom_pol(results);
TM_alpha_presentation_fields(reshape(results, [size(results), 1]));

end
