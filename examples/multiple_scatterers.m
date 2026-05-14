function multiple_scatterers(overrides)
% multiple_scatterers  Compare MoM, Filaments, and Polarizability theory for
% a Vogel-spiral array of dielectric cylinders in TM mode under rotation.
%
% Pass an `overrides` struct to tweak any default knob below, e.g.
%
%   multiple_scatterers(struct('n_scatterers', 20, 'eps_vec', 4))

if nargin < 1, overrides = struct(); end

defaults = struct( ...
    'eps_vec',              11.4, ...
    'radius_ratio',         2/100, ...
    'lambda',                1e-6, ...
    'Iz',                    1, ...
    'is_plane_wave',         0, ...
    'plane_wave_direction',  1, ...
    'len_n',                 200, ...
    'wid_n',                 200, ...
    'len',                   0.5e-6, ...
    'wid',                   0.5e-6, ...
    'n_scatterers',          15, ...
    'omega_factors',         1e-1*(-5e-4:5e-5:5e-4), ...
    'shift_factors',         0);
cfg = merge_defaults(defaults, overrides);

params.te = 0;
params.tm = 1;
params.lambda                = cfg.lambda;
params.Iz                    = cfg.Iz;
params.OMEGA                 = 0;
params.is_plane_wave         = cfg.is_plane_wave;
params.calc_full_sol         = 0;
params.plane                 = 0;
params.hit_plane             = 0;
params.len_n                 = cfg.len_n;
params.wid_n                 = cfg.wid_n;
params.len                   = cfg.len;
params.wid                   = cfg.wid;
params.plane_wave_direction  = cfg.plane_wave_direction;

VogelArrayXY = GA_generator(max(cfg.n_scatterers, 100));
params.sca_x = params.lambda * VogelArrayXY(1, 1:cfg.n_scatterers);
params.sca_y = params.lambda * VogelArrayXY(2, 1:cfg.n_scatterers);

disp(max(sqrt(params.sca_x.^2 + params.sca_y.^2)));

params = generate_parameters(params);

params.OMEGA_vec = params.omega * cfg.omega_factors;
params.shift_vec = cfg.shift_factors * params.lambda;

fprintf('OMEGA*rho_c/c        = %f\n', max(params.OMEGA_vec)*max(params.shift_vec)/3e8);
fprintf('OMEGA*max(rho_c)/c   = %f\n', max(params.sca_x.^2 + params.sca_y.^2)*max(params.OMEGA_vec)/3e8);

ni = length(params.shift_vec);
nj = length(params.OMEGA_vec);
results = cell(ni, nj);

for i = 1:ni
    params.er_in = cfg.eps_vec;
    params.n_in  = sqrt(params.mr_in * params.er_in);
    for j = 1:nj
        params.radius = cfg.radius_ratio * params.lambda;
        params.OMEGA  = params.OMEGA_vec(j);
        params.sca_x  = params.lambda * VogelArrayXY(1, 1:cfg.n_scatterers);
        params.sca_y  = params.lambda * VogelArrayXY(2, 1:cfg.n_scatterers);
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
