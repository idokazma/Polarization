function single_scatterer(overrides)
% single_scatterer  Single rotating dielectric cylinder in TM mode.
% Compares MoM, Filaments, and rotating-frame Polarizability theory across a
% sweep of rotation rates and source positions.
%
% Pass an `overrides` struct to tweak any default knob below, e.g.
%
%   single_scatterer(struct('lambda', 0.8e-6, 'omega_factors', 0))

if nargin < 1, overrides = struct(); end

defaults = struct( ...
    'eps_vec',                11.4, ...        % epsilon of the scatterer
    'radius_ratio',           1/100, ...       % radius in units of wavelength
    'lambda',                 1e-6, ...        % wavelength [m]
    'Iz',                     1, ...           % source current amplitude
    'is_plane_wave',          0, ...           % 1 -> plane wave, 0 -> line source
    'plane_wave_direction',   1, ...           % 1 = +x, 0 = +y (TE only)
    'len_n',                  61, ...          % grid samples per side
    'wid_n',                  61, ...
    'grid_extent_factor',     2.1, ...         % grid extent = radius * factor
    'omega_factors',          -3e-5:1e-5:3e-5, ... % OMEGA_vec = omega * this
    'shift_factors',          [1,2,4,10,20,50,100,200], ... % shift_vec = lambda * this
    'source_distance_lambda', 200, ...         % source ring radius in wavelengths
    'n_sources',              9);              % number of source positions
cfg = merge_defaults(defaults, overrides);

params.te = 0;
params.tm = 1;
params.lambda                = cfg.lambda;
params.Iz                    = cfg.Iz;
params.OMEGA                 = 0;
params.multiscatterer        = 0;
params.is_plane_wave         = cfg.is_plane_wave;
params.calc_full_sol         = 0;
params.plane                 = 0;
params.hit_plane             = 0;
params.len_n                 = cfg.len_n;
params.wid_n                 = cfg.wid_n;
params.len                   = cfg.radius_ratio * cfg.lambda * cfg.grid_extent_factor;
params.wid                   = cfg.radius_ratio * cfg.lambda * cfg.grid_extent_factor;
params.plane_wave_direction  = cfg.plane_wave_direction;

params = generate_parameters(params);

params.OMEGA_vec = params.omega * cfg.omega_factors;
params.shift_vec = cfg.shift_factors * params.lambda;

sources = cfg.source_distance_lambda * params.lambda * ...
          exp(1i*2*pi*linspace(0, 1, cfg.n_sources + 1));
sources = sources(1:end-1);

ni = length(params.shift_vec);
nj = length(params.OMEGA_vec);
nt = length(sources);
results = cell(ni, nj, nt);

counter    = 0;
total_runs = ni*nj*nt;
for t = 1:nt
    for i = 1:ni
        params.er_in = cfg.eps_vec;
        params.n_in  = sqrt(params.mr_in * params.er_in);
        for j = 1:nj
            params.radius       = cfg.radius_ratio * params.lambda;
            params.OMEGA        = params.OMEGA_vec(j);
            params.sca_x        = params.shift_vec(i);
            params.sca_y        = 0;
            params.source_loc_x = real(sources(t)) + params.sca_x;
            params.source_loc_y = imag(sources(t)) + params.sca_y;

            counter = counter + 1;
            fprintf('Calculating... Overall: %2.1f%% Shift: %2.1f%%, Omega: %2.1f%%, Source: %2.1f%%.\n', ...
                100*counter/total_runs, 100*i/ni, 100*j/nj, 100*t/nt);

            if params.is_plane_wave
                params = generate_plane_wave(params);
            else
                params = generate_source_wave(params);
            end

            results{i,j,t} = eval_TM_results(params);
        end
    end
end

TM_alpha_presentation_fields(results);

end
