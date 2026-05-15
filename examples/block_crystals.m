function block_crystals(overrides)
% block_crystals  Spectral sweep over a rotating block-crystal array.
% Loads scatterer positions from data/arraypoints.mat and computes the
% transmitted/scattered field at fixed observation points across wavelengths.
%
% Pass an `overrides` struct to tweak any default knob below, e.g.
%
%   block_crystals(struct('wavelengths', (0.5:0.01:1.5)*1e-6))

if nargin < 1, overrides = struct(); end

defaults = struct( ...
    'eps_vec',                11.4, ...
    'radius_ratio',           1/10, ...
    'Iz',                     1, ...
    'is_plane_wave',          0, ...
    'plane_wave_direction',   1, ...
    'len_n',                  81, ...
    'wid_n',                  81, ...
    'grid_extent_factor',     2.1, ...
    'array_scale',            0.75, ...        % multiplier applied to ArrayPoints
    'wavelengths',            (0.2:0.005:3.8)*1e-6, ...
    'observed_points',        [-3, 0, 3; 5, 5, 5] * 1e-6, ...
    'source_distance',        3e-6, ...
    'OMEGA',                  0);
cfg = merge_defaults(defaults, overrides);

data_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'arraypoints.mat');
if ~isfile(data_file)
    error('block_crystals:missingData', ...
          ['Required data file not found: %s\n' ...
           'See data/README.md for how to obtain or generate this file.'], data_file);
end
loaded = load(data_file);
ArrayPoints = loaded.ArrayPoints;

params.te = 0;
params.tm = 1;

params.Iz                    = cfg.Iz;
params.multiscatterer        = 1;
params.is_plane_wave         = cfg.is_plane_wave;
params.calc_full_sol         = 0;
params.plane                 = 0;
params.hit_plane             = 0;
params.len_n                 = cfg.len_n;
params.wid_n                 = cfg.wid_n;
params.len                   = cfg.radius_ratio * 1e-6 * cfg.grid_extent_factor;
params.wid                   = cfg.radius_ratio * 1e-6 * cfg.grid_extent_factor;
params.plane_wave_direction  = cfg.plane_wave_direction;

params.sca_x     = ArrayPoints(:, 1).' * 1e-6 * cfg.array_scale;
params.sca_y     = ArrayPoints(:, 2).' * 1e-6 * cfg.array_scale;
params.OMEGA_vec = cfg.OMEGA;
params.shift_vec = 0;
params.OMEGA     = cfg.OMEGA;

sources = cfg.source_distance * exp(1i*2*pi*linspace(0, 1, 5));
sources = sources(end-1);

observed_point = cfg.observed_points;
E_sol = zeros(size(observed_point, 2), length(cfg.wavelengths));
E_hit = zeros(size(observed_point, 2), length(cfg.wavelengths));

for w = 1:length(cfg.wavelengths)
    params.lambda = cfg.wavelengths(w);
    params        = generate_parameters(params);

    for t = 1:length(sources)
        params.source_loc_x = real(sources(t));
        params.source_loc_y = imag(sources(t));

        params.er_in  = cfg.eps_vec;
        params.n_in   = sqrt(params.mr_in * params.er_in);
        params.radius = cfg.radius_ratio * 1e-6;
        params.OMEGA  = cfg.OMEGA;

        fprintf('Wavelength %d/%d\n', w, length(cfg.wavelengths));

        if params.is_plane_wave
            params = generate_plane_wave(params);
        else
            params = generate_source_wave(params);
        end

        [Pvec, ~, alpha_Pol] = RotatingArray_2D_TM(params, params.E_inc_z, 0);
        disp(abs(alpha_Pol));

        for pp = 1:size(observed_point, 2)
            E_hit(pp, w) = scalar_green( ...
                [params.source_loc_x; params.source_loc_y], ...
                observed_point(:, pp), params, 0);
            for tt = 1:length(Pvec)
                E_sol(pp, w) = E_sol(pp, w) ...
                    - params.omega*params.mu0/4/(1i/4) * Pvec(tt) ...
                    * scalar_green([params.sca_x(tt); params.sca_y(tt)], observed_point(:, pp), params, 0);
            end
        end
    end
end

figure;
plot(cfg.wavelengths, 10*log10(abs(E_hit + E_sol) ./ abs(E_hit)));
xlabel('\lambda [m]'); ylabel('|E_{total}|/|E_{hit}| [dB]');
grid on; grid minor;

end
