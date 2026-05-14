function block_crystals()
% block_crystals  Spectral sweep over a rotating block-crystal array.
% Loads scatterer positions from data/arraypoints.mat and computes the
% transmitted/scattered field at fixed observation points across wavelengths.

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

eps_vec      = 11.4;
radius_ratio = 1/10;
params.Iz    = 1;

params.multiscatterer = 1;
params.is_plane_wave  = 0;
params.calc_full_sol  = 0;
params.plane          = 0;
params.hit_plane      = 0;

params.len_n = 81;
params.wid_n = 81;
params.len   = radius_ratio*1e-6*2.1;
params.wid   = radius_ratio*1e-6*2.1;

params.plane_wave_direction = 1;

params.sca_x = ArrayPoints(:,1).'*1e-6*0.75;
params.sca_y = ArrayPoints(:,2).'*1e-6*0.75;
params.OMEGA_vec = 0;
params.shift_vec = 0;
params.OMEGA     = 0;

dis     = 3e-6;
sources = dis * exp(1i*2*pi*linspace(0,1,5));
sources = sources(end-1);

params.wavelengths = (0.2:0.005:3.8)*1e-6;
observed_point = [-3, 0, 3; 5, 5, 5] * 1e-6;

E_sol = zeros(size(observed_point,2), length(params.wavelengths));
E_hit = zeros(size(observed_point,2), length(params.wavelengths));

for w = 1:length(params.wavelengths)
    params.lambda = params.wavelengths(w);
    params = generate_parameters(params);

    for t = 1:length(sources)
        params.source_loc_x = real(sources(t));
        params.source_loc_y = imag(sources(t));

        params.er_in  = eps_vec;
        params.n_in   = sqrt(params.mr_in*params.er_in);
        params.radius = radius_ratio*1e-6;
        params.OMEGA  = params.OMEGA_vec(1);

        fprintf('Wavelength %d/%d\n', w, length(params.wavelengths));

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
plot(params.wavelengths, 10*log10(abs(E_hit + E_sol) ./ abs(E_hit)));
xlabel('\lambda [m]'); ylabel('|E_{total}|/|E_{hit}| [dB]');
grid on; grid minor;

end
