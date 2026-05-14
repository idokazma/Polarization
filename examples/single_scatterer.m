function single_scatterer()
% single_scatterer  Single rotating dielectric cylinder in TM mode.
% Compares MoM, Filaments, and rotating-frame Polarizability theory across a
% sweep of rotation rates OMEGA and source distances.

params.te = 0;
params.tm = 1;

eps_vec      = 11.4;
radius_ratio = 1/100;
params.lambda = 1e-6;
params.Iz     = 1;
params.OMEGA  = 0;

params.multiscatterer = 0;
params.is_plane_wave  = 0;
params.calc_full_sol  = 0;
params.plane          = 0;
params.hit_plane      = 0;

params.len_n = 61;
params.wid_n = 61;
params.len   = radius_ratio*params.lambda*2.1;
params.wid   = radius_ratio*params.lambda*2.1;

params.plane_wave_direction = 1;

params = generate_parameters(params);

params.OMEGA_vec = params.omega*(-3e-5:1e-5:3e-5);
params.shift_vec = [1,2,4,10,20,50,100,200] * params.lambda;

dis     = 200;
sources = dis*params.lambda * exp(1i*2*pi*linspace(0,1,10));
sources = sources(1:end-1);

ni = length(params.shift_vec);
nj = length(params.OMEGA_vec);
nt = length(sources);
results = cell(ni, nj, nt);

counter    = 0;
total_runs = ni*nj*nt;
for t = 1:nt
    for i = 1:ni
        params.er_in = eps_vec;
        params.n_in  = sqrt(params.mr_in*params.er_in);
        for j = 1:nj
            params.radius = radius_ratio*params.lambda;
            params.OMEGA  = params.OMEGA_vec(j);
            params.sca_x  = params.shift_vec(i);
            params.sca_y  = 0;
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
