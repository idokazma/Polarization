function r = eval_TM_results(params)
% eval_TM_results  Run the three TM methods (MoM, Polarizability, Filaments)
% (plus Mie if applicable) for a single configuration and return a struct
% with all per-method outputs.
%
% Output struct fields:
%   r.params       snapshot of params (with effective_radius substituted)
%   r.mom          .E_sol, .mean_E, .mean_I, .alpha, .mean_E0    (or .ok=false on failure)
%   r.pol          .Pvec, .E_sol, .alpha
%   r.filaments    .E_sol, .mean_E, .mean_I, .alpha, .field_current_mat
%   r.mie          .E_sol, .mean_E       (only set when applicable)

% --- MoM ---
try
    [E_SOL_st_MoM, mean_E_MoM, mean_I_MoM, effective_radius, alpha_MoM, mean_E0_MoM] = MoM(params);
    r.mom = struct('ok', true, ...
                   'E_sol',   E_SOL_st_MoM, ...
                   'mean_E',  mean_E_MoM, ...
                   'mean_I',  mean_I_MoM, ...
                   'alpha',   alpha_MoM, ...
                   'mean_E0', mean_E0_MoM);
    params.radius = effective_radius;
catch err
    warning('eval_TM_results:MoM', 'MoM failed: %s', err.message);
    r.mom = struct('ok', false);
end

% --- Polarizability theory ---
[Pvec, E_SOL_st_POL, alpha_Pol] = RotatingArray_2D_TM(params, params.E_inc_z, 0);
r.pol = struct('Pvec', Pvec, 'E_sol', E_SOL_st_POL, 'alpha', alpha_Pol);

% --- Mie (analytical) - only valid for single scatterer at origin with plane-wave excitation ---
if (params.is_plane_wave == 1 && isscalar(params.sca_x) && params.sca_x == 0 && params.sca_y == 0)
    [E_SOL_st_mie, mean_E_mie] = Mie_Series_TM(params, params.E_inc_z);
    r.mie = struct('E_sol', E_SOL_st_mie, 'mean_E', mean(mean_E_mie));
    fprintf('done MIE\n');
end

% --- Filaments ---
[E_SOL_st_fil, mean_E_fil, mean_I_fil, alpha_Fil, field_current_mat] = filaments_TM_multiple(params, 0);
r.filaments = struct('E_sol', E_SOL_st_fil, ...
                     'mean_E', mean_E_fil, ...
                     'mean_I', mean_I_fil, ...
                     'alpha',  alpha_Fil, ...
                     'field_current_mat', {field_current_mat});

r.params = params;

end
