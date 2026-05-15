function params = generate_source_wave(params)
% generate_source_wave  Set the incident E/H fields on params to those of a
% line-current source located at (params.source_loc_x, params.source_loc_y).

if (params.tm)
    params.E_inc_z = scalar_green([params.source_loc_x; params.source_loc_y], ...
                                  [params.X(:), params.Y(:)]', params, 0);
    params.E_inc_z = reshape(params.E_inc_z, length(params.X), []);
    params.E_inc_x = 0*params.E_inc_z;
    params.E_inc_y = 0*params.E_inc_z;

    [params.H_inc_x, params.H_inc_y] = dyadic_green( ...
        [params.source_loc_x; params.source_loc_y], ...
        [params.X(:), params.Y(:)]', params, 0);
    params.H_inc_x = reshape(params.H_inc_x, length(params.X), []);
    params.H_inc_y = reshape(params.H_inc_y, length(params.X), []);
    params.H_inc_z = 0*params.H_inc_y;
end

if (params.te)
    params.H_inc_z = scalar_green([params.source_loc_x; params.source_loc_y], ...
                                  [params.X(:), params.Y(:)]', params, 0);
    params.H_inc_z = reshape(params.H_inc_z, length(params.X), []);
    params.H_inc_x = 0*params.H_inc_z;
    params.H_inc_y = 0*params.H_inc_z;

    [params.E_inc_x, params.E_inc_y] = dyadic_green( ...
        [params.source_loc_x; params.source_loc_y], ...
        [params.X(:), params.Y(:)]', params, 0);
    params.E_inc_x = reshape(1i/(params.omega*(params.e0*params.er_out)) * params.E_inc_x, length(params.X), []);
    params.E_inc_y = reshape(1i/(params.omega*(params.e0*params.er_out)) * params.E_inc_y, length(params.X), []);
    params.E_inc_z = 0*params.E_inc_y;
end

end
