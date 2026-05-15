function params = generate_plane_wave(params)
% generate_plane_wave  Set the incident E/H fields on params to those of a
% plane wave propagating along +x (TM mode) or with selectable direction (TE).

if (params.tm)
    params.E_inc_z = exp(-1i*params.k0*params.X);
    params.E_inc_x = 0*params.E_inc_z;
    params.E_inc_y = 0*params.E_inc_z;

    params.H_inc_y = 1i/(params.omega*(params.mu0*params.mr_out)) * ...
                     (-1i*params.k0*exp(-1i*params.k0*params.X));
    params.H_inc_x = 0*params.H_inc_y;
    params.H_inc_z = 0*params.H_inc_y;

    syms xsym ysym zsym k0sym
    params.Esym = [0, 0, exp(-1i*k0sym*xsym)];
    params.Xsym = [xsym ysym zsym];
    params.Hsym = 1/(1i*params.omega*(params.mu0*params.mr_out)) * curl(params.Esym, params.Xsym);
end

if (params.te)
    if (params.plane_wave_direction == 1)
        params.H_inc_z = exp(-1i*params.k0*params.X);
        params.H_inc_x = 0*params.H_inc_z;
        params.H_inc_y = 0*params.H_inc_z;

        params.E_inc_y = -1i/(params.omega*(params.e0*params.er_out)) * ...
                         (-1i*params.k0*exp(-1i*params.k0*params.X));
        params.E_inc_x = 0*params.E_inc_y;
        params.E_inc_z = 0*params.E_inc_y;

        syms xsym ysym zsym k0sym
        params.Hsym = [0, 0, exp(-1i*k0sym*xsym)];
        params.Xsym = [xsym ysym zsym];
        params.Esym = -1/(1i*params.omega*(params.e0*params.er_out)) * curl(params.Hsym, params.Xsym);
    else
        params.H_inc_z = exp(-1i*params.k0*params.Y);
        params.H_inc_x = 0*params.H_inc_z;
        params.H_inc_y = 0*params.H_inc_z;

        params.E_inc_x = 1i/(params.omega*(params.e0*params.er_out)) * ...
                         (-1i*params.k0*exp(-1i*params.k0*params.Y));
        params.E_inc_y = 0*params.E_inc_x;
        params.E_inc_z = 0*params.E_inc_x;

        syms xsym ysym zsym k0sym
        params.Hsym = [0, 0, exp(-1i*k0sym*ysym)];
        params.Xsym = [xsym ysym zsym];
        params.Esym = -1/(1i*params.omega*(params.e0*params.er_out)) * curl(params.Hsym, params.Xsym);
    end
end

end
