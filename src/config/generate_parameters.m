function params = generate_parameters(params)
% generate_parameters  Fill the params struct with the derived simulation
% constants, the X/Y grids, the media properties and the filament defaults.
%
% Inputs (required fields on params):
%   params.lambda                  wavelength [m]
%   params.len, params.wid         physical extent of the X/Y grid [m]
%   params.len_n, params.wid_n     number of grid samples
%
% Outputs (fields added/overwritten on params):
%   x, y, X, Y, x_mom, y_mom, e0, mu0, c, f, k0, omega,
%   mr_in, mr_out, er_in, er_out, n_in, n_out,
%   R_in, R_out, N_filaments, N_testpoints, max_m

params.y = linspace(0, params.len, params.len_n) - params.len/2;
params.x = linspace(0, params.wid, params.wid_n) - params.wid/2;

[params.X, params.Y] = meshgrid(params.x, params.y);

params.y_mom = linspace(0, params.len, params.len_n) - params.len/2;
params.x_mom = linspace(0, params.wid, params.wid_n) - params.wid/2;

% Physical constants
params.e0  = 8.8541878128e-12;
params.mu0 = 1.25663706212e-6;

% Media defaults; caller may override before this call
if ~isfield(params, 'mr_in'),  params.mr_in  = 1;    end
if ~isfield(params, 'mr_out'), params.mr_out = 1;    end
if ~isfield(params, 'er_in'),  params.er_in  = 11.4; end
if ~isfield(params, 'er_out'), params.er_out = 1;    end

params.c     = 1/sqrt(params.e0*params.mu0);
params.f     = params.c/params.lambda;
params.k0    = 2*pi/params.lambda;
params.omega = params.f*2*pi;

params.n_in  = sqrt(params.mr_in *params.er_in);
params.n_out = sqrt(params.mr_out*params.er_out);

% Filament method defaults
params.R_out        = 1.2;
params.R_in         = 0.8;
params.N_filaments  = 30;
params.N_testpoints = 60;
params.max_m        = 50;

end
