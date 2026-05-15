function main(target, overrides)
% MAIN  Entry point for the Polarization simulations.
%
%   main()                            run the default scenario ('single')
%   main(target)                      run scenario with built-in defaults
%   main(target, overrides)           run scenario with overrides on top of defaults
%
% Available targets:
%   'single'        Single rotating scatterer comparison (MoM / Filaments / Pol)
%   'multiple'      Vogel-array multiple scatterers comparison
%   'blockcrystals' Spectral sweep over a rotating block-crystal array
%                   (requires data/arraypoints.mat)
%
% Overrides is a struct whose fields replace same-named defaults inside the
% chosen scenario. See `help single_scatterer`, `help multiple_scatterers`,
% `help block_crystals` (or the `defaults = struct(...)` block at the top of
% each example) for the available knobs. Example:
%
%   main('single', struct('lambda', 0.8e-6, 'omega_factors', 0))

if nargin < 1, target    = 'single';   end
if nargin < 2, overrides = struct();   end

add_paths();

switch lower(target)
    case 'single'
        single_scatterer(overrides);
    case 'multiple'
        multiple_scatterers(overrides);
    case {'blockcrystals', 'block_crystals', 'block'}
        block_crystals(overrides);
    otherwise
        error('main:unknownTarget', ...
              'Unknown target ''%s''. Expected ''single'', ''multiple'', or ''blockcrystals''.', target);
end

end

function add_paths()
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'src')));
addpath(fullfile(here, 'examples'));
end
