function main(target)
% MAIN  Entry point for the Polarization simulations.
%
%   main('single')        Single rotating scatterer comparison (MoM / Filaments / Pol)
%   main('multiple')      Vogel-array multiple scatterers comparison
%   main('blockcrystals') Spectral sweep over a rotating block-crystal array
%                         (requires data/arraypoints.mat)
%
% With no argument, runs the single-scatterer example.

if nargin < 1, target = 'single'; end

add_paths();

switch lower(target)
    case 'single'
        single_scatterer();
    case 'multiple'
        multiple_scatterers();
    case {'blockcrystals', 'block_crystals', 'block'}
        block_crystals();
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
