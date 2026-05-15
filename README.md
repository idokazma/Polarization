# Polarization

**2D simulation of polarizability in a rotating reference frame.**

MATLAB code accompanying a Master's thesis by Ido Kazma at Tel Aviv University,
under the supervision of Prof. Ben Z. Steinberg. The numerics support the
analytic results derived in
[Steinberg, Shamir & Boag, "Rotating Green's function for 2D problems" (PRE)](https://www.eng.tau.ac.il/~steinber/papers_ps/RotatingG_PRE_SteinbergShamirBoag.pdf).

The repository simulates electromagnetic scattering from one or more
dielectric cylinders that rotate together with the reference frame, and
extracts the effective polarizability tensor — comparing a closed-form
"polarizability-theory" model against two independent numerical solvers
(Method of Moments and the Filament method) and, where applicable, against
the analytical Mie series.

## Requirements

- MATLAB R2018b or newer (uses `isfile`, `isscalar`, the Symbolic Math Toolbox
  for the symbolic source representations in `src/sources/`).
- The Symbolic Math Toolbox is only needed by the symbolic blocks in
  `generate_plane_wave` / `generate_source_wave`; the numerics work without it.

## Quick start

From MATLAB at the repository root:

```matlab
>> main('single')         % single rotating scatterer (default)
>> main('multiple')       % Vogel-spiral array of scatterers
>> main('blockcrystals')  % spectral sweep over a fixed block-crystal layout
```

`main.m` puts `src/` (recursively) and `examples/` on the path and dispatches
to the chosen scenario. To start from a fresh MATLAB session, just `cd` into
the repo root and call `main(...)` — no manual `addpath` needed.

The `blockcrystals` target requires `data/arraypoints.mat`, which is not
checked in. See `data/README.md` for how to provide or regenerate it.

### Overriding parameters

Each scenario exposes its defaults at the top of its file. Pass a struct as a
second argument to `main` (or directly to the scenario function) to override
any of them:

```matlab
>> main('single', struct('lambda', 0.8e-6, 'omega_factors', 0))
>> main('multiple', struct('n_scatterers', 20, 'eps_vec', 4))
>> main('blockcrystals', struct('wavelengths', (0.5:0.01:1.5)*1e-6))
```

Fields you don't list keep their default values. The full list of knobs is
the `defaults = struct(...)` block at the top of each `examples/*.m` file.

## What the code does

For each configuration (wavelength, scatterer layout, rotation rate Ω,
incident field), `src/runners/eval_TM_results.m` runs the four methods and
returns a struct of results. The example scripts sweep over Ω (and other
parameters), accumulate the per-configuration structs into a cell array, and
hand the array to the plotting and presentation functions.

The central object of study is the **effective polarizability tensor** α
linking the induced dipole moments to the local (H, E) field at each
scatterer. Under rotation, α acquires off-diagonal terms that the
polarizability-theory model predicts in closed form; the two numerical
methods (MoM and Filaments) serve as independent benchmarks.

| Method                        | File                                  | Scope                                            |
|-------------------------------|---------------------------------------|--------------------------------------------------|
| Mie series (analytical)       | `src/methods/Mie_Series_TM.m`         | Single scatterer at origin, plane-wave only      |
| Method of Moments (MoM)       | `src/methods/MoM.m`                   | Multiple scatterers, volumetric discretisation   |
| Filament method               | `src/methods/filaments_TM_multiple.m` | Multiple scatterers, boundary discretisation     |
| Rotating-frame polarizability | `src/methods/RotatingArray_2D_TM.m`   | Multiple scatterers, dipole theory (the novel contribution) |

All four share the 2D rotating-frame Green's functions in `src/green/`
(`scalar_green` for the E_z component, `dyadic_green` for the (H_x, H_y)
components).

## Repository layout

```
main.m                  Dispatcher entry point
examples/               Scenario scripts called by main.m
  single_scatterer.m
  multiple_scatterers.m
  block_crystals.m
src/
  config/               Build the params struct + grid + constants
  geometry/             Scatterer layouts (e.g. Vogel spiral)
  green/                2D rotating-frame Green's functions
  methods/              Computational methods (Mie, MoM, Filaments, Polarizability)
  plotting/             Comparison plots + polarizability-tensor presentation
  runners/              Per-configuration orchestration (eval_TM_results)
  sources/              Incident-field generators (plane wave, line source)
experimental/TE/        Partial TE-polarization implementation (off the default path)
data/                   Large binary inputs (gitignored)
docs/figures/           Figures generated for the thesis / paper
```

Each top-level subdirectory has its own README with more detail.

## The `params` struct

A single struct named `params` carries all simulation state — physical
constants, the X/Y grid, media properties, scatterer positions, source
location, the incident field, and the sweep vectors. The intended setup
pattern is:

```matlab
params.lambda = 1e-6;          % wavelength [m]
params.len    = 0.5e-6;        % grid extent
params.wid    = 0.5e-6;
params.len_n  = 200;           % grid samples
params.wid_n  = 200;
params.tm     = 1;             % TM polarization
params.te     = 0;
% ... source / scatterer / sweep config ...

params = generate_parameters(params);   % fills constants, X/Y, n_in/out, ...
params = generate_source_wave(params);  % or generate_plane_wave(params)
result = eval_TM_results(params);       % runs all four methods
```

See `src/README.md` for the full conventions.

## Adding a new scenario

1. Create `examples/my_scenario.m` as a function with no arguments.
2. Populate `params`, call `generate_parameters`, then loop over your sweep
   calling `eval_TM_results` and accumulating into a cell array.
3. Hand the cell array to `plot_fil_mom_pol` and/or
   `TM_alpha_presentation_fields`.
4. Add a case branch to `main.m` so the scenario can be dispatched by name.

`examples/README.md` walks through this in more detail.

## TE polarization

The TE counterparts of the methods (`Mie_Series_TE`, `filaments_TE_multiple`,
`eval_TE_results`) are incomplete and parked under `experimental/TE/`. They
are not on the default path and the active entry points do not call them.
See `experimental/TE/README.md` for the known issues and what's needed to
re-activate them.

## Known limitations

- The rotating-frame polarizability method is **TM only**. There is no TE
  counterpart of `RotatingArray_2D_TM`.
- The Mie series is **single-scatterer only**, used as a sanity check
  against MoM for the single-cylinder case.
- `data/arraypoints.mat` is not checked in (`*.mat` is gitignored), so the
  `blockcrystals` scenario won't run out of the box on a fresh clone.

## License

No license file is currently present. The code is research output from a
Master's thesis; if you intend to reuse it, please contact the authors.
