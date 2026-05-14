# Polarization

Simulation and analysis of polarizability in a rotating reference frame —
MATLAB code accompanying a Master's thesis under the supervision of
Prof. Ben Z. Steinberg at Tel Aviv University. The numerics support the
analytic results derived in
[Steinberg, Shamir & Boag, "Rotating Green's function for 2D problems" (PRE)](https://www.eng.tau.ac.il/~steinber/papers_ps/RotatingG_PRE_SteinbergShamirBoag.pdf).

## How to run

From MATLAB at the repository root:

```matlab
>> main('single')         % single rotating scatterer (default)
>> main('multiple')       % Vogel-array multiple scatterers
>> main('blockcrystals')  % spectral sweep (needs data/arraypoints.mat)
```

`main.m` adds `src/` and `examples/` to the path and dispatches to the
chosen example.

## Layout

```
main.m                  Dispatcher entry point
examples/               Top-level scenario scripts
src/
  config/               Build the params struct + grid + constants
  geometry/             Scatterer layouts (e.g. Vogel spiral)
  green/                2D rotating-frame Green's functions
  methods/              Computational methods (Mie, MoM, Filaments, Polarizability)
  plotting/             Comparison plots and polarizability-tensor presentation
  runners/              Per-configuration orchestration (eval_TM_results)
  sources/              Incident field generators (plane wave, line source)
experimental/TE/        Partial TE-polarization implementation (not on default path)
data/                   Large binary inputs (gitignored; see data/README.md)
docs/figures/           Figures generated for the thesis / paper
```

## Methods

| Method                        | File                                  | Scope                                                    |
|-------------------------------|---------------------------------------|----------------------------------------------------------|
| Mie series (analytical)       | `src/methods/Mie_Series_TM.m`         | Single scatterer at origin                               |
| Method of Moments (MoM)       | `src/methods/MoM.m`                   | Multiple scatterers, volumetric                          |
| Filament method               | `src/methods/filaments_TM_multiple.m` | Multiple scatterers, boundary                            |
| Rotating-frame polarizability | `src/methods/RotatingArray_2D_TM.m`   | Multiple scatterers, dipole theory (the novel contribution) |

## TE polarization

The TE counterparts of the methods are incomplete and parked under
`experimental/TE/` with their own README. They are not on the default path.
