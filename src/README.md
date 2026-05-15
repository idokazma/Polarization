# src/

Library code grouped by role. Every `.m` file here is a function with explicit
inputs and outputs — nothing reads or writes the base workspace.

| Folder       | Role                                                                                          |
|--------------|-----------------------------------------------------------------------------------------------|
| `config/`    | Build the `params` struct: physical constants, X/Y grids, default media properties, filament defaults. |
| `geometry/`  | Scatterer layouts. Currently just `GA_generator` (Vogel / golden-angle spiral).               |
| `green/`     | 2D rotating-frame Green's functions used by every method (`scalar_green`, `dyadic_green`, `scalar_green_mom`). |
| `methods/`   | The four computational methods: Mie series, MoM, Filaments, and rotating-frame Polarizability theory. |
| `sources/`   | Incident-field generators: plane wave and line-current source.                                |
| `runners/`   | Per-configuration orchestration. `eval_TM_results(params)` calls every method and returns a results struct. |
| `plotting/`  | Comparison plots (`plot_fil_mom_pol`) and the polarizability-tensor presentation (`TM_alpha_presentation_{fields,currents}`). |

## Typical call order

```
generate_parameters(params)
        │
        ▼
generate_plane_wave(params)  or  generate_source_wave(params)
        │
        ▼
eval_TM_results(params)
        │  ├─ MoM
        │  ├─ RotatingArray_2D_TM
        │  ├─ Mie_Series_TM         (only single scatterer at origin + plane wave)
        │  └─ filaments_TM_multiple
        ▼
results{i,j,t}  ──► plot_fil_mom_pol / TM_alpha_presentation_*
```

## Conventions

- The `params` struct is the single carrier of state. Anything a function
  needs (lambda, grid, media, OMEGA, scatterer positions, source location,
  incident field) is a field on `params`.
- Functions return new outputs; they do not mutate the caller's `params`
  (MATLAB pass-by-value).
- Green's functions take `(source, test, params, inside)` where `inside`
  selects between `n_in/mu_in` and `n_out/mu_out`.
- `scalar_green_mom` is a legacy variant of `scalar_green` with an explicit
  `(n, k0, c, OMEGA)` signature, kept because `MoM.m` still calls it.

For the partial TE implementation (off the default path), see
`experimental/TE/README.md`.
