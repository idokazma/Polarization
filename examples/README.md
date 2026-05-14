# examples/

Each file is a self-contained scenario that pulls together pieces from `src/`
and produces the final comparison plots. They are invoked via the top-level
`main.m` dispatcher:

```matlab
>> main('single')         % examples/single_scatterer.m
>> main('multiple')       % examples/multiple_scatterers.m
>> main('blockcrystals')  % examples/block_crystals.m
```

| Example                    | What it does                                                                                            | Sweep                                    | Notes |
|----------------------------|---------------------------------------------------------------------------------------------------------|------------------------------------------|-------|
| `single_scatterer.m`       | One rotating dielectric cylinder; compares MoM, Filaments, and rotating-frame Polarizability theory.    | `shift_vec` × `OMEGA_vec` × source ring  | Source-wave excitation. |
| `multiple_scatterers.m`    | Vogel-spiral array of cylinders; same three-method comparison.                                          | `OMEGA_vec` (shift fixed at 0)           | Uses `GA_generator` for positions. |
| `block_crystals.m`         | Spectral sweep over a fixed scatterer layout; computes the relative field at three observation points.  | `wavelengths` (no rotation)              | Requires `data/arraypoints.mat` — see `data/README.md`. |

## Adding a new example

1. Create `examples/my_scenario.m` as a function (no input args).
2. Populate the `params` struct (lambda, geometry, sca_x/sca_y, source, sweep
   vectors) and call `params = generate_parameters(params)`.
3. Loop over your sweep, calling `generate_plane_wave` or `generate_source_wave`
   followed by `eval_TM_results(params)`, accumulating into a results cell array.
4. Hand the results to `plot_fil_mom_pol(results)` and/or
   `TM_alpha_presentation_fields(results)`.
5. Add a case branch to `main.m` so it can be dispatched by name.
