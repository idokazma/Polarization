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

## Overriding defaults

Every scenario function takes an optional `overrides` struct. The defaults
live in a single `defaults = struct(...)` block at the top of each file;
`merge_defaults(defaults, overrides)` shallowly overlays the caller's values
on top. So:

```matlab
>> main('single', struct('lambda', 0.8e-6))
>> single_scatterer(struct('omega_factors', 0, 'n_sources', 1))
```

Fields not present in `overrides` keep their default values. To see what
knobs a scenario exposes, open the file and read the `defaults` block — it
is the canonical list.

## Adding a new example

1. Create `examples/my_scenario.m` as `function my_scenario(overrides)`.
2. Build a `defaults = struct(...)` with every tweakable value, then
   `cfg = merge_defaults(defaults, overrides);`.
3. Populate the `params` struct from `cfg`, call
   `params = generate_parameters(params)`, then loop over your sweep calling
   `eval_TM_results(params)` and accumulating into a results cell array.
4. Hand the results to `plot_fil_mom_pol(results)` and/or
   `TM_alpha_presentation_fields(results)`.
5. Add a case branch to `main.m` so it can be dispatched by name.
