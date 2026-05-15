# Experimental: TE polarization

These files implement (partially) the TE counterpart of the methods used in
`src/`. They are **not** placed on the default MATLAB path by `main.m` and are
not exercised by any of the entry points in `examples/`.

## Contents

- `Mie_Series_TE.m` — Mie analytic solution for TE polarization. Single
  scatterer at the origin only. Last reviewed: incomplete (e.g. line 40 is
  missing a trailing semicolon, and several output branches are stubs).
- `filaments_TE_multiple.m` — Filament method for TE polarization. Uses the
  old `(source, test, n, k0, c, OMEGA)` signature of `dyadic_green`, which
  the current `src/green/dyadic_green.m` no longer accepts.
- `eval_TE_results.m` — Was the TE counterpart to `src/runners/eval_TM_results.m`
  in workspace-script form. Still operates on `{i,j,t}` cell arrays in the
  caller workspace; needs to be ported to the new function/struct contract
  before it can be wired back in.

## Known issues, in priority order

1. `filaments_TE_multiple.m` calls `dyadic_green(source, test, n, k0, c, OMEGA)`,
   but `src/green/dyadic_green.m` is `dyadic_green(source, test, params, inside)`.
   Either add a TE-aware code path to the main `dyadic_green`, or restore a
   `dyadic_green_te.m` helper alongside `scalar_green_mom.m`.
2. The previous post-processor for TE results (`TE_alpha_presentation`) was
   referenced from the entry scripts but never existed in the repository.
3. There is no TE counterpart of `RotatingArray_2D_TM.m` — the central
   theoretical contribution is currently TM-only.

## To re-activate

Once the issues above are resolved, add this folder to the MATLAB path in
`main.m` and add a `te` target to the dispatcher.
