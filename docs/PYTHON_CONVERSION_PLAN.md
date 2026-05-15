# Python Conversion Plan

A staged plan for porting the MATLAB `Polarization` codebase to Python while
preserving the numerical behaviour of the four TM methods (Mie, MoM,
Filaments, rotating-frame Polarizability) and the example-driven CLI.

## 1. Scope and non-goals

In scope
- All TM code currently exercised by `main('single')`, `main('multiple')`,
  and `main('blockcrystals')`.
- The four solvers in `src/methods/`, the Green's functions in `src/green/`,
  the parameter / source / geometry helpers, and the comparison plots
  (`plot_fil_mom_pol`, `TM_alpha_presentation_fields`).
- Round-trip parity tests against MATLAB reference outputs for the three
  scenarios.

Out of scope (initial port)
- `experimental/TE/` — the TE methods are documented as incomplete; port
  them only after the TM pipeline is reproduced.
- `TM_alpha_presentation_currents.m` if it turns out to be unused on the
  default paths (verify before porting).
- Re-deriving the symbolic plane-wave / line-source field expressions from
  scratch; the analytical closed forms are already inlined alongside the
  symbolic blocks in `generate_plane_wave.m` / `generate_source_wave.m`, so
  we will port those directly and drop the symbolic toolbox dependency.

## 2. Target stack

| Concern              | MATLAB                       | Python replacement                         |
|----------------------|------------------------------|--------------------------------------------|
| Arrays / linear alg. | builtin, `\`, `linsolve`     | `numpy`, `numpy.linalg.solve`              |
| Sparse solve         | `sparse`, `\`                | `scipy.sparse`, `scipy.sparse.linalg.spsolve` |
| Bessel / Hankel      | `besselj`, `besselh(m,1,z)`  | `scipy.special.jv`, `scipy.special.hankel1` |
| Polar conversion     | `cart2pol`                   | `numpy.arctan2`, `numpy.hypot`             |
| `.mat` I/O           | `load`                       | `scipy.io.loadmat`                         |
| Plotting             | `figure`, `imagesc`, `scatter` | `matplotlib.pyplot`                       |
| Symbolic (plane / source wave) | Symbolic Math Toolbox | none — inline the closed-form expressions |
| `params` struct      | MATLAB struct                | `dataclass` (`SimParams`) with explicit fields |
| `cell` arrays of results | nested cells           | nested `list` (or `numpy.empty(shape, dtype=object)`) |
| CLI dispatch         | `main.m` switch              | `argparse` in `polarization/__main__.py`   |

Python version: 3.11+. Pin via `pyproject.toml` with `numpy`, `scipy`,
`matplotlib` as runtime deps; `pytest`, `pytest-regressions` for tests.

## 3. Proposed package layout

Mirror the MATLAB layout one-to-one so reviewers can diff side by side:

```
polarization/
  __init__.py
  __main__.py                   # argparse dispatcher (replaces main.m)
  config/
    parameters.py               # generate_parameters
    defaults.py                 # merge_defaults equivalent (just dict merge)
  geometry/
    vogel.py                    # GA_generator
  green/
    scalar.py                   # scalar_green, scalar_green_mom
    dyadic.py                   # dyadic_green
  methods/
    mie_tm.py                   # Mie_Series_TM
    mom.py                      # MoM
    filaments_tm.py             # filaments_TM_multiple
    rotating_array_tm.py        # RotatingArray_2D_TM
  sources/
    plane_wave.py               # generate_plane_wave
    line_source.py              # generate_source_wave
  runners/
    eval_tm.py                  # eval_TM_results
  plotting/
    compare.py                  # plot_fil_mom_pol
    alpha_fields.py             # TM_alpha_presentation_fields
examples/
  single_scatterer.py
  multiple_scatterers.py
  block_crystals.py
tests/
  test_green.py
  test_geometry.py
  test_mie.py
  test_mom.py
  test_filaments.py
  test_rotating_array.py
  test_examples_smoke.py
  fixtures/                     # small reference outputs from MATLAB
pyproject.toml
README.md (Python section appended)
```

## 4. The `params` struct → `SimParams` dataclass

A single `@dataclass` (or `pydantic` model — vote for dataclass to avoid the
dep) collects every field the MATLAB struct accumulates. Group fields by
provenance, mark derived ones as `field(init=False)`, fill them in a
`finalize()` method that is the analogue of `generate_parameters`. Replace
in-place `params.x = ...` with `dataclasses.replace` or explicit mutation —
the MATLAB code is pass-by-value, but we don't need to mimic that since
each scenario rebuilds `params` each loop iteration anyway.

Important: keep numeric values as `numpy.complex128` where the MATLAB code
relies on complex arithmetic (Green's functions, fields). Use explicit
`dtype=complex` on `np.zeros` allocations to avoid silent upcasts.

## 5. Indexing and shape conventions

Every loop and slice needs an indexing audit:
- `1:N` → `range(N)`; `linspace(0, L, N)` is identical except for endpoint
  inclusion — MATLAB and NumPy both include the endpoint, so safe.
- `theta = linspace(0, 2*pi, N+1); theta = theta(1:end-1)` →
  `np.linspace(0, 2*np.pi, N, endpoint=False)`.
- `[X, Y] = meshgrid(x, y)` → `X, Y = np.meshgrid(x, y)` (default
  `indexing='xy'` matches MATLAB's row-major-ish meshgrid layout).
- `A(:)` → `A.ravel(order='F')` (MATLAB is column-major). For grid arrays
  produced by `meshgrid` where downstream code only treats them as a flat
  list, we can use `order='C'` — but only after confirming with a parity
  test, since `MoM.m` builds `scatterers_x, scatterers_y` by stacking
  `X_mom(loc)` and indexes into the same linear order later.
- Source / test point arrays passed to Green's functions are MATLAB
  2×N matrices (`[x; y]`). Either keep that shape in Python (`(2, N)` ndarray)
  or transpose to `(N, 2)` — pick one convention package-wide. Recommend
  `(N, 2)` because NumPy is row-major, and update the Green's function
  signatures accordingly.
- `besselh(m, 1, z)` is the Hankel function of the first kind:
  `scipy.special.hankel1(m, z)`. Watch the call sites where `m` is a
  numpy array (Mie series uses `m=0:max_m`) — `hankel1` and `jv` broadcast
  natively.

## 6. Module-by-module port order

Each step is mergeable on its own and has a parity test against MATLAB
output captured ahead of time (see §8).

1. **`config` + `geometry` + `green`** — pure functions, no plotting, no
   solvers. Highest test-to-effort ratio; everything else depends on them.
   - `generate_parameters` (48 lines) → `SimParams.finalize()`.
   - `merge_defaults` (19 lines) → trivial dict overlay (or just use
     dataclass `replace`).
   - `GA_generator` (20 lines) → `vogel.golden_angle_array(N)`.
   - `scalar_green`, `scalar_green_mom`, `dyadic_green` — small (≤40 lines
     each), all vectorised, port verbatim. Add a docstring that lists the
     `(x, y)` broadcasting expectations because the MATLAB code reshapes
     1×N row vectors implicitly.

2. **`sources/`** — `generate_plane_wave` and `generate_source_wave` (~50
   lines each). Drop the `syms`/`curl` blocks: the closed-form `Esym/Hsym`
   they produce is the same expression as the analytic line above them.
   Audit each example to confirm `params.Esym` / `params.Hsym` are not
   consumed downstream — `filaments_TM_multiple.m:88` only references them
   inside a `if (false)` block, and `MoM.m:112` declares the symbols but
   never uses the result. Both can be dropped on the Python side.

3. **`methods/mie_tm.py`** — single-scatterer analytical reference (94
   lines, vectorised). Easiest solver to validate (it's a closed form).
   Use it as the first end-to-end correctness anchor.

4. **`methods/rotating_array_tm.py`** — 203 lines, dense system over
   `N_scatterers`, dominated by the `GRTM` Green's helper which duplicates
   `scalar_green` with the wrong `omega*mu0` prefactor. Decide explicitly
   whether to call `green.scalar_green` here or keep the local `GRTM`
   variant; the MATLAB code keeps both, and tests must lock in whichever
   we choose. Recommend: keep `GRTM` inline, document the prefactor
   convention, do not try to "clean up" until parity is locked in.

5. **`methods/mom.py`** — 374 lines, the hairiest module. Significant
   commented-out scaffolding; port only the live code path
   (`MoM(params)` → `(SOL, mean_E, mean_I, effective_radius, alpha_mom,
   mean_E0)`). The `calc_full_sol == 1` branch is gated and currently
   unreachable from the example scripts; port it but mark it
   `@pytest.mark.skip(reason="no example exercises this branch")`.

6. **`methods/filaments_tm.py`** — 591 lines. Largest module; expect this
   to be the longest leg of the port. Split into helpers
   (`_build_filaments`, `_assemble_system`, `_evaluate_fields`) as we go,
   driven by parity tests, not by refactoring instinct.

7. **`runners/eval_tm.py`** — straightforward dispatcher; mirrors the
   MATLAB `eval_TM_results` exactly (same struct-of-method-outputs shape).

8. **`plotting/`** — `plot_fil_mom_pol` (68 lines) is direct.
   `TM_alpha_presentation_fields` (656 lines) contains the alpha-tensor
   regression and the final figures; port last and validate the figures by
   eye against MATLAB-rendered PNGs.

9. **`examples/` + `__main__.py`** — scenario scripts and CLI dispatcher.
   Each example is a function that builds a `SimParams`, sweeps, and
   feeds results to the plotting layer. CLI: `python -m polarization
   single --lambda 0.8e-6 --omega-factors 0`. Use `argparse` subparsers,
   one per scenario, with `--override key=value` for arbitrary fields to
   match the MATLAB `overrides` struct ergonomics.

## 7. CLI / overrides ergonomics

MATLAB callers can do `main('single', struct('lambda', 0.8e-6))`. The
closest Python equivalent that does not require code edits is JSON or
`key=value` pairs on the command line:

```
python -m polarization single \
    --override lambda=8e-7 \
    --override omega_factors='[-3e-5, 0, 3e-5]'
```

Values are parsed as JSON so lists / floats / strings round-trip cleanly.
We can additionally accept `--config path/to/overrides.json` for the
case where users want to keep parameter sets under version control.

## 8. Testing and parity strategy

The goal is to keep the port honest, not to invent a new test suite from
zero. Approach:

1. **Capture MATLAB reference outputs once.** For each method, write a
   small MATLAB harness (or use the existing examples with a fixed seed)
   to save `(params, output)` pairs to `.mat` for representative inputs.
   Check the `.mat` files into `tests/fixtures/` (override `.gitignore`
   for that subtree).
   - Per-method numeric tests: e.g. `tests/fixtures/mie_single.mat`,
     `tests/fixtures/mom_15_scatterers.mat`,
     `tests/fixtures/filaments_15_scatterers.mat`.
   - Per-scenario end-to-end snapshot: the cell array returned by
     `single_scatterer` collapsed to a stack of E_sol / alpha arrays.

2. **`pytest` parity tests** load the fixture via `scipy.io.loadmat`,
   re-run the Python implementation on the same inputs, and assert
   `np.allclose(out_py, out_matlab, rtol=1e-10, atol=1e-12)`. Use
   per-method tolerances; Filaments and MoM may need looser `rtol`
   (1e-8) because they involve solving a linear system.

3. **Smoke tests** run each scenario on a tiny grid (`len_n=11`, two
   omegas, one source) and assert finite outputs, expected shapes, and
   that the plotting layer does not raise. Mark with
   `@pytest.mark.slow` for the full-fidelity runs.

4. **Plot validation** is by eye, not assertion. Generate PNGs in CI
   under `docs/figures/python/` for each scenario at fixed parameters
   and diff against the MATLAB-rendered PNGs in `docs/figures/`.

## 9. Risks and known gotchas

- **Column-major flattening.** MATLAB's `A(:)` is column-major and the
  MoM code relies on it implicitly when it stacks `X_mom(loc)` and later
  reshapes. Every `.ravel()` / `.flatten()` call in the port must specify
  `order='F'` until proven otherwise.
- **`besselh` branch / sign.** SciPy's `hankel1(0, z)` agrees with
  MATLAB's `besselh(0, 1, z)` for `Re(z) > 0`. Confirm with a unit test
  at a handful of complex `z` values before trusting bulk results.
- **`linsolve` vs sparse solve.** `RotatingArray_2D_TM` builds `MM`, then
  `sparse(MM)`, then `MM\rhs`. The matrix is dense (the sparsification
  is cosmetic). Use `numpy.linalg.solve` and ignore the sparse wrapper.
- **`scalar_green` vs `scalar_green_mom`.** Two near-identical functions
  with different signatures. MATLAB keeps both. Keep both in Python too —
  fixing this is a separate refactor after parity is locked in.
- **Symbolic toolbox.** Already discussed; the symbolic blocks are dead
  code on the active paths. Removing them removes a heavy dependency.
- **`data/arraypoints.mat`.** Required by `blockcrystals`, not in the
  repo. Document the fallback (skip the scenario with a clear error
  message) until the file is provided.
- **MATLAB 1-based slicing.** Every `i = 1:N` loop and every
  `arr(start_ind:(start_ind+L-1))` slice in `MoM.m:226-234` needs a
  manual translation. Single most common porting bug; budget a parity
  test for it.
- **Cell array `{...}` semantics.** In `eval_TM_results.m:45`, the
  filaments output is stored as `{field_current_mat}` (a 1-cell wrapper)
  and unpacked later with `cell2mat`. The Python port should not
  replicate the wrapping — return the array directly and update callers.

## 10. Milestones

| # | Milestone                                                  | Exit criterion                              |
|---|------------------------------------------------------------|---------------------------------------------|
| 1 | Package skeleton + `green/` + `geometry/` + `config/`       | Unit tests pass at 1e-12 against fixtures.  |
| 2 | `sources/` + Mie + parity test                              | `Mie_Series_TM` parity at 1e-10.            |
| 3 | `RotatingArray_2D_TM` + parity test                         | Single-scatterer alpha matches Mie to 1e-3. |
| 4 | `MoM` + parity test                                         | 15-scatterer MoM matches MATLAB at 1e-8.    |
| 5 | `filaments_TM_multiple` + parity test                       | 15-scatterer Filaments matches at 1e-8.     |
| 6 | `eval_TM_results` + `__main__.py` + smoke tests            | All three scenarios run end to end.         |
| 7 | `plot_fil_mom_pol` + `TM_alpha_presentation_fields`         | Side-by-side figures match MATLAB output.   |
| 8 | README rewrite, MATLAB code marked legacy                   | Docs explain both stacks; CI green.         |

Estimated effort: ~3-5 days of focused work for milestones 1-4, ~2-3
days for milestones 5-6 (Filaments is the bottleneck), ~1-2 days for
plotting and docs. Total: ~7-10 working days for a complete TM port at
parity, excluding any back-and-forth on the alpha-tensor presentation
plots.

## 11. Open questions

These should be resolved before starting the implementation:

1. **MATLAB co-existence.** Keep `main.m` and the `.m` files in the repo
   indefinitely as a reference, or remove them once parity is achieved?
   Recommend keeping them under `legacy/matlab/` for one release cycle.
2. **TE port.** Defer until TM is done, or port in parallel? Recommend
   defer — the TE methods are flagged as incomplete and we should fix
   them in MATLAB or Python, not both.
3. **Plotting target.** Static matplotlib figures, or interactive
   (`plotly`, `bokeh`)? Recommend matplotlib for parity with the thesis
   figures; revisit later.
4. **Distribution.** Is this published as a PyPI package, or does it
   live only in this repo? Affects whether we need a `src/` layout and a
   `console_scripts` entry point.
