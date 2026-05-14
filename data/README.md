# data/

This folder holds large binary inputs that are **not** committed to git
(see top-level `.gitignore`: `*.mat`).

## Required files

- **`arraypoints.mat`** — needed by `examples/block_crystals.m`. Must contain a
  variable `ArrayPoints` of size `Nx2` (x,y coordinates of the scatterers in
  units of µm). Without this file, `main('blockcrystals')` aborts with a
  clear error message.

If you have the original thesis dataset, drop the `.mat` file in this folder.
Otherwise the array can be regenerated from `src/geometry/GA_generator.m` —
e.g.:

```matlab
ArrayPoints = GA_generator(100).';
save('data/arraypoints.mat', 'ArrayPoints');
```
