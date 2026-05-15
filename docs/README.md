# docs/

Non-code artifacts associated with the project.

## figures/

Figures produced for the thesis / paper. Currently:

- `epsFig.eps` — large EPS figure (not referenced by any `.m` file). Kept here
  rather than at the repository root so it doesn't clutter the entry point.

Note: `.eps`, `.svg`, `.fig` are gitignored by default at the repo root. If you
add new figures here that should be tracked, either commit them with
`git add -f` or narrow the `.gitignore` rules.
