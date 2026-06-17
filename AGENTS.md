# Repository Guidelines

## Project Structure & Module Organization

iFEM is a MATLAB finite element package. Core routines are grouped by role: `mesh/` for mesh generation and refinement, `fem/` for finite element assembly and error routines, `solver/` for direct and multigrid solvers, `dof/` for degree-of-freedom maps, `transfer/` for interpolation and transfer operators, and `tool/` for plotting and utilities. PDE data files live in `data/`. Runnable demonstrations and regression-style checks are under `example/`, with focused experiments in `research/` and ad hoc diagnostics in `debug/`. Documentation sources are in `docs/`; generated or legacy MATLAB documentation is in `ifemdoc/`.

## Build, Test, and Development Commands

There is no compiled build step for the MATLAB code. From MATLAB, run:

```matlab
cd /path/to/ifem
setpath
```

`setpath` adds repository subfolders to the MATLAB path and removes `docs/` from it. Use examples as smoke tests, for example:

```matlab
Poissonfemrate
Poisson3femrate
Maxwell3ND0femrate
```

For documentation site work, run commands from `docs/`:

```sh
bundle install
bundle exec jekyll serve
```

This serves the site locally at the configured `/ifem/` base path.

## Coding Style & Naming Conventions

Write MATLAB functions in `.m` files using clear function names that match the filename, such as `uniformrefine3.m` or `Poisson3P2.m`. Follow the existing vectorized sparse-matrix style rather than adding slow element-by-element loops when avoidable. Use descriptive numerical names already common in the codebase: `node`, `elem`, `bdFlag`, `pde`, `option`, `soln`, `eqn`, and `info`. Keep comments concise and include MATLAB help text for new public functions.

## Testing Guidelines

Prefer small, runnable MATLAB scripts in the relevant `example/`, `debug/`, or `research/` subdirectory. Name test or rate scripts after the method being exercised, for example `PoissonRT0mfemrate.m` or `transferedgetest.m`. When changing numerical routines, run the closest 2D and, if applicable, 3D example and check convergence rates, residuals, or displayed zero-difference diagnostics.

## Commit & Pull Request Guidelines

Git history uses short, imperative summaries such as `fix update of circumcenter` and `add ivem`. Keep commits focused on one numerical method, bug fix, or documentation update. Pull requests should describe the mathematical or API change, list MATLAB examples run, note affected folders, and include screenshots only for plotting or documentation-site changes.

## Security & Configuration Tips

Do not commit machine-specific MATLAB paths, generated caches, or private datasets. Keep large benchmark data in `data/` only when it is required by checked-in examples.
