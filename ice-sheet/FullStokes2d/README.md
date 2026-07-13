# FullStokes2d

This directory contains two-dimensional nonlinear full-Stokes ice-flow examples and inversion scripts.  The main solver is a P2--P1 Taylor--Hood discretization with Glen viscosity, periodic side boundaries, top traction, and basal sliding/no-slip variants.

## Layout

```text
FullStokes2d/
  NonlinearStokesP2P1.m        % main 2-D full-Stokes solver
  FullStokesAdjInvSlab.m       % nondimensional slab-bed adjoint inversion
  FullStokesAdjInvSin.m        % nondimensional sinusoidal-bed adjoint inversion
  FullStokesAdjInvSlabISM.m    % PHGISM/ISM internal-unit slab inversion
  HOM/                         % ISMIP-HOM B/D forward benchmarks
  test/                        % MMS tests, diagnostic scripts, older inversions
  output/                      % generated EPS figures
```

## Quick Start

From MATLAB:

```matlab
cd('/path/to/ifem')
setpath
seticepath
cd ice-sheet/FullStokes2d
```

Recommended smoke tests:

```matlab
run('test/NSConverRate.m')
FullStokesAdjInvSlab
FullStokesAdjInvSlabISM
```

ISMIP-HOM forward examples live in `HOM/`:

```matlab
cd HOM
HOM_B
HOM_D
```

## Units

Most legacy test scripts are nondimensional.  `FullStokesAdjInvSlabISM.m` is the PHGISM/ISM internal-unit version:

- coordinates: km;
- velocity: m/year;
- basal friction coefficient: Pa yr m^{-1};
- pressure output: PHGISM pressure DOF `p/1e5`.

The ISM script scales the coefficients passed to `NonlinearStokesP2P1.m` with `LEN_SCALING = 1e3`, `EQU_SCALING = 1e-8`, and `PRES_SCALING = 1e5`.

## Main Solver

### `NonlinearStokesP2P1.m`

```matlab
[soln,eqn,info] = NonlinearStokesP2P1(node,elem,bdFlag,pde,option);
```

Important `pde` fields:

```matlab
pde.A
pde.n
pde.beta
pde.m
pde.f
pde.rho
pde.gravity
pde.g_N
```

Important `option` fields:

```matlab
option.periodic
option.periodic_x
option.eps_reg
option.maxIt
option.tol
option.residual_tol
option.damping
option.quadorder
option.assemble_tangent
option.pressure_constraint
```

Pressure constraint modes:

```matlab
'auto'
'mean-zero'
'none'
```

## Inversion Scripts

### `FullStokesAdjInvSlab.m`

Adjoint beta inversion on a sloping slab bed.  This is the compact nondimensional benchmark for the boundary-integral objective, adjoint gradient, and matrix-free Gauss--Newton update.

### `FullStokesAdjInvSin.m`

Same inversion framework as `FullStokesAdjInvSlab.m`, but with an ISMIP-HOM-B-like sinusoidal bed geometry.

### `FullStokesAdjInvSlabISM.m`

PHGISM/ISM internal-unit slab-bed beta inversion.  Use this script when comparing units and scaling with the PHGISM ice-sheet code.

Generated EPS figures are written under:

```text
output/<scriptName>/
```

## ISMIP-HOM Benchmarks

### `HOM/HOM_B.m`

2-D flowline version of ISMIP-HOM experiment B.  The script compares this solver only with the official full-Stokes curve in `HOM_B_official_fs_curve.csv`.

### `HOM/HOM_D.m`

2-D flowline version of ISMIP-HOM experiment D.  The script compares this solver only with the official full-Stokes curve in `HOM_D_official_fs_curve.csv`.

### `HOM/HOM_B_L5Fields.m`, `HOM/HOM_D_L5Fields.m`

Field plots for the `L = 5 km` benchmark cases.

## Tests and Diagnostics

The `test/` directory contains MMS data, convergence tests, derivative checks, finite-difference inversion prototypes, and older diagnostic scripts:

```text
test/NSConverRate.m
test/NSEpsContinuation.m
test/NSDerivativeComparison.m
test/NSFDInversion.m
test/NSAdjInvTikhonov.m
test/NSRegression.m
```

After changing the solver, run at least:

```matlab
run('test/NSConverRate.m')
run('test/NSEpsContinuation.m')
run('test/NSDerivativeComparison.m')
```

Check convergence flags, nonlinear residuals, MMS error decay, divergence error, and adjoint derivative errors.
