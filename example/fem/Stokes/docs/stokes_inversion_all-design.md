# Unified Stokes inversion script

## Understanding

- Add `stokes_inversion_all.m` while preserving the three existing scripts.
- Select behavior by calling `stokes_inversion_all('slope')`,
  `stokes_inversion_all('figure')`, or
  `stokes_inversion_all('rectangle')`.
- Keep one shared copy of the forward, adjoint, finite-difference, Hessian, and update code.
- Restrict mode branches to geometry, reference Robin parameter, and plotting behavior.
- Preserve the existing `scheme = 4/5/6` experiments and numerical defaults.

## Assumptions

- `slope` reproduces `stokes_inversion.m`.
- `figure` reproduces the plotting/export intent of `stokes_inversion2.m`.
- `rectangle` reproduces the horizontal geometry and constant `m0` of
  `stokes_inversion3.m`.
- Existing helper functions and MATLAB path setup remain unchanged.
- The original scripts remain available as regression references.

## Decision log

1. Use one script with a mode switch instead of three copied algorithm branches.
   This keeps numerical fixes synchronized.
2. Use the mode names `slope`, `figure`, and `rectangle` as requested.
3. Keep numerical parameters near the top of the function for direct MATLAB use.
4. Keep plotting differences behind boolean configuration flags.
5. Do not delete or redirect the original scripts.
6. Expose the mode as a function argument; calling without an argument defaults
   to `slope`.

## Validation

- Check each supported mode is accepted and invalid modes fail early.
- Compare mode-specific constants and plotting calls with the corresponding
  original script.
- Run a MATLAB-compatible syntax check if MATLAB or Octave is available.
