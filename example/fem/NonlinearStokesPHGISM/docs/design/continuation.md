# Nonlinear Stokes regularization continuation

## Purpose

Add continuation examples for the manufactured-solution and physical
ice-slab problems.  The examples reduce `eps_reg` from `1e-1` to `1e-4`
while using the previous velocity as the next nonlinear initial guess.

## Design

- Add one standalone script for each problem.
- Keep the mesh fixed during continuation so `option.u0` remains compatible.
- Use `epsList = 10.^(-1:-1:-4)`.
- For the manufactured solution, rebuild the PDE data at every stage so its
  body force and traction match the current regularization.
- Record convergence, iteration counts, viscosity ranges, and MMS errors.
- Stop after the first failed nonlinear solve.
- Plot the solution from the last converged stage.

## Assumptions

- These are small local MATLAB experiments; parallel execution is unnecessary.
- The solver accepts a velocity initial guess but not a pressure initial guess.
- Existing single-run examples remain unchanged.

## Decision log

- Use two explicit scripts instead of a shared helper to keep the examples
  readable and independently configurable.
- Use a fixed mesh rather than combining mesh refinement with continuation.
- Do not propagate a failed stage into the next regularization level.
