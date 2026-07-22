# Finite-difference Newton inversion design

## Understanding

- Rebuild `NSFDInversion.m` as a strict comparison with
  `NonlinearStokesAdjInvSlabBed.m`.
- Match its geometry, mesh, PDE coefficients, observations, true parameter,
  initial parameter, objective, regularization, LM damping, line search, and
  stopping criteria.
- Compute both the reduced-objective gradient and full Hessian by finite
  differences.
- Do not use adjoint equations, tangent equations, or a `J'*J`
  Gauss--Newton approximation.
- Keep the implementation intended for the current small inverse space
  (`Nm=4`).

## Assumptions

- The only intended algorithmic difference from the adjoint script is the
  derivative implementation.
- Accuracy and comparison clarity are more important than runtime.
- A second-order centered stencil is used for every derivative.
- Hessian symmetry is enforced numerically before the Newton solve.

## Decision log

1. Use the adjoint script as the structural and configuration reference.
2. Differentiate the complete reduced objective, including regularization.
3. Use centered differences for the gradient, Hessian diagonal, and Hessian
   mixed entries.
4. Solve `(H + lambda*I)*step = -gradient`.
5. Retain the adjoint script's decreasing-objective backtracking and LM
   update.
6. Report gradient norm, Hessian symmetry error, and the number of nonlinear
   forward evaluations used by each inverse iteration.
7. Print Hessian eigenvalues, the number of negative eigenvalues, and a
    spectral condition estimate to expose weakly observable directions.

## Validation

- MATLAB `checkcode` must parse the modified file.
- The finite-difference Hessian must be symmetric to roundoff after explicit
  symmetrization.
- A shortened smoke test must complete at least one inverse iteration.
- Compare the configuration values against the adjoint script.
