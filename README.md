# Numerical Fluids in MATLAB

MATLAB solvers for numerical fluid dynamics, built from scratch to study discretization, solution methods, and validation. All projects use finite-difference methods and are checked against exact or analytical solutions where available.

## Projects

| Project | Focus | Key Result |
|---|---|---|
| [01 Potential Flow Over a Cylinder](01_potential_flow_cylinder) | 2D Laplace equation on a stretched polar mesh, sparse direct solve | Stream function matches analytical solution to within 5×10⁻⁴; surface Cp converges toward exact with mesh refinement |
| [02 Advection Schemes](02_advection_schemes) | Trapezoidal and RK4 time integration with CD2 and compact (CMPT4) spatial schemes | Dispersion and dissipation errors compared across grids and time steps |
| [03 Quasi-1D Duct Flow](03_quasi1d_duct_flow) | Newton solver for frictional compressible flow, forward and inverse problems | Parametric study of inlet Mach number vs. exit pressure |

## Requirements

MATLAB R2020a or newer (scripts use local functions and graphics). No additional toolboxes needed.

## Notes

Each project folder has its own README with the problem description, method, and results.
