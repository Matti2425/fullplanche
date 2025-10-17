# Convergence Study Implementation

This document describes the new convergence study functionality added to the CFD simulation.

## Overview

The convergence study validates the numerical CFD solution against analytical solutions by running simulations with progressively refined meshes (8×8, 16×16, 32×32, 64×64 nodes).

## New Files Created

### 1. `compute_analytical_psi.m`
Computes the analytical stream function for three cases:

- **CASE 1**: Uniform channel flow
  - Equation: ψ(x, y) = v_in·y
  
- **CASE 2**: Flow around static cylinder (with doublet)
  - Equation: ψ(x, y) = v_in·y + (μ/2π)·(y-y₀)/((x-x₀)² + (y-y₀)²)
  
- **CASE 3**: Flow around rotating cylinder (with circulation)
  - Equation: ψ(x, y) = ψ_static - (Γ/2π)·ln(√((x-x₀)² + (y-y₀)²))

### 2. `convergence_study.m`
Main script that:
1. Loops over mesh sizes (N, M = 8, 16, 32, 64)
2. Runs simulations for each case
3. Computes analytical solutions
4. Calculates error metrics:
   - **L₂-norm error**: √(Σ(ψ_num - ψ_analytical)²)
   - **L∞-norm error**: max|ψ_num - ψ_analytical|
5. Analyzes convergence rates
6. Generates convergence plots
7. Saves results to `convergence_study_results.mat`

## Running the Convergence Study

### Quick Start
```matlab
convergence_study
```

This will execute the complete convergence study for all three cases and all mesh refinements.

## Output

### Console Output
- Progress information for each case and mesh size
- Convergence rates (should be close to 2.0 for 2nd-order accuracy)
- Summary table with all metrics

### Files Generated
- **convergence_study_results.mat**: Contains all convergence data
  - `mesh_sizes`: [8, 16, 32, 64]
  - `L2_error`: L₂-norm errors for each case and mesh
  - `Linf_error`: L∞-norm errors for each case and mesh
  - `dx`: Grid spacing for each mesh
  - `psi_numerical`: Numerical solutions stored for each case/mesh
  - `psi_analytical`: Analytical solutions for each case/mesh
  - `iterations`: Solver iterations needed
  - `convergence_error`: Final solver residuals

### Plots Generated
A 2×2 figure showing:
1. **L₂-norm Convergence**: Error vs grid spacing with 1st and 2nd order reference lines
2. **L∞-norm Convergence**: Maximum error vs grid spacing
3. **Solver Iterations**: Number of iterations vs mesh size
4. **Solver Residual**: Final convergence residual vs mesh size

## Expected Results

For a well-implemented 2nd-order finite-difference scheme:
- **Convergence Rate**: Should approach 2.0 as mesh refines
- **L₂-norm Error**: Should decrease by ~4× when mesh is halved (dx → dx/2)
- **L∞-norm Error**: Should follow similar trend

### Typical Output Example
```
Case | L2 Rate | Linf Rate | Expected (2nd order)
-----|---------|-----------|--------------------
1    |  1.95   |   1.98   |         2.0
2    |  1.92   |   1.95   |         2.0
3    |  1.88   |   1.91   |         2.0
```

## Accessing Results

After running convergence_study, you can load and analyze results:

```matlab
% Load convergence data
load('convergence_study_results.mat', 'convergence_data');

% Access specific data
L2_errors = convergence_data.L2_error;      % 3×4 matrix
mesh_sizes = convergence_data.mesh_sizes;   % [8, 16, 32, 64]

% Extract solutions for a specific case
psi_num_case2 = convergence_data.psi_numerical{2, 3};  % Case 2, 32×32 mesh
psi_ana_case2 = convergence_data.psi_analytical{2, 3};
```

## Modifying the Study

### Change Mesh Sizes
Edit line 11 in `convergence_study.m`:
```matlab
mesh_sizes = [8, 16, 32, 64];  % Change these values
```

### Change Cases Studied
Edit line 12:
```matlab
num_cases = 3;  % Change to 1, 2, or 3 to run specific cases only
```

### Modify Tolerances
In `readInputData.m`, change `Numerical.d` for solver tolerance.

## Technical Details

### Stream Function Relationship
The velocity field is derived from the stream function using:
- u_x = ∂ψ/∂y
- u_y = -∂ψ/∂x

### Domain
- Length (L): 2 m
- Height (H): 1 m
- Cylinder diameter (D): 0.2 m
- Cylinder center: (1 m, 0.5 m)
- Inlet velocity (v_in): 1 m/s

### Boundary Conditions
- **Inlet (x=0)**: ψ = v_in·y
- **Outlet (x=L)**: Fully developed (dψ/dx = 0)
- **Top wall (y=H)**: ψ = v_in·H
- **Bottom wall (y=0)**: ψ = 0
- **Cylinder surface**: ψ = psi_body (constant)

## Troubleshooting

### Issue: Too many iterations for large meshes
- **Solution**: Increase `Numerical.i_max` in `readInputData.m`

### Issue: Convergence not improving with mesh refinement
- **Solution**: Check boundary conditions and numerical scheme implementation

### Issue: Memory error with 64×64 mesh
- **Solution**: Run convergence study in parts for each case separately

## References

The analytical solutions are based on potential flow theory:
- Uniform flow: Standard channel flow
- Cylinder (static): Doublet superimposed on uniform flow
- Cylinder (rotating): Doublet + Rankine vortex (Magnus effect)