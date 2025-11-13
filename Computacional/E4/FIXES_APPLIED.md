# Fixes Applied to Match JoelCampo_Burgers.m

## Critical Parameter Fixes

### 1. C1 (Time Step Coefficient)
- **Before**: `C1 = 0.025`
- **After**: `C1 = 0.001` 
- **Impact**: This was causing a **25x difference** in time step size, leading to completely different convergence behavior and results
- **Formula**: `time_step = C1 * Re / N^2`

### 2. Algorithm Matching
- Removed overly restrictive bounds checking in convective term loop
- Removed `real()` wrapper from energy calculation to match exactly: `E_k = u_after .* conj(u_after)`
- Convective term now matches Joel's code exactly

## Important: Clear Cache

**You MUST clear the cache** to get fresh results with the corrected parameters:

```matlab
clearCache()
```

Or manually delete the `cache` directory.

## Parameter Summary

For **Re=100** (and any Re value), the correct parameters are:
- `C1 = 0.001` (NOT 0.025)
- `tolerance = 1e-05`
- `maxIt = 10^5`
- `m = 2`

## Verification

After clearing cache and running with Re=100, your DNS (N=20) results should match Joel's `case_study = 1` results.

