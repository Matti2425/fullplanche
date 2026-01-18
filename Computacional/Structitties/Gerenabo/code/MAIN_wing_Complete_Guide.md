# MAIN_wing.m - Complete Study Guide

**Comprehensive explanation of the complete wing structural model combining beam and shell elements.**

---

## 1. Overview: What Does This Code Do?

This script performs **Finite Element Analysis (FEA)** of a **complete aircraft wing structure** by combining the beam and shell elements from Parts I and II.

### Purpose
1. **Static Analysis**: Calculate deformations under various loads (bending, torsion, pressure + gravity)
2. **Modal Analysis**: Determine natural frequencies and vibration modes
3. **Stress Analysis**: Compute Von Mises stresses in shell components
4. **Model Order Reduction (MOR)**: Optional reduced-order solution using modal superposition

### Key Innovation
**Hybrid Modeling**: Combines different element types in a single structure:
- **Beams** (Timoshenko): Stringers (longitudinal stiffeners)
- **Shells** (Flat 4-node quad): Wing box, ribs, skin

---

## 2. Wing Structure Components

### Physical Wing Structure

A typical wing box consists of:
1. **Wing Box (Shells)**: Upper/lower skins + front/rear spars (main load-bearing box)
2. **Stringers (Beams)**: Longitudinal stiffeners attached to skins (resist bending)
3. **Ribs (Shells)**: Transverse plates that maintain airfoil shape
4. **Skin (Shells)**: Outer aerodynamic surface (can be modeled separately or as part of wing box)

```
         ┌─────────────────────┐
         │      Skin (SK)      │
         │  ┌───────────────┐  │
         │  │  Ribs (RB)    │  │
         │  │  |||||||||||  │  │ ← Stringers (ST)
         │  │  Wing Box(WB) │  │
         │  └───────────────┘  │
         └─────────────────────┘
              Root → Tip
```

---

## 3. Code Structure Walkthrough

### Part I: Solution Controls (Lines 7-17)

**Three key switches**:
```matlab
problem_type = 3;  % 1: Bending, 2: Torsion, 3: Pressure + Weight
solution_type = 2; % 1: Direct Method, 2: Model Order Reduction
Nm = 10;           % Number of modes to compute
```

**Component flags** (enables/disables parts):
```matlab
WB = 1; % Wingbox (shells)
ST = 1; % Stringers (beams)
RB = 1; % Ribs (shells)
SK = 0; % Skin (shells) - usually 0 if wing box includes skin
```

**Why toggle components?**
- Study contribution of each component to stiffness
- Computational efficiency (skip unwanted parts)
- Compare simplified vs. detailed models

---

### Part II: Data & Preprocessing (Lines 19-115)

**Material Properties** (Lines 20-24):
```matlab
E = 75e9;       % Young's Modulus (Pa) - Aluminum
v_poisson = 0.3;
rho = 2200;     % Density (kg/m³)
```

**Loads** (Lines 34-88):
```matlab
% Problem Type 1: Bending (point forces at tip)
Fe_wb = [Fa_shear indPointA 3; Fb_shear indPointB 3];

% Problem Type 2: Torsion (couple forces creating torque)
Fe_wb = [Fa_torque indPointA 3; Fb_torque indPointB 3];

% Problem Type 3: Pressure + Gravity (realistic flight condition)
Pe_sk = ComputePressureState(...);  % Aerodynamic pressure on skin
Be = [-g0 * ones(N_nodes,1), node_indices, 3*ones(N_nodes,1)];  % Gravity
```

**Element Properties**:

**Wing Box (Lines 92-95)**:
```matlab
h_vector_wb = [17.5e-3; 15e-3; 6e-3];  % Thickness varies by zone
% Different materials can be assigned via Tm_wb
```

**Stringers (Lines 98-105)**:
```matlab
r_st = 5e-3;  % Circular cross-section radius
A_st = pi*r_st^2;          % Area
I_st = (pi/4)*r_st^4;      % Moment of inertia
J_st = (pi/2)*r_st^4;      % Torsional constant
ky_vector_st = 5/6;        % Shear correction factor
```

**Ribs & Skin (Lines 108-114)**:
```matlab
h_vector_rb = 5e-3;  % Rib thickness
h_vector_sk = 2e-3;  % Skin thickness (thinner than wing box)
```

---

### Part III: Solver - Matrix Assembly (Lines 117-139)

**Key Innovation**: Separate computation functions called for each component type:

```matlab
% Wing Box (Shells) - Uses MAIN_shell.m logic
[K_wb, M_wb, f_wb, Bb_wb, Bmn_wb, Bmt_wb, Bs_wb, R_wb] = ...
    ComputeKM_shells(xn, N_nodes, size(Tn_wb,1), Tn_wb, Tm_wb, ...
                     h_vector_wb, E_vector_wb, v_vector_wb, rho_vector_wb, ...
                     Fe_wb, Be, Pe_wb);

% Stringers (Beams) - Uses MAIN_beam.m logic
[K_st, M_st, f_st] = ...
    ComputeKM_beams(xn, N_nodes, size(Tn_st,1), Tn_st, Tm_st, ...
                    E_vector_st, G_vector_st, rho_vector_st, j_orientations, ...
                    A_vector_st, Iy_vector_st, Iz_vector_st, J_vector_st, ...
                    ky_vector_st, kz_vector_st, kt_vector_st, ...
                    Fe_st, Be, Qe_st);

% Ribs (Shells)
[K_rb, M_rb, f_rb, ...] = ComputeKM_shells(...);

% Skin (Shells)
[K_sk, M_sk, f_sk, ...] = ComputeKM_shells(...);
```

**Why separate functions?**
- Modular code (easier to debug/modify)
- Reusable (same shell function for wing box, ribs, and skin)
- Each component can have different properties

**Output matrices**:
- `K_*`: Stiffness matrix for component
- `M_*`: Mass matrix for component
- `f_*`: Force vector for component
- `Bb_*, Bmn_*, ...`: Strain-displacement matrices (saved for stress calculation)
- `R_*`: Rotation matrices (for stress transformation)

---

### Part IV: Global Assembly (Lines 141-150)

**Component Superposition**:
```matlab
K = WB*K_wb + ST*K_st + RB*K_rb + SK*K_sk;
M = WB*M_wb + ST*M_st + RB*M_rb + SK*M_sk;
f = WB*f_wb + ST*f_st + RB*f_rb + SK*f_sk;
```

**Key Insight**: All components share the **same nodal mesh** (`xn`), so their matrices have the same dimensions and can be added directly!

**Active DOFs** (Line 149):
```matlab
active_dofs = find(diag(K) ~= 0 | diag(M) ~= 0);
```
This filters out DOFs that aren't connected to any elements (important when components are toggled off).

---

### Part V: Modal Analysis (Lines 152-172)

**Standard Eigenvalue Solution**:
```matlab
[V, D] = eigs(K(free_dofs, free_dofs), M(free_dofs, free_dofs), Nm_eff, 'sm');
[lambda_sorted, idx] = sort(diag(D));
freq = sqrt(lambda_sorted) / (2*pi);  % Convert to Hz
```

**Mass Normalization** (Lines 168-171):
```matlab
for k = 1:Nm_eff
    norm_factor = sqrt(Phi(:,k)' * M * Phi(:,k));
    Phi(:,k) = Phi(:,k) / norm_factor;  % Φ^T * M * Φ = I
end
```

**Why normalize?**
- Standard practice in structural dynamics
- Makes modes orthonormal with respect to mass matrix
- Simplifies modal superposition formulas

---

### Part VI: Static Solution (Lines 174-197)

**Method 1: Direct Solution** (Lines 179-182):
```matlab
f_eff = f(free_dofs) - K(free_dofs, prescribed_dofs) * u(prescribed_dofs);
u(free_dofs) = K(free_dofs, free_dofs) \ f_eff;
```

Standard FEM: Solve $\mathbf{K}_{ff}\mathbf{u}_f = \mathbf{f}_f$

**Method 2: Model Order Reduction (MOR)** (Lines 184-196):
```matlab
% Superposition: u ≈ Σ Φ_i * α_i
for k = Im_eff
    modal_force = Phi(:,k)' * f;        % Generalized force
    stiffness_val = lambda_sorted(k);   % Modal stiffness (ω_k²)
    alpha = modal_force / stiffness_val; % Modal amplitude
    u_reduced = u_reduced + Phi(:,k) * alpha;
end
```

**MOR Theory**:
$$\mathbf{u} \approx \sum_{i=1}^{m} \phi_i \alpha_i, \quad \alpha_i = \frac{\phi_i^T \mathbf{f}}{\omega_i^2}$$

**Advantages of MOR**:
- Much faster for large systems (only uses $m$ modes instead of all DOFs)
- Physical insight (see which modes dominate response)
- Works well when first few modes capture most of the deformation

**When to use Direct vs. MOR?**
- **Direct**: More accurate, required for high-frequency/localized loads
- **MOR**: Faster, good for global deformation with dominant low-frequency modes

---

### Part VII: Stress Calculation (Lines 199-205)

**Von Mises Stress** (for shell components):
```matlab
if WB, sigVM_wb = Compute_stress_shells(..., u, ..., Bb_wb, Bmn_wb, Bmt_wb, Bs_wb, R_wb); end
```

**How it works**:
1. Extract displacements `u` for element nodes
2. Use saved `B` matrices to compute strains: $\boldsymbol{\varepsilon} = \mathbf{B}\mathbf{u}_e$
3. Compute stresses: $\boldsymbol{\sigma} = \mathbf{C}\boldsymbol{\varepsilon}$
4. Calculate Von Mises: $\sigma_{VM} = \sqrt{\sigma_x^2 - \sigma_x\sigma_y + \sigma_y^2 + 3\sigma_{xy}^2}$

**Note**: Stresses computed separately for membrane, bending, and shear, then combined.

---

## 4. Post-Processing (Lines 207-299)

### 2D Spanwise Plots (Lines 210-267)

**Extract twist and deflection**:
```matlab
for i = 1:n_points
    uz_front = u(6*(node_front-1) + 3);
    uz_rear  = u(6*(node_rear-1) + 3);
    
    % Twist angle from differential deflection
    theta_x(i) = (uz_rear - uz_front) / dist_y;
    
    % Deflection at shear center
    u_z(i) = uz_front + theta_x(i) * (yc - y1);
end
```

**Physical Meaning**:
- Front and rear spar nodes at each spanwise station
- Difference in Z-displacement → twist angle
- Interpolate to shear center (where lift acts) for meaningful deflection

### 3D Visualization (Lines 269-293)

**Deformed shape with stress contours**:
```matlab
plotDeformed('wing', xn, Tn_wb, u, scale, sigVM_wb);
```

Shows:
- Undeformed geometry (wireframe)
- Deformed geometry (colored by stress)
- Separate plots for wing box, ribs, skin

### Mode Shapes (Lines 295-297)

```matlab
plotModes('wing', Phi, freq, 1:min(9, Nm_eff));
```

Visualizes first 9 vibration modes with frequencies.

---

## 5. Key Variable Reference

| Variable | Dimensions | Description |
|:---------|:-----------|:------------|
| `xn` | $N_{nodes} \times 3$ | **Shared** nodal coordinates for all components |
| `Tn_wb` | $N_{wb} \times 4$ | Wing box element connectivity (shells) |
| `Tn_st` | $N_{st} \times 2$ | Stringer element connectivity (beams) |
| `K_wb`, `K_st`, etc. | $6N \times 6N$ | Component stiffness matrices |
| `K`, `M`, `f` | $6N \times 6N$ | **Global** assembled system |
| `WB`, `ST`, `RB`, `SK` | Scalar (0 or 1) | Component toggle flags |
| `Phi` | $6N \times N_m$ | Mode shapes (columns = modes) |
| `Im` | Vector | Indices of modes to use for MOR |

---

## 6. Exam Q&A

### Q1: Why can we add K_wb + K_st + K_rb + K_sk directly?
**A**: All components share the **same global node numbering** in `xn`. Each component matrix is assembled into the same $6N \times 6N$ global space. When a component doesn't connect to a node, those entries are zero, so addition is valid.

### Q2: What's the point of the component flags (WB, ST, RB, SK)?
**A**: Study individual contributions:
- Set `WB=1, ST=0, RB=0, SK=0` → See wing box-only stiffness
- Compare deflections with/without stringers → Quantify stringer effectiveness
- Debugging: Isolate issues to specific components

### Q3: How does MOR work and when should I use it?
**A**: **Modal superposition**: Approximate solution as $\mathbf{u} \approx \sum_{i=1}^{m} \phi_i \alpha_i$

**When to use**:
- Large systems where direct solve is slow
- Loads dominated by low-frequency response
- First $m=9$ modes usually capture \u003e95% of global deformation

**When NOT to use**:
- High-frequency/localized loads (need more modes)
- Impact/transient analysis (direct better)

### Q4: Why are stringers modeled as beams but skin as shells?
**A**: **Geometry matches element type**:
- **Stringers**: Long, slender → 1D beam ideal
- **Skin/Wing box**: Surface → 2D shell ideal

Using the right element type gives accuracy with fewer DOFs.

### Q5: What's the difference between the 3 problem types?
**A**:
| Type | Load | Purpose |
|:-----|:-----|:--------|
| 1: Bending | Point force at tip | Test vertical stiffness |
| 2: Torsion | Torque (couple) at tip | Test torsional stiffness |
| 3: Pressure + Weight | Distributed pressure + gravity | Realistic flight condition |

### Q6: Why compute stress only for shells, not beams?
**A**: The code focuses on shell stresses because:
- Shells have complex 2D stress states (membrane + bending)
- Stringer stresses can be computed from beam stress formulas if needed
- Von Mises stress is most relevant for thin-walled structures (shells)

### Q7: Can I change stringer cross-section from circular to rectangular?
**A**: Yes! Modify lines 98-101:
```matlab
% Rectangular: width w, height h
A_st = w * h;
Iy_st = (w * h^3) / 12;
Iz_st = (h * w^3) / 12;
J_st = ???;  % Complex for rectangle, use torsion constant formula
```

### Q8: What determines the number of modes needed for MOR?
**A**: Trial and error:
- Run with `Im = 1:3` (first 3 modes)
- Compare to direct solution
- If error \u003e tolerance, increase to `Im = 1:6`, etc.
- Usually first 6-9 modes sufficient for global deformation

---

## 7. Important Theoretical Concepts

### Hybrid Modeling

**Principle**: Different elements for different structural components based on geometry and behavior.

| Component | Element Type | Reason |
|:----------|:------------|:-------|
| Stringers | Beam (2-node) | Slender, resist bending/axial |
| Wing Box | Shell (4-node) | Surface, resist membrane+bending |
| Ribs | Shell (4-node) | Maintain shape, transfer shear |
| Skin | Shell (4-node) | Aerodynamic surface, carry pressure |

### Model Order Reduction (MOR)

**Full System**:
$$\mathbf{K}\mathbf{u} = \mathbf{f} \quad \text{(solve for } N \text{ DOFs)}$$

**Reduced System**:
$$\mathbf{u} \approx \mathbf{\Phi}_m \boldsymbol{\alpha}, \quad \boldsymbol{\alpha} = (\mathbf{\Phi}_m^T\mathbf{K}\mathbf{\Phi}_m)^{-1}(\mathbf{\Phi}_m^T\mathbf{f})$$

where $\mathbf{\Phi}_m$ contains first $m$ mode shapes.

**Computational Cost**:
- Direct: $O(N^3)$ for $N$ DOFs
- MOR: $O(m^3)$ for $m \ll N$ modes → **Much faster!**

### Component Contribution Analysis

**Stiffness contribution**:
```matlab
% Wing box only
K1 = K_wb;  deflection_wb_only = ...

% Wing box + stringers
K2 = K_wb + K_st;  deflection_with_stringers = ...

% Stringer effectiveness
improvement = (deflection_wb_only - deflection_with_stringers) / deflection_wb_only * 100
```

---

## 8. Comparison: Beam vs Shell vs Wing

| Aspect | MAIN_beam.m | MAIN_shell.m | MAIN_wing.m |
|:-------|:------------|:-------------|:------------|
| **Element Types** | Beam only | Shell only | **Both** |
| **Structure** | Single beam | Single shell | Complete wing assembly |
| **Components** | 1 | 1 | 4 (WB, ST, RB, SK) |
| **DOFs/node** | 6 | 6 | 6 (shared mesh) |
| **Matrices** | K, M | K, M | K_wb+K_st+K_rb+K_sk |
| **Solution** | Direct only | Direct only | Direct **or** MOR |
| **Stress** | Not computed | Von Mises | Von Mises (shells only) |
| **Load Types** | Point force/moment | Point/pressure | Point/pressure/gravity |

---

## 9. Summary Flowchart

```
1. Set controls (problem_type, solution_type, component flags)
2. Load shared mesh (xn) and connectivity for each component
3. Define BCs (clamped root) and loads based on problem_type
4. FOR each component IF flag=1:
   ├─ Wing Box (WB):  Call ComputeKM_shells → K_wb, M_wb, f_wb
   ├─ Stringers (ST): Call ComputeKM_beams  → K_st, M_st, f_st
   ├─ Ribs (RB):      Call ComputeKM_shells → K_rb, M_rb, f_rb
   └─ Skin (SK):      Call ComputeKM_shells → K_sk, M_sk, f_sk
5. Assemble global: K = Σ K_i, M = Σ M_i, f = Σ f_i
6. Modal analysis: Solve K*Φ = ω²*M*Φ (get modes for MOR)
7. Static solution:
   ├─ Direct: Solve K*u = f directly
   └─ MOR: u ≈ Σ Φ_i * (Φ_i'*f / ω_i²)
8. Compute stresses in shells: σ_VM = f(u, B_matrices)
9. Post-process:
   ├─ Spanwise plots (deflection, twist)
   ├─ 3D deformed shape with stress contours
   └─ Mode shape visualization
```

---

## 10. Advanced Topics

### Load Distribution

**Pressure loading** (Problem Type 3):
```matlab
Pe_sk = ComputePressureState(xn, n_u, n_l, b, c, p_inf, alpha);
```

This function likely:
1. Computes aerodynamic pressure distribution (based on airfoil theory)
2. Integrates pressure over shell element area
3. Returns equivalent nodal forces

**Gravity loading**:
```matlab
Be = [-g0 * ones(N_nodes,1), node_indices, 3*ones(N_nodes,1)];
```

Applied as body force (force per unit volume) distributed to all nodes.

### Stress Recovery

**Why save B matrices?**
```matlab
[K_wb, M_wb, f_wb, Bb_wb, Bmn_wb, Bmt_wb, Bs_wb, R_wb] = ComputeKM_shells(...);
```

During matrix assembly, we compute strain-displacement matrices. Instead of recomputing during post-processing, we **save** them:
- `Bb`: Bending B-matrix
- `Bmn`, `Bmt`: Membrane B-matrices (normal, tangential)
- `Bs`: Shear B-matrix
- `R`: Rotation matrix

Then in `Compute_stress_shells`, we use: $\boldsymbol{\varepsilon} = \mathbf{B}\mathbf{u}_e \rightarrow \boldsymbol{\sigma} = \mathbf{C}\boldsymbol{\varepsilon}$

### Modal Truncation Error

**How many modes is enough?**

Check **modal participation factors**:
```matlab
participation = abs(Phi' * f);  % How much force projects onto each mode
cumulative = cumsum(participation.^2) / sum(participation.^2);
% Use modes until cumulative > 0.95 (95% of response)
```

---

## 11. Practical Tips for Exam

### Toggling Components
```matlab
% Which component contributes most to stiffness?
% Run 4 times with different flags:
WB=1, ST=0, RB=0, SK=0  → deflection1
WB=0, ST=1, RB=0, SK=0  → deflection2
WB=0, ST=0, RB=1, SK=0  → deflection3
WB=0, ST=0, RB=0, SK=1  → deflection4
% Largest deflection = least stiff component
```

### Choosing Solution Method
```matlab
% For exam efficiency:
solution_type = 2;  % MOR is faster
Im = 1:6;           % Usually sufficient

% For accuracy (checking answer):
solution_type = 1;  % Direct method
```

### Interpreting Results
- **u_z plot**: Should be smooth, max at tip (cantilever)
- **theta_x plot**: Torsion → linear increase; Bending → near zero
- **Frequencies**: First mode typically 1-10 Hz for wing (bending); higher modes involve torsion/ribs
