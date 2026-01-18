# MAIN_shell.m - Complete Study Guide

**Comprehensive explanation of the Finite Element Analysis implementation for 4-node Quadrilateral Flat Shell elements.**

---

## 1. Overview: What Does This Code Do?

This script performs **Finite Element Analysis (FEA)** of shell structures using **4-node Quadrilateral Flat Shell Elements** as presented in **CE06 Flat shells.pdf**.

### Purpose
1. **Static Analysis**: Calculates deformations under applied loads (bending or torsion)
2. **Modal Analysis**: Determines natural frequencies and vibration mode shapes

### Key Equation (PDF page 23)
The flat shell element combines **three components**:
$$\mathbf{K}_e = \mathbf{K}_m + \mathbf{K}_b + \mathbf{K}_s$$

Where:
- $\mathbf{K}_m$: **Membrane** stiffness (in-plane loading)
- $\mathbf{K}_b$: **Bending** stiffness (plate behavior)
- $\mathbf{K}_s$: **Shear** stiffness (transverse shear deformation)

Plus drilling rotation stabilization:
- $\mathbf{K}_r$: Penalty stiffness for rotation about normal (prevents singularity)

---

## 2. Code Structure Walkthrough

### Part I: Initialization & Data (Lines 9-30)

**Lines 10-11**: Problem setup
```matlab
problem_type = 'bending'; % bending or torsion
Nm = 6;                   % Number of modes to extract
```

**Lines 14-16**: Material properties
```matlab
E = 75e9;         % Young's Modulus (Pa)
v_poisson = 0.3;  % Poisson's ratio
rho = 2200;       % Density (kg/m³)
```

**Lines 23-30**: Geometry
- Shell thickness varies: `h1`, `h2`, `h3` (different zones)
- Geometric parameters for load application

### Part II: Preprocessing (Lines 32-86)

**Lines 36**: Load mesh
```matlab
load('shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');
% xn: Nodal coordinates [N_nodes × 3]
% Tn: Element connectivity [N_elements × 4] (4-node quads)
% Tm: Material IDs (maps to thickness)
```

**Lines 45-55**: Boundary conditions (Clamped root)
```matlab
% Fix all 6 DOFs at root nodes
Up = [0, indRoot(i), dof];  % Fully clamp root edge
```

**Lines 64-73**: Load application
- **Bending**: Forces in Z-direction
- **Torsion**: Forces creating a couple about beam axis

---

## 3. Theoretical Foundation → Code Implementation

### Shell Element Theory (PDF pages 23-25)

A flat shell element has **5 local DOFs per node**:
- $u, v, w$: Translations in local $x', y', z'$
- $\theta_x, \theta_y$: Rotations about local $x', y'$ axes
- (Rotation about $z'$ is handled separately via drilling DOF)

But in **global 3D**, each node has **6 DOFs** → requires rotation transformation

---

### A. Local Coordinate System (Lines 105-144)

**Code**:
```matlab
% Element nodes 1-2-3-4 (counterclockwise)
x_e = xn(Tn(e,:), :); % 4×3 coordinates

% Compute local axes
v13 = x_e(3,:) - x_e(1,:);
v24 = x_e(4,:) - x_e(2,:);
Normal = cross(v13, v24); 
k_dir = (Normal / norm(Normal))';  % Local z' (normal to element)

% Local x' direction
v_bisect = (x_e(2,:) + x_e(3,:) - x_e(4,:) - x_e(1,:)) / 2;
i_dir = (v_bisect / norm(v_bisect))';
j_dir = cross(k_dir, i_dir);  % Local y'

% Rotation matrix
Rp = [i_dir, j_dir, k_dir];  % 3×3 rotation matrix
```

**Creating Element Rotation Matrix (20×24)**:
```matlab
% Map 5 local DOFs to 6 global DOFs per node
R_node = zeros(5,6);
R_node(1:3, 1:3) = Rp';           % Rotate translations
R_node(4:5, 4:6) = Rp(:,1:2)';    % Rotate rotations (only θ_x, θ_y)

% Assemble for all 4 nodes
R_ele = zeros(20,24);  % 4 nodes: 5 local × 4 = 20 rows
                       %          6 global × 4 = 24 columns
for n = 1:4
    R_ele( (n-1)*5+(1:5), (n-1)*6+(1:6) ) = R_node;
end
```

**Why 5 local DOFs?** The local formulation doesn't include rotation about the element normal ($\theta_{z'}$) - this "drilling" DOF is added later as a penalty.

---

### B. Shape Functions (Lines 96-97, 202-210)

**PDF Theory (page 24)**:

Bilinear shape functions for 4-node quadrilateral:
$$N_i(\xi,\eta) = \frac{1}{4}(1 + a_i\xi)(1 + b_i\eta)$$

where:
| Node $i$ | $a_i$ | $b_i$ |
|:---------|:------|:------|
| 1 | -1 | -1 |
| 2 | +1 | -1 |
| 3 | +1 | +1 |
| 4 | -1 | +1 |

**Derivatives in natural coordinates**:
$$N_{,\xi}^i = \frac{1}{4}a_i(1 + b_i\eta), \quad N_{,\eta}^i = \frac{1}{4}b_i(1 + a_i\xi)$$

**Code**:
```matlab
a_coef = [-1 1 1 -1];
b_coef = [-1 -1 1 1];

% At Gauss point (xi, nu):
val = (1 + a_coef(n)*xi)*(1 + b_coef(n)*nu)/4;
dval_dxi = a_coef(n)*(1 + b_coef(n)*nu)/4;
dval_dnu = b_coef(n)*(1 + a_coef(n)*xi)/4;
```

**Jacobian and Physical Derivatives**:
```matlab
J = sum( dN_dxi(:,n) * x_e(n,:) * [i_dir, j_dir] );  % 2×2 Jacobian
dN_phys = J \ dN_dxi;  % Transform to physical coordinates
```

---

### C. Shear Stiffness $\mathbf{K}_s$ (Lines 148-179) ⚠️ **CRITICAL!**

**PDF Theory**: 
$$\mathbf{K}_s^e = \mathbf{R}_e^T \int_{S_e} \mathbf{B}_s^{\prime e T} \mathbf{C}_s^{\prime e} \mathbf{B}_s^{\prime e} \, dS \, \mathbf{R}_e$$

**PDF Warning (page 25)**:
> "To improve the performance of the quadrilateral shell element, it is recommended to **sub-integrate** (with just **1 Gauss point** instead of 4) the shear components."

**Code**:
```matlab
% CENTER POINT EVALUATION (xi=0, nu=0) - Reduced integration!
J_center = zeros(2,2);
for n=1:4
    dNi_dxi = [a_coef(n); b_coef(n)]/4;
    J_center = J_center + dNi_dxi * x_e(n,:) * [i_dir, j_dir];
end
invJ = inv(J_center);

% Build shear B-matrix
Bs = zeros(2, 20);  % 2 shear strains, 20 DOFs
for n=1:4
    dN_phys = invJ * [a_coef(n); b_coef(n)]/4;
    N_val = 1/4;  % Shape function at center
    
    col_idx = (n-1)*5 + (1:5);
    %              u v  w           θ_x    θ_y
    Bs(1, col_idx) = [0 0 dN_phys(1)  0       N_val];   % γ_xz = ∂w/∂x + θ_y
    Bs(2, col_idx) = [0 0 dN_phys(2) -N_val   0    ];   % γ_yz = ∂w/∂y - θ_x
end

% Shear constitutive matrix
Cs = eye(2) * (5/6) * E*h / (2*(1+v));  % Shear modulus with correction factor

% Assemble (multiply by area ≈ 4*detJ)
Ks = (4 * detJ) * (R_ele' * Bs' * Cs * Bs * R_ele);
```

**Physical Meaning**:
- Shear strain relates out-of-plane displacement derivative to rotation
- **1-point integration** at element center prevents **shear locking** (same issue as beams!)
- The factor $5/6$ is the shear correction factor

---

### D. Membrane Stiffness $\mathbf{K}_m$ (Lines 186-187, 216-225, 242)

**PDF Theory**:
$$\mathbf{K}_m^e = \mathbf{R}_e^T \int_{S_e} \mathbf{B}_m^{\prime e T} \mathbf{C}_m^{\prime e} \mathbf{B}_m^{\prime e} \, dS \, \mathbf{R}_e$$

where $\mathbf{C}_m$ is the plane stress constitutive matrix:
$$\mathbf{C}_m = \frac{Eh}{1-\nu^2} \begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2} \end{bmatrix}$$

**Code**:
```matlab
% Constitutive matrix (plane stress)
C_membrane = E*h/(1-v^2) * [1 v 0; v 1 0; 0 0 (1-v)/2];

% At each Gauss point (4 total):
for k = 1:4
    % Build membrane B-matrix
    Bm = zeros(3, 20);
    for n=1:4
        col = (n-1)*5 + (1:5);
        % Membrane strains: [ε_x; ε_y; γ_xy]
        Bm(1, col) = [dN_phys(1,n) 0 0 0 0];         % ε_x = ∂u/∂x
        Bm(2, col) = [0 dN_phys(2,n) 0 0 0];         % ε_y = ∂v/∂y
        Bm(3, col) = [dN_phys(2,n) dN_phys(1,n) 0 0 0];  % γ_xy = ∂u/∂y + ∂v/∂x
    end
    
    dV = w4(k) * detJ_k;  % Gauss weight × Jacobian
    Km = Km + dV * (R_ele' * Bm' * C_membrane * Bm * R_ele);
end
```

**Physical Meaning**:
- Handles in-plane forces and strains
- Uses **4 Gauss points** (2×2 quadrature) for exact integration

---

### E. Bending Stiffness $\mathbf{K}_b$ (Lines 189-190, 227-234, 243)

**PDF Theory**:
$$\mathbf{K}_b^e = \mathbf{R}_e^T \int_{S_e} \mathbf{B}_b^{\prime e T} \mathbf{C}_b^{\prime e} \mathbf{B}_b^{\prime e} \, dS \, \mathbf{R}_e$$

where:
$$\mathbf{C}_b = \frac{Eh^3}{12(1-\nu^2)} \begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2} \end{bmatrix}$$

**Code**:
```matlab
% Constitutive matrix (bending/plate theory)
C_bending = E*h^3/(12*(1-v^2)) * [1 v 0; v 1 0; 0 0 (1-v)/2];

% At each Gauss point:
for n=1:4
    col = (n-1)*5 + (1:5);
    % Curvatures: [κ_x; κ_y; κ_xy]
    Bb(1, col) = [0 0 0  0             dN_phys(1,n)];   % κ_x = ∂θ_y/∂x
    Bb(2, col) = [0 0 0 -dN_phys(2,n)  0           ];   % κ_y = -∂θ_x/∂y
    Bb(3, col) = [0 0 0 -dN_phys(1,n)  dN_phys(2,n)];   % κ_xy = ∂θ_y/∂y - ∂θ_x/∂x
end

Kb = Kb + dV * (R_ele' * Bb' * C_bending * Bb * R_ele);
```

**Physical Meaning**:
- Handles out-of-plane bending moments
- Curvature is derivative of rotation (like Timoshenko theory)
- Note the $h^3$ factor - bending stiffness grows cubically with thickness!

---

### F. Mass Matrix $\mathbf{M}$ (Lines 192, 236-244)

**PDF Theory**:
$$\mathbf{M}^e = \mathbf{R}_e^T \int_{S_e} \mathbf{N}_e^T \boldsymbol{\rho}'^e \mathbf{N}_e \, dS \, \mathbf{R}_e$$

**Code**:
```matlab
% Inertia matrix (5×5 for local DOFs)
rho_mat_base = rho*h*diag([1 1 1 h^2/12 h^2/12]);
%                          ↑ ↑ ↑    ↑       ↑
%                          u v w   θ_x    θ_y

% At each Gauss point:
N_mat = zeros(5, 20);
for n=1:4
    col = (n-1)*5 + (1:5);
    N_mat(:, col) = N(n) * eye(5);  % Shape function interpolation
end

Me = Me + dV * (R_ele' * N_mat' * rho_mat_base * N_mat * R_ele);
```

**Physical Meaning**:
- Translational inertia: $\rho h$ (mass per area)
- Rotational inertia: $\rho h^3/12$ (same as moment of inertia for plate)
- Uses **4 Gauss points** for exact integration

---

### G. Drilling Rotation Stabilization $\mathbf{K}_r$ (Lines 260-307) ⚠️ **IMPORTANT!**

**PDF Warning (pages 27-28)**:
> "Another issue appearing in 3D configurations is the problem with **coplanar nodes**. The model considers only 5 degrees of freedom in the local reference frame, while in the global 3D setting, we will have 6. This means for each coplanar node, we will end up having an additional equation that is either 0 or linearly dependent with another, **making the global stiffness matrix singular**."

**PDF Solution**:
> "An additional stiffness component associated to the $\theta_{z'}$ degree of freedom is added:
> $$K_{\theta_{z'}} = \frac{Eh \cdot S_e}{small \ factor}$$"

**Code**:
```matlab
% Compute element normals
for e=1:N_elements
    v13 = xn(Tn(e,3),:) - xn(Tn(e,1),:);
    v24 = xn(Tn(e,4),:) - xn(Tn(e,2),:);
    n_ele = cross(v13, v24);
    n_ele = n_ele / norm(n_ele);  % Unit normal
    
    % Average normals at nodes (to detect flat regions)
    for i=1:4
        Normal_Node(Tn(e,i), :) = Normal_Node(Tn(e,i), :) + n_ele;
    end
end

% Add penalty stiffness
Kr = sparse(Ndof, Ndof);
alpha_coeff = 1e-3 * E * h1;  % Small coefficient

for e=1:N_elements
    % ...calculate n_ele, Area...
    for i=1:4
        node = Tn(e,i);
        n_node = Normal_Node(node, :) / norm(Normal_Node(node, :));
        
        % If element is flat (normal aligns with average):
        alignment = dot(n_ele, n_node);
        if abs(alignment) > 0.9
            % Add penalty to drilling rotation about n_ele
            K_rot_penalty = (alpha_coeff * Area) * (n_ele' * n_ele);
            
            idx_rot = 6*(node-1) + (4:6);  % Global rotation DOFs
            Kr(idx_rot, idx_rot) = Kr(idx_rot, idx_rot) + K_rot_penalty;
        end
    end
end

K = K + Kr;  % Add to global stiffness
```

**Physical Meaning**:
- Without this, the matrix is singular (6th DOF has no stiffness)
- Adds small artificial stiffness to rotation about element normal
- Doesn't affect results (it's orthogonal to physical rotations)
- Only applied to flat/coplanar regions

---

## 4. Solution Methods

### Part IV: Static Analysis (Lines 316-341)

**Code**:
```matlab
% Build force vector
F_global = zeros(Ndof, 1);
for i = 1:size(Fe,1)
    node = Fe(i,2);
    dof = Fe(i,3);
    val = Fe(i,1);
    idx = 6*(node-1) + dof;
    F_global(idx) = F_global(idx) + val;
end

% Partitioning
prescribed_dofs = 6*(Up(:,2)-1) + Up(:,3);
free_dofs = setdiff(1:Ndof, prescribed_dofs);

% Solve (accounting for prescribed displacements)
F_eff = F_global(free_dofs) - K(free_dofs, prescribed_dofs) * U(prescribed_dofs);
U(free_dofs) = K(free_dofs, free_dofs) \ F_eff;
```

### Part V: Modal Analysis (Lines 343-357)

**Code**:
```matlab
n_modes_calc = min(Nm + 4, length(free_dofs));
[V_red, D_red] = eigs(K(free_dofs, free_dofs), M(free_dofs, free_dofs), n_modes_calc, 'sm');

[frequencies, sort_idx] = sort(sqrt(diag(D_red))/(2*pi));
V_red = V_red(:, sort_idx);

% Reconstruct full eigenvectors
Phi = zeros(Ndof, n_modes_calc);
Phi(free_dofs, :) = V_red;
```

---

## 5. Key Variable Reference

| Variable | Dimensions | Description |
|:---------|:-----------|:------------|
| `xn` | $N_{nodes} \times 3$ | Nodal coordinates $(x, y, z)$ |
| `Tn` | $N_{elem} \times 4$ | Element connectivity [node1, node2, node3, node4] |
| `K` | $6N \times 6N$ | Global stiffness matrix |
| `M` | $6N \times 6N$ | Global mass matrix |
| `R_ele` | $20 \times 24$ | Element rotation matrix (5 local DOFs × 4 nodes → 6 global DOFs × 4 nodes) |
| `Bs` | $2 \times 20$ | Shear strain-displacement matrix |
| `Bm` | $3 \times 20$ | Membrane strain-displacement matrix |
| `Bb` | $3 \times 20$ | Bending curvature-displacement matrix |

---

## 6. Exam Q&A

### Q1: Why does the shell element have 5 LOCAL DOFs but 6 GLOBAL DOFs?
**A**: In the local frame aligned with the shell:
- 3 translations: $u, v, w$
- 2 rotations: $\theta_x, \theta_y$ (bending rotations)
- Rotation about the normal ($\theta_z$) is the **drilling DOF** - it's not part of classical shell theory

In 3D global coordinates, all 6 DOFs exist, so we need the drilling stabilization to prevent singularity.

### Q2: Why is the shear stiffness integrated with 1 Gauss point?
**A**: **Reduced integration** to prevent **shear locking**:
- Full integration (4 points) makes thin shells too stiff
- 1-point integration (at element center) gives accurate results
- Same principle as Timoshenko beams! (PDF page 25 explicitly recommends this)

### Q3: How many Gauss points are used for each component?
**A**:
- **Shear** ($K_s$): **1 point** (center) - reduced integration
- **Membrane** ($K_m$): **4 points** (2×2 Gauss) - full integration
- **Bending** ($K_b$): **4 points** (2×2 Gauss) - full integration
- **Mass** ($M$): **4 points** (2×2 Gauss) - full integration

### Q4: What is the constitutive matrix difference between membrane and bending?
**A**:
- **Membrane**: $C_m = \frac{Eh}{1-\nu^2} [...]$ - **proportional to $h$**
- **Bending**: $C_b = \frac{Eh^3}{12(1-\nu^2)} [...]$ - **proportional to $h^3$**

This means thick shells are MUCH stiffer in bending!

### Q5: Why do we need drilling rotation stabilization?
**A**: Without it, the stiffness matrix is **singular** because:
- Shell theory uses 5 DOFs (no rotation about normal)
- Global 3D system has 6 DOFs per node
- The 6th DOF (drilling rotation) has zero stiffness → singular matrix

Solution: Add small artificial stiffness to drilling DOF (~$10^{-3} \times Eh$)

### Q6: What's the difference between this and a beam element?
**A**:
| Feature | Beam | Shell |
|:--------|:-----|:------|
| Nodes per element | 2 | 4 |
| Geometry | 1D (line) | 2D (surface) |
| Components | $K_a + K_b + K_s + K_t$ | $K_m + K_b + K_s$ |
| Local DOFs | 6 | 5 (+1 drilling) |

### Q7: If I double the thickness $h$, what happens to stiffness?
**A**:
- **Membrane**: Doubles ($K_m \propto h$)
- **Bending**: Increases by **8× ($K_b \propto h^3$)**
- **Shear**: Doubles ($K_s \propto h$)

Bending stiffness is VERY sensitive to thickness!

### Q8: What is `R_ele` and why is it 20×24?
**A**: Maps between local and global DOFs:
- **20 rows**: 5 local DOFs × 4 nodes
- **24 columns**: 6 global DOFs × 4 nodes

The drilling DOF (6th global) is excluded from local formulation.

---

## 7. Important Theoretical Concepts

### Shell Theory Components

| Component | Physical Behavior | Thickness Dependence |
|:----------|:------------------|:---------------------|
| **Membrane** | In-plane stretching | $K_m \propto h$ |
| **Bending** | Out-of-plane curvature | $K_b \propto h^3$ |
| **Shear** | Transverse shear | $K_s \propto h$ |

### Integration Strategy

Similar to beams, shells suffer from **locking**:
- **Shear locking**: Use 1-point integration for $K_s$
- **Membrane locking**: Not present in this formulation (bilinear elements OK for moderate distortion)

### Coordinate Transformations

- **Local**: $(x', y', z')$ aligned with element surface
- **Global**: $(x, y, z)$ fixed in space
- **Transformation**: $\mathbf{K}_{global} = \mathbf{R}^T \mathbf{K}_{local} \mathbf{R}$

---

## 8. Summary Flowchart

```
1. Load mesh (xn, Tn, Tm) - 4-node quads
2. Define BCs (clamped root) and loads
3. FOR each element:
   ├─ Compute local axes (i, j, k) from node positions
   ├─ Build rotation matrix R_ele (20×24)
   ├─ SHEAR stiffness Ks:
   │  └─ 1 Gauss point (center) - reduced integration
   ├─ MEMBRANE stiffness Km:
   │  └─ 4 Gauss points (2×2 quadrature)
   ├─ BENDING stiffness Kb:
   │  └─ 4 Gauss points (2×2 quadrature)
   ├─ MASS matrix Me:
   │  └─ 4 Gauss points (2×2 quadrature)
   └─ Assemble: K_el = Ks + Km + Kb
4. Add DRILLING stabilization Kr:
   └─ Penalty stiffness for θ_z to prevent singularity
5. Apply BCs → partition into free/prescribed DOFs
6. STATIC: Solve K*U = F
7. MODAL: Solve K*Φ = ω²*M*Φ
8. Post-process and visualize
```

---

## 9. Comparison: Beam vs Shell

| Aspect | MAIN_beam.m | MAIN_shell.m |
|:-------|:------------|:-------------|
| **Element Type** | 2-node line | 4-node quadrilateral |
| **Theory** | Timoshenko beam | Flat shell (membrane + plate) |
| **DOFs/node** | 6 | 6 (5 physical + 1 drilling) |
| **Stiffness** | $K_a + K_b + K_s + K_t$ | $K_m + K_b + K_s + K_r$ |
| **Shape Functions** | Linear (1D) | Bilinear (2D) |
| **Shear Integration** | 1 point (center) | 1 point (center) |
| **Other Integration** | 1 point (exact) | 4 points (2×2 Gauss) |
| **Locking Issue** | Shear locking | Shear locking |
| **Special Feature** | Torsion stiffness | Drilling stabilization |

Both codes share the **same fundamental strategy**: Use reduced integration on shear terms to prevent locking!
