# MAIN_beam.m - Complete Study Guide

**Comprehensive explanation of the Finite Element Analysis implementation for 3D Timoshenko Beam elements.**

---

## 1. Overview: What Does This Code Do?

This script performs **Finite Element Analysis (FEA)** of a 3D beam structure using **Linear Timoshenko Beam Elements** as presented in **CE05 Beams.pdf**. 

### Purpose
1. **Static Analysis**: Calculates deformations under applied loads (bending or torsion)
2. **Modal Analysis**: Determines natural frequencies and vibration mode shapes

### Key Equation (PDF page 27)
$$\mathbf{M}\ddot{\mathbf{u}} + (\mathbf{K}_a + \mathbf{K}_b + \mathbf{K}_s + \mathbf{K}_t)\mathbf{u} = \mathbf{f}$$

Where:
- $\mathbf{K}_a$: Axial stiffness (stretching)
- $\mathbf{K}_b$: Bending stiffness (curvature)
- $\mathbf{K}_s$: Shear stiffness (transverse deformation)
- $\mathbf{K}_t$: Torsion stiffness (twisting)
- $\mathbf{M}$: Mass matrix (inertia)

---

## 2. Code Structure Walkthrough

### Part I: Initialization & Data (Lines 12-35)

**Lines 12-13**: Problem setup
```matlab
problem_type = "torsion"; % bending or torsion
Nm = 6;                   % Number of modes to compute
```

**Lines 16-19**: Material properties
```matlab
E = 75e9;              % Young's Modulus (Pa)
v_poisson = 0.3;       % Poisson's ratio
G = E/(2*(1+v_poisson)); % Shear Modulus (Pa)
rho = 2200;            % Density (kg/m³)
```

**Lines 20-35**: Geometric and cross-sectional properties
- Beam dimensions: chord `c`, span `b`
- Cross-section properties: Area `A`, moments of inertia `Iy`, `Iz`, torsional constant `J`
- **Correction factors**: `ky`, `kz`, `kt` - these are crucial for Timoshenko theory to adjust shear and torsion behavior

### Part II: Preprocessing (Lines 37-73)

**Lines 40-46**: Load mesh data
```matlab
load('beam.mat','xn','Tn','Tm');
% xn: Nodal coordinates [N_nodes × 3]
% Tn: Element connectivity [N_elements × 2]
% Tm: Material IDs [N_elements × 1]
```

**Lines 51-56**: Boundary conditions (Cantilever beam)
```matlab
Up = [0 1 1;  % Node 1 is CLAMPED (all 6 DOFs = 0)
      0 1 2;
      0 1 3;
      0 1 4;
      0 1 5;
      0 1 6];
```

**Lines 59-64**: External loads
```matlab
switch(problem_type)
    case "bending"
        Fe = [1 N_nodes 3];  % Force in Z-direction at tip
    case "torsion"
        Fe = [1 N_nodes 4];  % Moment about X-axis at tip
end
```

### Part III: Solver - Matrix Assembly (Lines 90-196)

This is the core loop that builds the global stiffness `K` and mass `M` matrices.

---

## 3. Theoretical Foundation → Code Implementation

### General FEM Formula (PDF page 27)
$$\mathbf{K}_{\text{component}}^e = \mathbf{R}_e^T \int_{\ell_e} \mathbf{B}_{\text{component}}^{\prime e T} \bar{\mathbf{C}}_{\text{component}}^{\prime e} \mathbf{B}_{\text{component}}^{\prime e} \, d\ell \, \mathbf{R}_e$$

**Components**:
- $\mathbf{R}_e$: Rotation matrix (local → global coordinate transformation)
- $\mathbf{B}$: Strain-displacement matrix (relates node displacements to strains)
- $\bar{\mathbf{C}}$: Constitutive matrix (material properties: $E$, $G$, etc.)

---

### A. Local Coordinate System (Lines 103-116)

**Code**:
```matlab
% Compute local axes
x1 = xn(Tn(e,1),:)';
x2 = xn(Tn(e,2),:)';
Le = norm(x2 - x1);

i = (x2 - x1)/Le;              % x' axis (along beam)
j = j_orientations(Tm(e,1),:)'; % y' axis (user-defined)
k = cross(i,j);                % z' axis (right-hand rule)

% Build rotation matrices
R3 = [i j k];           % 3×3 rotation matrix
Rp = blkdiag(R3, R3);   % 6×6 (for one node: 3 translations + 3 rotations)
Re = blkdiag(Rp, Rp);   % 12×12 (for two nodes)
```

**Explanation**: 
- Each beam element has its own local coordinate system aligned with the beam axis
- `blkdiag` creates block diagonal matrices to rotate both translations and rotations
- The $12 \times 12$ matrix `Re` handles both nodes of the element

---

### B. Shape Functions & Derivatives (Line 119)

**PDF Theory (page 28)**:

Natural coordinates $\xi \in [-1, 1]$:
$$N_1(\xi) = \frac{1}{2}(1 - \xi), \quad N_2(\xi) = \frac{1}{2}(1 + \xi)$$

Derivatives in local coordinates:
$$N_{,x'}^1 = \frac{-1}{\ell_e}, \quad N_{,x'}^2 = \frac{+1}{\ell_e}$$

**Code**:
```matlab
Nx = [-1/Le; 1/Le];  % Shape function derivatives
```

---

### C. Axial Stiffness $\mathbf{K}_a$ (Lines 118-125)

**PDF Theory**: 
$$\mathbf{K}_a^e = \mathbf{R}_e^T \int_{\ell_e} \mathbf{B}_a^{\prime e T} (EA) \mathbf{B}_a^{\prime e} \, d\ell \, \mathbf{R}_e$$

**Code**:
```matlab
% Strain-displacement matrix (1×12)
Ba_row = [Nx(1) 0 0 0 0 0 Nx(2) 0 0 0 0 0];
Ba = sparse(Ba_row);

% Constitutive matrix (scalar)
Ca = E_vector(Tm(e,1))*A_vector(Tm(e,1));  % EA

% Element stiffness
Ka = Le * (Re' * (Ba' * Ca * Ba) * Re);
```

**Physical Meaning**: 
- $\mathbf{B}_a$ relates nodal displacements to axial strain $\varepsilon = \frac{du}{dx}$
- Operates only on DOFs 1 and 7 (axial displacements)
- Integration is exact because $\mathbf{B}_a$ is constant → multiply by `Le`

---

### D. Bending Stiffness $\mathbf{K}_b$ (Lines 127-134)

**PDF Theory**: 
$$\mathbf{K}_b^e = \mathbf{R}_e^T \int_{\ell_e} \mathbf{B}_b^{\prime e T} \begin{bmatrix} EI_y & 0 \\ 0 & EI_z \end{bmatrix} \mathbf{B}_b^{\prime e} \, d\ell \, \mathbf{R}_e$$

**Code**:
```matlab
% Strain-displacement matrix (2×12)
Bb = sparse([...
    0 0 0 0 Nx(1) 0   0 0 0 0 Nx(2) 0;     % κ_z = dθ_y/dx
    0 0 0 0 0 Nx(1)  0 0 0 0 0 Nx(2)]);   % κ_y = dθ_z/dx

% Constitutive matrix (2×2 diagonal)
Cb = E_vector(Tm(e,1)) * diag([Iy_vector(Tm(e,1)) Iz_vector(Tm(e,1))]);

% Element stiffness
Kb = Le * (Re' * (Bb' * Cb * Bb) * Re);
```

**Physical Meaning**:
- $\mathbf{B}_b$ relates nodal rotations to curvature $\kappa = \frac{d\theta}{dx}$
- Operates on DOFs 5,6,11,12 (rotations about $y'$ and $z'$ axes)
- Note: In Timoshenko theory, curvature depends on rotation derivatives (not displacement second derivatives)

---

### E. Shear Stiffness $\mathbf{K}_s$ (Lines 136-146) ⚠️ **CRITICAL FOR EXAM!**

**PDF Theory**: 
$$\mathbf{K}_s^e = \mathbf{R}_e^T \int_{\ell_e} \mathbf{B}_s^{\prime e T} \begin{bmatrix} k_y GA & 0 \\ 0 & k_z GA \end{bmatrix} \mathbf{B}_s^{\prime e} \, d\ell \, \mathbf{R}_e$$

**PDF Warning (page 28)**:
> "For the shear component, even though 2 Gauss points provide exact integration, **the element suffers from shear-locking for slender beams**. A **sub-integration with only 1 Gauss point avoids this problem**."

**Code**:
```matlab
N_G = 1/2;  % Shape function value at midpoint (ξ=0)

% Strain-displacement matrix (2×12)
Bs = sparse([...
    0 Nx(1) 0 0 0 -N_G  0 Nx(2) 0 0 0 -N_G;  % γ_y = ∂u_y/∂x - θ_z
    0 0 Nx(1) 0 N_G 0   0 0 Nx(2) 0 N_G 0]); % γ_z = ∂u_z/∂x - θ_y

% Constitutive matrix (2×2 diagonal)
Cs = G_vector(Tm(e,1))*A_vector(Tm(e,1))* ...
     diag([ky_vector(Tm(e,1)) kz_vector(Tm(e,1))]);

% Element stiffness
Ks = Le * (Re' * (Bs' * Cs * Bs) * Re);
```

**Physical Meaning**:
- Shear strain in Timoshenko: $\gamma = \frac{\partial w}{\partial x} - \theta$ (difference between slope of deflection and rotation)
- `N_G = 1/2` implements **reduced integration** (1-point Gauss at element center)
- This prevents **shear locking** - a numerical artifact that makes slender beams artificially stiff

**Exam Question**: Why not use 2 Gauss points for exact integration?
**Answer**: Full integration causes shear locking. The reduced integration (1-point) gives much more accurate results for thin beams, as explicitly recommended in the PDF.

---

### F. Torsion Stiffness $\mathbf{K}_t$ (Lines 148-154)

**PDF Theory**: 
$$\mathbf{K}_t^e = \mathbf{R}_e^T \int_{\ell_e} \mathbf{B}_t^{\prime e T} (k_t GJ) \mathbf{B}_t^{\prime e} \, d\ell \, \mathbf{R}_e$$

**Code**:
```matlab
% Strain-displacement matrix (1×12)
Bt_row = [0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0 0];
Bt = sparse(Bt_row);

% Constitutive matrix (scalar)
Ct = G_vector(Tm(e,1))*J_vector(Tm(e,1))*kt_vector(Tm(e,1));

% Element stiffness
Kt = Le * (Re' * (Bt' * Ct * Bt) * Re);
```

**Physical Meaning**:
- Relates nodal rotations to twist rate $\beta = \frac{d\theta_x}{dx}$
- Operates on DOFs 4 and 10 (rotations about beam axis)
- $J$ is the torsional constant (polar moment of inertia for circular sections)

---

### G. Mass Matrix $\mathbf{M}$ (Lines 159-183)

**PDF Theory**: 
$$\mathbf{M}^e = \mathbf{R}_e^T \int_{\ell_e} \mathbf{N}_e^T \boldsymbol{\rho}'^e \mathbf{N}_e \, d\ell \, \mathbf{R}_e$$

**PDF Recommendation (page 28)**: "For the mass matrix, **2 Gauss points provide exact integration**."

**Code**:
```matlab
% Gauss quadrature setup
Xi = [-1/sqrt(3); 1/sqrt(3)];  % Gauss points: ±0.577
w  = [1; 1];                   % Weights

% Inertia matrix (6×6)
rho_mat = rho * blkdiag(A*speye(3), J, Iy, Iz);

Me = sparse(12,12);

% Integrate using 2 Gauss points
for kk = 1:2
    N = [(1 - Xi(kk))/2; (1 + Xi(kk))/2];  % Shape functions
    
    % Build N matrix (6×12)
    Nmat = sparse(6,12);
    for k = 1:6
        Nmat(k, k)     = N(1);
        Nmat(k, 6 + k) = N(2);
    end

    Me = Me + (w(kk)*Le/2) * (Re' * (Nmat' * rho_mat * Nmat) * Re);
end
```

**Physical Meaning**:
- Creates a **consistent mass matrix** (not lumped)
- Uses 2-point Gauss quadrature for exact integration
- Captures full inertial distribution of the beam

---

### H. Global Assembly (Lines 186-196)

**Code**:
```matlab
% Compute DOF indices for this element
Idof = zeros(12,1);
for j = 1:6
    Idof(j)     = 6*(Tn(e,1)-1) + j;  % Node 1 DOFs
    Idof(6 + j) = 6*(Tn(e,2)-1) + j;  % Node 2 DOFs
end

% Add element matrices to global system
K(Idof,Idof) = K(Idof,Idof) + Ka + Kb + Ks + Kt;
M(Idof,Idof) = M(Idof,Idof) + Me;
```

**Explanation**: Element matrices are placed in the global system at positions corresponding to the element's DOFs.

---

## 4. Solution Methods

### Part IV: Static Analysis (Lines 200-216)

**Code**:
```matlab
% Build force vector
F = zeros(Ndof, 1);
F_global = (F_node - 1) * 6 + F_dof;
F(F_global) = Fe(1,1);

% Apply boundary conditions
prescribed_dofs = (Up(:,2) - 1) * 6 + Up(:,3);
free_dofs = setdiff(1:Ndof, prescribed_dofs);

% Solve: K*U = F
U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);
```

**Explanation**: 
- Only solve for "free" DOFs (prescribing $U = 0$ where clamped)
- The backslash `\` operator solves the linear system $\mathbf{K}_{ff}\mathbf{U}_f = \mathbf{F}_f$

### Part V: Modal Analysis (Lines 218-226)

**Code**:
```matlab
n_modes = min(Nm, length(free_dofs));
[Phi, Omega2] = eigs(K(free_dofs, free_dofs), M(free_dofs, free_dofs), n_modes, 'sm');
frequencies = sqrt(diag(Omega2)) / (2*pi);
```

**Theory**: Solves the generalized eigenvalue problem:
$$\mathbf{K}\boldsymbol{\phi} = \omega^2 \mathbf{M}\boldsymbol{\phi}$$

Where:
- $\boldsymbol{\phi}$: Mode shape (eigenvector)
- $\omega^2$: Squared natural frequency (eigenvalue)
- `'sm'`: Finds smallest magnitude eigenvalues (fundamental modes)

---

## 5. Post-Processing (Lines 234-333)

### Extract Results
```matlab
u_x = U(1:6:end);       % Axial displacement
u_y = U(2:6:end);       % Lateral displacement
u_z = U(3:6:end);       % Vertical displacement
theta_x = U(4:6:end);   % Twist
theta_y = U(5:6:end);   % Bending rotation
theta_z = U(6:6:end);   % Bending rotation
```

### Plotting
- **Lines 246-295**: Static deformation plots
- **Lines 306-333**: Mode shape visualization

---

## 6. Key Variable Reference

| Variable | Dimensions | Description |
|:---------|:-----------|:------------|
| `xn` | $N_{nodes} \times 3$ | Nodal coordinates $(x, y, z)$ |
| `Tn` | $N_{elem} \times 2$ | Element connectivity [node1, node2] |
| `K` | $6N \times 6N$ | Global stiffness matrix (sparse) |
| `M` | $6N \times 6N$ | Global mass matrix (sparse) |
| `U` | $6N \times 1$ | Displacement vector $[u_x, u_y, u_z, \theta_x, \theta_y, \theta_z]$ |
| `Nx` | $2 \times 1$ | Shape function derivatives $[-1/L, 1/L]$ |
| `Re` | $12 \times 12$ | Rotation matrix (local → global) |

---

## 7. Exam Q&A

### Q1: What happens if I change `problem_type` from "torsion" to "bending"?
**A**: The force application changes:
- **Bending**: Force on DOF 3 ($u_z$) → vertical deflection
- **Torsion**: Moment on DOF 4 ($\theta_x$) → twisting

### Q2: How is this beam constrained?
**A**: **Cantilever** (clamped-free):
- Node 1: All 6 DOFs fixed ($u_x = u_y = u_z = \theta_x = \theta_y = \theta_z = 0$)
- Last node: Free to deform

To make it **simply supported** (pinned-pinned): Fix translations at both ends, leave rotations free.

### Q3: Why use `N_G = 1/2` for shear stiffness?
**A**: This is the **reduced integration fix for shear locking**:
- Full integration (2 Gauss points) makes slender beams too stiff
- 1-point integration (at element center where $\xi=0$, $N=1/2$) prevents this
- Explicitly recommended in PDF page 28

### Q4: What confirms this is Timoshenko (not Euler-Bernoulli) theory?
**A**: Three indicators:
1. Presence of shear stiffness matrix $\mathbf{K}_s$
2. Shear correction factors $k_y, k_z$
3. Independent rotation variables (not kinematically constrained to displacement derivatives)

### Q5: What does `blkdiag(R3, R3)` do?
**A**: Creates block diagonal matrix:
$$\begin{bmatrix} R_3 & 0 \\ 0 & R_3 \end{bmatrix}$$
Applies the same rotation to both translations and rotations (each node has 6 DOFs).

### Q6: How many Gauss points are used?
**A**:
- **Stiffness** ($K_a, K_b, K_s, K_t$): 1 point (constant strain → integrated as `Le * ...`)
- **Mass** ($M$): 2 points at $\xi = \pm 1/\sqrt{3}$

### Q7: If I double Young's Modulus ($E$), what happens?
**A**:
- **Static deflections**: Halved (stiffness doubles)
- **Natural frequencies**: Increase by $\sqrt{2}$ (since $\omega \propto \sqrt{K/M}$)

### Q8: What if I refine the mesh (more elements)?
**A**: The code adapts automatically:
- `N_elements` and `N_nodes` are computed from loaded data
- Convergence: Results become more accurate (FEM converges from below for displacements)

---

## 8. Important Theoretical Concepts

### Timoshenko vs. Euler-Bernoulli
| Feature | Timoshenko | Euler-Bernoulli |
|:--------|:-----------|:----------------|
| Shear deformation | ✓ Included | ✗ Neglected |
| Best for | Thick/short beams | Slender beams |
| DOFs per node | 6 (3D) | 4 (2D) |
| Rotation-displacement relation | Independent | $\theta = -w'$ |

### Why Reduced Integration?
- **Shear locking**: Over-stiff numerical behavior in thin beams when using full integration
- **Solution**: Underintegrate the shear term (1 Gauss point instead of 2)
- **Trade-off**: Slight accuracy loss, but massive improvement for slender elements

### Coordinate Systems
- **Local** ($x', y', z'$): Aligned with beam axis, simplifies formulation
- **Global** ($x, y, z$): Fixed reference frame
- **Transformation**: $\mathbf{K}_{global} = \mathbf{R}^T \mathbf{K}_{local} \mathbf{R}$

---

## 9. Summary Flowchart

```
1. Load mesh (xn, Tn, Tm)
2. Define BCs (Up) and loads (Fe)
3. FOR each element:
   ├─ Compute local axes (i, j, k)
   ├─ Build rotation matrix Re
   ├─ Compute element matrices:
   │  ├─ Ka (axial)
   │  ├─ Kb (bending)
   │  ├─ Ks (shear - reduced integration!)
   │  ├─ Kt (torsion)
   │  └─ Me (mass - 2-point Gauss)
   └─ Assemble into global K and M
4. Apply BCs → partition into free/prescribed DOFs
5. STATIC: Solve K*U = F
6. MODAL: Solve K*Φ = ω²*M*Φ
7. Post-process and plot results
```
