import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import solvers
import re

# =========================================================
#                Convergence Analysis
# =========================================================

def convergence_map(solver, E0_func, e_vals, M_vals, tol=1e-10):
    """Compute number of iterations for convergence across a grid."""
    counts = np.zeros((len(M_vals), len(e_vals)))
    for i, M in enumerate(M_vals):
        for j, e in enumerate(e_vals):
            E0 = E0_func(M, e)
            _, n_iter, _ = solver(M, e, E0, tol=tol)
            counts[i, j] = n_iter
    return counts


def plot_convergence_map(e_vals, M_vals_deg, iterations, title, filename):
    """
    Plot a convergence map and save it safely (works on Windows too).
    """
    # Sanitize filename (replace forbidden characters)
    safe_filename = re.sub(r'[<>:"/\\|?*+–]', '_', str(filename))

    plt.figure(figsize=(8, 5))
    plt.imshow(
        iterations,
        extent=[M_vals_deg.min(), M_vals_deg.max(), e_vals.min(), e_vals.max()],
        origin='lower',
        aspect='auto',
        cmap='jet'
    )
    plt.colorbar(label='Iterations to converge')
    plt.xlabel('Mean anomaly M [deg]')
    plt.ylabel('Eccentricity e')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(safe_filename, dpi=300)
    plt.close()

# =========================================================
#                Comparative Iteration Plot
# =========================================================

def plot_iterations(M, e, E0, solvers, names, filename):
    plt.figure(figsize=(8, 6))

    for solver, name in zip(solvers, names):
        _, _, hist = solver(M, e, E0)
        plt.semilogy(range(len(hist)), np.abs(hist - hist[-1]), label=name)

    plt.xlabel("Iteration")
    plt.ylabel("Absolute Error |Eₙ - E*|")
    plt.title(f"Convergence comparison (e={e:.3f}, M={np.degrees(M):.1f}°)")
    plt.legend()
    plt.grid(True, which='both', ls='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


# =========================================================
#                Main Execution
# =========================================================

if __name__ == "__main__":
    results_dir = Path("./results")
    results_dir.mkdir(exist_ok=True)

    e_vals = np.linspace(0, 1, 500)
    M_vals = np.linspace(0, 2*np.pi, 360)
    M_vals_deg = np.degrees(M_vals)

    solvers = [solvers.newton_solver, solvers.fixed_point_solver, solvers.halley_solver]
    solver_names = ["Newton–Raphson", "Fixed-Point", "Halley"]

    E0_funcs = {
        "E0 = M": lambda M, e: M,
        "E0 = π": lambda M, e: np.pi * np.ones_like(M),
        "E0 = M + e*cos(M)": lambda M, e: M + e * np.cos(M),
    }

    # --- Generate convergence maps ---
    for guess_name, E0_func in E0_funcs.items():
        for solver, name in zip(solvers, solver_names):
            counts = convergence_map(solver, E0_func, e_vals, M_vals)
            plot_convergence_map(
                e_vals, M_vals_deg, counts,
                f"{name} convergence ({guess_name})",
                results_dir / f"map_{name.replace(' ', '_')}_{guess_name.replace('=','_')}.png"
            )

    # --- Generate comparative convergence plot ---
    e_test = 0.95
    M_test = np.radians(350)
    E0_test = M_test  # Example starting point
    plot_iterations(
        M_test, e_test, E0_test,
        solvers, solver_names,
        results_dir / "comparison_convergence.png"
    )

    print(f"All results saved in: {results_dir.resolve()}")
