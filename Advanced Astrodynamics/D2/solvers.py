# =========================================================
#                Kepler Equation Solvers
# =========================================================

import numpy as np

def newton_solver(M, e, E0, tol=1e-12, max_iter=30):
    """
    Newton–Raphson method for Kepler's equation with numerical safeguards.
    """
    E = E0
    history = [E0]
    for i in range(max_iter):
        f = E - e * np.sin(E) - M
        f_prime = 1 - e * np.cos(E)

        # Safety check: avoid division by near-zero
        if abs(f_prime) < 1e-14:
            # small random perturbation or fallback
            f_prime = np.sign(f_prime) * 1e-14 if f_prime != 0 else 1e-14

        dE = -f / f_prime
        E += dE
        history.append(E)

        if abs(dE) < tol:
            return E, i + 1, np.array(history)
    return E, max_iter, np.array(history)


def fixed_point_solver(M, e, E0, tol=1e-12, max_iter=200):
    """Fixed-point iteration solver: E = M + e*sin(E)."""
    E = E0
    history = [E0]
    for i in range(max_iter):
        E_next = M + e * np.sin(E)
        history.append(E_next)
        if abs(E_next - E) < tol:
            return E_next, i + 1, np.array(history)
        E = E_next
    return E, max_iter, np.array(history)

def halley_solver(M, e, E0, tol=1e-12, max_iter=30):
    """
    Halley's method for solving Kepler's equation with stability control.
    """
    E = E0
    history = [E0]
    for i in range(max_iter):
        f = E - e * np.sin(E) - M
        f1 = 1 - e * np.cos(E)
        f2 = e * np.sin(E)

        denom = f1 - 0.5 * f * f2 / f1 if abs(f1) > 1e-14 else 1e-14

        # Avoid dividing by extremely small values
        if abs(denom) < 1e-14:
            denom = np.sign(denom) * 1e-14 if denom != 0 else 1e-14

        dE = -f / denom
        E += dE
        history.append(E)

        if abs(dE) < tol:
            return E, i + 1, np.array(history)
    return E, max_iter, np.array(history)

