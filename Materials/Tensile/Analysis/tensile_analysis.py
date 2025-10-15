# tensile_analysis_energy_v2.py
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# ---------------- CONFIG ----------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "../data")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "../results")
os.makedirs(RESULTS_DIR, exist_ok=True)  # crea la carpeta si no existe
PATTERN = os.path.join(DATA_DIR, "*.csv")

# Geometría por nombre de fichero (ajusta A0 [mm^2], L0 [mm])
geometry = {
    "PP":          {"A0": 40.0, "L0": 80.0},
    "PP_60":       {"A0": 40.0, "L0": 80.0},
    "Twintex_TD":  {"A0": 57.5, "L0": 170.0},
    "Twintex_MD":  {"A0": 57.5, "L0": 170.0},
}

# Rango para calcular E (por muestra). Si una muestra no está aquí, usamos default_fit_range.
fit_ranges = {
    "PP": (0.0, 0.020),
    "PP_60": (0.0, 0.008),
    "Twintex_MD": (0.0, 0.035),
    "Twintex_TD": (0.0, 0.035),
}
default_fit_range = (0.001, 0.01)

# Suavizado Savitzky-Golay
WINDOW = 31
POLY = 3

# Offset para sigma_0.2
OFFSET = 0.002  # 0.2%

# Salidas
OUT_CSV = os.path.join(SCRIPT_DIR, "tensile_results_energy.csv")
OUT_PLOT = os.path.join(SCRIPT_DIR, "stress_strain_energy.png")
# ----------------------------------------

def detect_columns(df):
    load_col = None
    elong_col = None
    for c in df.columns:
        cl = str(c).lower()
        if ("load" in cl) or ("force" in cl):
            load_col = c
        if ("elong" in cl) or ("displ" in cl) or ("extens" in cl):
            elong_col = c
    return load_col, elong_col

def compute_sigma_eps(df, load_col, elong_col, A0, L0):
    df2 = df[[load_col, elong_col]].dropna().reset_index(drop=True)
    sigma = df2[load_col].values / A0        # N/mm² → MPa
    eps = df2[elong_col].values / L0         # adimensional
    return df2, sigma, eps

def smooth_sigma(sigma, window=WINDOW, poly=POLY):
    if len(sigma) > window:
        w = window if window % 2 == 1 else window - 1
        w = min(w, len(sigma) - 1 if (len(sigma) - 1) % 2 == 1 else len(sigma) - 2)
        if w < 3:
            return sigma
        return savgol_filter(sigma, window_length=w, polyorder=poly)
    else:
        return sigma

def compute_E(eps, sigma_s, fit_range):
    eps_min, eps_max = fit_range
    mask = (eps >= eps_min) & (eps <= eps_max)
    X = eps[mask]
    y = sigma_s[mask]
    if len(X) >= 3:
        coef = np.polyfit(X, y, 1)
        E = float(coef[0])
        intercept = float(coef[1])
        return E, intercept, (X.min(), X.max())
    else:
        return np.nan, np.nan, (np.nan, np.nan)

def compute_offset_yield(eps, sigma_s, E, intercept, offset=OFFSET):
    if np.isnan(E):
        return (np.nan, np.nan)
    sigma_offset = E * (eps - offset) + intercept
    diff = sigma_s - sigma_offset
    idx = np.where(diff >= 0)[0]
    if len(idx) == 0:
        return (np.nan, np.nan)
    k = idx[0]
    if k == 0:
        return (np.nan, np.nan)
    d1 = diff[k-1]; d2 = diff[k]
    if (d1 - d2) == 0:
        t = 0.0
    else:
        t = d1 / (d1 - d2)
    eps_y = eps[k-1] + t*(eps[k] - eps[k-1])
    sigma_y = sigma_s[k-1] + t*(sigma_s[k] - sigma_s[k-1])
    return (float(eps_y), float(sigma_y))

def integrate_energy(eps, sigma_s):
    if len(eps) < 2:
        return np.nan
    U = np.trapezoid(sigma_s, eps)   # MPa ~ MJ/m³
    return float(U)

# --- main ---
files = sorted(glob.glob(PATTERN))
if not files:
    raise FileNotFoundError(f"No .csv files found in {DATA_DIR}")

rows = []
plt.figure(figsize=(9,6))
colors = plt.cm.tab10.colors

for i, path in enumerate(files):
    name = os.path.splitext(os.path.basename(path))[0]
    print(f"\nProcesando: {name}")
    df = pd.read_csv(path, sep=';', engine='python', comment='#')
    df = df.loc[:, df.notna().any()]

    load_col, elong_col = detect_columns(df)
    if load_col is None or elong_col is None:
        print(f"  -> columnas 'Load'/'Elongation' no detectadas en {path}. Columnas: {df.columns.tolist()}")
        continue

    geom = geometry.get(name, {"A0": 30.0, "L0": 50.0})
    A0, L0 = geom["A0"], geom["L0"]

    df2, sigma, eps = compute_sigma_eps(df, load_col, elong_col, A0, L0)
    sigma_s = smooth_sigma(sigma, window=WINDOW, poly=POLY)

    fit_range = fit_ranges.get(name, default_fit_range)
    E, intercept, used_range = compute_E(eps, sigma_s, fit_range)

    eps_y, sigma_y = compute_offset_yield(eps, sigma_s, E, intercept, offset=OFFSET)
    U = integrate_energy(eps, sigma_s)
    U_J = U * A0 * L0 * 1e-3 if not np.isnan(U) else np.nan

    idx_u = int(np.nanargmax(sigma_s))
    sigma_u = float(sigma_s[idx_u])
    eps_u = float(eps[idx_u])
    eps_f = float(eps[-1])

    rows.append({
        "Material": name,
        "E (GPa)": E/1000 if not np.isnan(E) else np.nan,
        "σ₀.₂ (MPa)": sigma_y,
        "ε_y (%)": eps_y*100 if not np.isnan(eps_y) else np.nan,
        "σ_u (MPa)": sigma_u,
        "ε_u (%)": eps_u*100,
        "ε_f (%)": eps_f*100,
        "U (MJ/m³)": U,
        "U_total (J)": U_J
    })

    # Graficar curva principal
    c = colors[i % len(colors)]
    plt.plot(eps, sigma_s, label=name, color=c, linewidth=1.2)
    plt.fill_between(eps, sigma_s, 0.0, alpha=0.12, color=c)

    # Tramo de regresión solo dentro del rango usado
    if not np.isnan(E) and not np.isnan(intercept):
        mask_plot = (eps >= used_range[0]) & (eps <= used_range[1])
        Xplot = eps[mask_plot]
        if len(Xplot) > 1:
            plt.plot(Xplot, E*Xplot + intercept, '--', color=c, linewidth=1.3)

    # Línea offset 0.2% y punto yield
    if not np.isnan(E) and not np.isnan(intercept):
        Xlocal = np.linspace(used_range[0], used_range[1], 100)
        plt.plot(Xlocal, E*(Xlocal - OFFSET) + intercept, ':', color=c, linewidth=0.9)
    if not np.isnan(eps_y) and not np.isnan(sigma_y):
        plt.plot(eps_y, sigma_y, 'x', color=c, markersize=6)

# --- Gráfico final ---
plt.xlabel("Strain ε [-]")
plt.ylabel("Stress σ [MPa]")
plt.title("Stress–Strain curves (smoothed) with elastic fit and energy area")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(OUT_PLOT, dpi=300)
plt.show()

# --- Tabla resumen ---
table = pd.DataFrame(rows)
table.to_csv(OUT_CSV, index=False)

# Mostrar resultados bonitos en terminal
pd.set_option("display.float_format", "{:,.3f}".format)
print("\n=== RESULTADOS ===")
print(table[["Material","E (GPa)","σ₀.₂ (MPa)","ε_y (%)","σ_u (MPa)","ε_u (%)","ε_f (%)","U (MJ/m³)","U_total (J)"]].to_string(index=False))
print("\nResultados guardados en:")
print(f"  → {OUT_CSV}")
print(f"  → {OUT_PLOT}")
