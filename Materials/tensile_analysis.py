import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import glob
import os

# ----------------------------------------------------------
# CONFIGURACIÓN
# ----------------------------------------------------------
# Patrón de archivos a leer (todos los CSV del directorio actual)
file_pattern = "*.csv"

# Define geometría (A0 [mm²], L0 [mm]) para cada probeta
# Ajusta estos valores según tus muestras reales
geometry = {
    "PP":          {"A0": 40.0, "L0": 80.0},
    "PP_60":       {"A0": 40.0, "L0": 80.0},
    "Twintex_TD":  {"A0": 57.5, "L0": 170.0},
    "Twintex_MD":  {"A0": 57.5, "L0": 170.0},
}

# Parámetros de suavizado y rango elástico
window = 31
poly = 3
eps_min, eps_max = 0.0002, 0.005
offset = 0.002  # 0.2%
# ----------------------------------------------------------

# Lista de archivos
files = glob.glob(file_pattern)
if not files:
    raise FileNotFoundError("No se encontraron archivos .csv en el directorio actual.")

results = []
plt.figure(figsize=(8,6))

for file in files:
    # Nombre base del archivo sin extensión
    name = os.path.splitext(os.path.basename(file))[0]
    print(f"\nProcesando: {name}")

    # Leer datos (usa delimitador ';')
    df = pd.read_csv(file, sep=';', engine='python')

    # Buscar columnas automáticamente
    load_col, elong_col = None, None
    for col in df.columns:
        col_lower = str(col).lower()
        if "load" in col_lower or "force" in col_lower:
            load_col = col
        if "elong" in col_lower or "extens" in col_lower or "displ" in col_lower:
            elong_col = col

    if load_col is None or elong_col is None:
        raise ValueError(f"No se encontraron columnas adecuadas en {file}. Columnas detectadas: {list(df.columns)}")

    # Recuperar geometría
    if name in geometry:
        A0 = geometry[name]["A0"]
        L0 = geometry[name]["L0"]
    else:
        print(f"⚠️ No hay datos de geometría para {name}, se usarán valores por defecto (A0=30, L0=50).")
        A0, L0 = 30.0, 50.0

    # Calcular σ y ε
    df["sigma"] = df[load_col] / A0
    df["eps"] = df[elong_col] / L0

    # Suavizado para reducir ruido
    if len(df) > window:
        df["sigma_s"] = savgol_filter(df["sigma"], window_length=window, polyorder=poly)
    else:
        df["sigma_s"] = df["sigma"]

    # --- Módulo E (ajuste lineal inicial) ---
    mask = (df["eps"] >= eps_min) & (df["eps"] <= eps_max)
    X = df.loc[mask, "eps"].values
    y = df.loc[mask, "sigma_s"].values
    if len(X) > 2:
        E, b = np.polyfit(X, y, 1)
    else:
        E, b = np.nan, np.nan

    # --- Cedencia 0.2% offset ---
    sigma_offset = E * (df["eps"] - offset)
    diff = df["sigma_s"] - sigma_offset
    k = np.where(diff >= 0)[0][0] if np.any(diff >= 0) else None
    if k is not None and k > 0:
        t = diff.iloc[k-1] / (diff.iloc[k-1] - diff.iloc[k])
        eps_y = df["eps"].iloc[k-1] + t * (df["eps"].iloc[k] - df["eps"].iloc[k-1])
        sigma_y = df["sigma_s"].iloc[k-1] + t * (df["sigma_s"].iloc[k] - df["sigma_s"].iloc[k-1])
    else:
        eps_y, sigma_y = np.nan, np.nan

    # --- Resistencia última ---
    idx_u = df["sigma_s"].idxmax()
    sigma_u = df.loc[idx_u, "sigma_s"]
    eps_u = df.loc[idx_u, "eps"]

    # --- Deformación a rotura ---
    eps_f = df["eps"].iloc[-1]

    # Guardar resultados
    results.append({
        "Material": name,
        "E (GPa)": E/1000,
        "σ₀.₂ (MPa)": sigma_y,
        "ε_y (%)": eps_y*100,
        "σ_u (MPa)": sigma_u,
        "ε_u (%)": eps_u*100,
        "ε_f (%)": eps_f*100
    })

    # Graficar curva
    plt.plot(df["eps"], df["sigma_s"], label=name)

# --- Gráfico final ---
plt.xlabel("Strain ε [-]")
plt.ylabel("Stress σ [MPa]")
plt.title("Stress–Strain Curves")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("stress_strain_curves.png", dpi=300)
plt.show()

# --- Tabla resumen ---
table = pd.DataFrame(results)
table.to_csv("tensile_results.csv", index=False)
print("\n=== Resultados ===")
print(table.to_string(index=False))
print("\n📁 Resultados guardados en 'tensile_results.csv' y 'stress_strain_curves.png'")
