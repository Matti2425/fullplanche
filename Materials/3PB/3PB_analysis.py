import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pathlib import Path

# ===============================
# CONFIGURACIÓN INICIAL
# ===============================

# Carpetas de entrada/salida
input_folder = Path("data")
output_folder = Path("results")
output_folder.mkdir(exist_ok=True)

# Dimensiones de las muestras (mm)
sample_dimensions = {
    "PP": {"L": 80, "b": 10, "h": 4},
    "PP_Short_GF": {"L": 80, "b": 10, "h": 4},
    "Twintex": {"L": 80, "b": 20, "h": 2.3},
    "Sandwich": {"L": 80, "b": 40, "h": 19.5},
}

# Rangos de strain para el ajuste lineal (ε_min, ε_max)
fit_ranges = {
    "PP": (0.0015, 0.0090),
    "PP_Short_GF": (0.0, 0.008),
    "Twintex": (0.0005, 0.0025),
    "Sandwich": (0.005, 0.015),
}

# ===============================
# FUNCIÓN DE PROCESADO
# ===============================
def process_three_point_bending(file_path, name, L, b, h, eps_min, eps_max):
    # Leer CSV
    df = pd.read_csv(file_path)
    df.columns = df.columns.str.strip()
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Variables base
    F = df["Load"]        # N
    δ = df["Deflection"]  # mm

    # Cálculo de esfuerzo y deformación
    σ = 3 * F * L / (2 * b * h**2)      # MPa si F en N y dimensiones en mm
    ε = (6 * h * δ) / (L**2)            # adimensional

    df["Stress"] = σ
    df["Strain"] = ε

    # Selección de rango lineal
    mask = (ε >= eps_min) & (ε <= eps_max)
    slope, intercept, r, p, stderr = linregress(ε[mask], σ[mask])
    E = slope  # Módulo de Young (MPa)

    # Energía absorbida (área bajo la curva completa)
    energy_absorbed = np.trapezoid(σ, ε)  # MJ/m³ si σ en MPa (1 MPa = 1 N/mm²)

    # Mostrar resultados
    print(f"\n{'='*60}")
    print(f"Resultados - {name}")
    print(f"{'='*60}")
    print(f"L = {L} mm, b = {b} mm, h = {h} mm")
    print(f"Rango de ajuste: ε ∈ [{eps_min:.4f}, {eps_max:.4f}]")
    print(f"Módulo de Young (E): {E:.2f} MPa")
    print(f"Energía absorbida: {energy_absorbed:.4f} MJ/m³")
    print(f"R² = {r**2:.5f}")
    print()

    # Gráfica
    plt.figure(figsize=(8,6))
    plt.plot(ε, σ, label="Experimental data", color="steelblue")

    # Recta de regresión solo en el rango correspondiente
    ε_fit = np.linspace(eps_min, eps_max, 100)
    σ_fit = slope * ε_fit + intercept
    plt.plot(ε_fit, σ_fit, "r--", label="Linear regression")

    plt.xlabel("Strain (ε)")
    plt.ylabel("Stress (MPa)")
    plt.title(f"3-point bending curve - {name}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(output_folder / f"{name}_curve.png", dpi=300)
    plt.close()

    # Guardar CSV procesado
    df.to_csv(output_folder / f"{name}_processed.csv", index=False)

    return {
        "Sample": name,
        "E (MPa)": E,
        "R²": r**2,
        "Energy (MJ/m³)": energy_absorbed
    }

# ===============================
# PROCESAR TODOS LOS ARCHIVOS
# ===============================
results = []
for name in sample_dimensions:
    dims = sample_dimensions[name]
    eps_min, eps_max = fit_ranges[name]
    file_path = input_folder / f"{name}.csv"

    results.append(
        process_three_point_bending(
            file_path, name,
            dims["L"], dims["b"], dims["h"],
            eps_min, eps_max
        )
    )

# ===============================
# GRÁFICO COMPARATIVO DE TODAS LAS MUESTRAS
# ===============================
plt.figure(figsize=(8,6))

for name in sample_dimensions:
    df = pd.read_csv(output_folder / f"{name}_processed.csv")
    plt.plot(df["Strain"], df["Stress"], label=name)

plt.xlabel("Strain (ε)")
plt.ylabel("Stress (MPa)")
plt.title("Global Stress–Strain Comparison")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(output_folder / "comparison_stress_strain.png", dpi=300)
plt.show()

# ===============================
# EXPORTAR TABLA FINAL
# ===============================
summary = pd.DataFrame(results)
summary.to_csv(output_folder / "summary_results.csv", index=False)

print("Resumen general:")
print(summary.to_string(index=False))
print("\nArchivos exportados en carpeta 'results/'")
