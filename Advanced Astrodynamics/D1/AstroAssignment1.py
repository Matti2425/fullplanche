import numpy as np

# =======================
# Step 1 and 2: Julian Date (JD) and Julian Century (JC)
# =======================
def JD_from_UTC(year, month, day):
    if month <=2:
        year = year - 1
        month = month + 12
    if month > 2:
        year = year
        month = month
    JD = int(365.25*year) + int(30.6001*(month+1)) + day + 1720981.5
    return JD

def JC(JD):
    JC = (JD - 2451545)/36525
    return JC

print("JD (15.09.2017):", JD_from_UTC(2017, 9, 15))
print("JC (15.09.2017):",JC(JD_from_UTC(2017, 9, 15)))
print("----------------------------------------------")

# =======================
# Step 3: Propagate six orbital elements
# =======================
#Aus bekannten Bahnelementen zu einem Ausgangsepoche (z. B. J2000) die Bahnelemente zu einem neuen Zeitpunkt berechnen.

def deg_correction(x): #every angle in 360 degrees system - Winkel in [0,360)
    return x % 360

def propagate_elements(X_0 , X_punkt, JC):
    a = X_0["a"] + X_punkt["a"] * JC  #AU
    e = X_0["e"] + X_punkt["e"] * JC  # -
    i = deg_correction(X_0["i"] + X_punkt["i"] * JC)  # rest is in deg
    omega = deg_correction(X_0["omega"] + X_punkt["omega"] * JC)
    w_bar = deg_correction(X_0["w_bar"] + X_punkt["w_bar"] * JC)
    L = deg_correction(X_0["L"] + X_punkt["L"] * JC)
    return dict(a = a, e = e, i = i, omega = omega, w_bar = w_bar, L = L)

# ---- Tabellenwerte --- (slide: Planetary orbital elements)
X_0_earth = dict(a=1.00000261, e=0.01671123, i=-0.00001531, omega=0.0, w_bar=102.93768193, L=100.46457166)
X_punkt_earth = dict(a=0.00000562, e=-0.00004392, i=-0.01294668, omega=0.0, w_bar=0.32327364, L=35999.37244981)

X_0_saturn = dict(a=9.53667594, e=0.05386179, i=2.48599187, omega=113.66242448, w_bar=92.59887831, L=49.95424423)
X_punkt_saturn = dict(a=-0.00125060, e=-0.00050991, i=0.00193609, omega=-0.28867794, w_bar=-0.41897216, L=1222.49362201)


pe_earth = propagate_elements(X_0_earth, X_punkt_earth, JC(JD_from_UTC(2017, 9, 15)))
pe_saturn = propagate_elements(X_0_saturn, X_punkt_saturn, JC(JD_from_UTC(2017, 9, 15)))
print( "Propagated orbital elements:")
print("Earth: ", pe_earth)
print("Saturn: ", pe_saturn)
print("----------------------------------------------")

# =======================
# Step 4: Compute angular momentum h using semi-major axis a [AU] and eccentricity e [-]
# =======================
#Gravitationsparameter der Sonne mu=1.32712440018*10^20 m^3/s^2  - ask the Prof
k_gauss= 0.01720209895
mu_sonne_AU = k_gauss**2

def angular_momentum (a,e):
    h = np.sqrt(mu_sonne_AU*a*(1-e**2))
    return h

h_earth = angular_momentum(pe_earth["a"], pe_earth["e"])
h_saturn = angular_momentum(pe_saturn["a"], pe_saturn["e"])
print("Angular momentum h in [AU^2/day]:")
print("Earth: ", h_earth, "AU^2/day")
print("Saturn: ", h_saturn, "AU^2/day")
print("----------------------------------------------")

# =======================
# Step 5: Determine the argument of the perihelion w and mean anomaly M at JC
# =======================
# w = w_bar -omega and M=L-w
def w (w_bar, omega):
    return deg_correction(w_bar-omega)
def M (L,w_bar):
    return deg_correction(L-w_bar)

w_earth = w(pe_earth["w_bar"], pe_earth["omega"])
M_earth =M(pe_earth["L"], pe_earth["w_bar"])
w_saturn = w(pe_saturn["w_bar"], pe_saturn["omega"])
M_saturn = M(pe_saturn["L"], pe_saturn["w_bar"])

print("Argument of the perihelion w and mean anomaly M at JC:")
print("Earth: w =", w_earth, "degrees", " AND M =", M_earth, "degrees" )
print("Saturn: w =", w_saturn, "degrees", " AND M =", M_saturn, "degrees" )
print("----------------------------------------------")

# =======================
# Step 6: Calculate eccentric anomaly E solving the Kepler equation M = E - e*sin(E)
# STANDARD SERIES
# =======================
# solving Kepler equation with Standard series - error = O(e^4)

def rad_correction(x):
    return x % (2*np.pi)

def kepler_E_series(M_deg, e):
    M= rad_correction(np.deg2rad(M_deg)) #Es wickelt den Winkel auf den Hauptbereich [0,2π) zurück
    E = M + e*np.sin(M) + 0.5*e**2*np.sin(2*M)+ (e**3/8.0)*(3*np.sin(3*M) - np.sin(M)) + (1/3)*e**4*np.sin(4*M)
    return rad_correction(E)

E_earth = kepler_E_series(M_earth, pe_earth["e"])
E_saturn = kepler_E_series(M_saturn, pe_saturn["e"])

print("Eccentric anomaly: ")
print("E_earth: ", E_earth, "Rad")
print("E_saturn: ", E_saturn, "Rad")
print("----------------------------------------------")

# =======================
# Step 7: Compute true anomaly
# =======================

def true_anomaly(E,e):
    sin_t = np.sqrt(1 - e**2) * np.sin(E) / (1 - e*np.cos(E))
    cos_t = (np.cos(E) - e) / (1 - e*np.cos(E))
    true_anomaly = rad_correction(np.arctan2(sin_t, cos_t))
    return true_anomaly

theta_earth = true_anomaly(E_earth, pe_earth["e"])
theta_saturn = true_anomaly(E_saturn, pe_saturn["e"])
print("True anomaly: ")
print("true_anomaly_earth (theta): ", theta_earth, "Rad")
print("true_anomaly_saturn (theta): ", theta_saturn, "Rad")
print("----------------------------------------------")

# =======================
# Step 8: Compute vector (r,v)
# =======================
# 1) in the perifocal system (Bahnebene, Origin= Sun) using theta

def r_v_pf(h, e, theta):  # h [AU^2/day], e [-], theta [rad]
    # Radius r = (h^2/μ) / (1 + e cosθ)
    r = (h**2 / mu_sonne_AU) / (1.0 + e*np.cos(theta))
    # Positionsvektor (perifokal, PQW)
    r_pf = np.array([r*np.cos(theta), r*np.sin(theta), 0.0])
    # Geschwindigkeitsvektor (perifokal, PQW) – deine kompakte Formel
    v_pf = (mu_sonne_AU / h) * np.array([-np.sin(theta), e + np.cos(theta), 0.0])
    return r_pf, v_pf, r

# neu auswerten:
r_pf_E, v_pf_E, _ = r_v_pf(h_earth,  pe_earth["e"],  theta_earth)
r_pf_S, v_pf_S, _ = r_v_pf(h_saturn, pe_saturn["e"], theta_saturn)
print("Vectors (r,v) in perifocal system (r_pf in [AU] and v_pf in [AU/day]):")
print("earth: r_pf_E =", r_pf_E, "v_pf_E = ", v_pf_E)
print("saturn: r_pf_S =", r_pf_S, "v_pf_S = ", v_pf_S)
print("----------------------------------------------")

# 2) Rotation to the heliocentric ecliptic J2000 system (X to aries)
# rotation matrices
def R1(alpha):  # rotation around x-axis with angle alpha (rad)
    ca, sa = np.cos(alpha), np.sin(alpha) #abkürzungen
    return np.array([[1,0,0],
                     [0,ca,-sa],
                     [0,sa,ca]])
def R3(phi): # rotation around z-axis with angle phi (rad)
    cp, sp = np.cos(phi), np.sin(phi)
    return np.array([[cp,-sp,0],
                     [sp,cp,0],
                     [0,0,1]])
# convert to radiant and use Rotation matrices
def pf_to_ecliptic(elem, r_pf, v_pf):
    """
    elem: dict mit a,e,i, omega(=Ω), w(=varpi), L  [Winkel in Grad]
    r_pf, v_pf: Vektoren im perifokalen System (PQW)
    Rückgabe: r, v im heliozentrisch-ekliptischen J2000 (x→Aries)
    """
    i_rad = np.deg2rad(elem["i"])
    omega_rad = np.deg2rad(elem["omega"])  # Ω
    w_deg = deg_correction(elem["w_bar"]-elem["omega"])  # because i have \bar{ω}= w + omega (ω = varpi - Ω)
    w_rad = np.deg2rad(w_deg)

    Rotation = R3(omega_rad) @ R1(i_rad) @ R3(w_rad)  #matrix multiplication (also possible with np.matmul)
    r_ecliptic = Rotation @ r_pf
    v_ecliptic = Rotation @ v_pf
    return r_ecliptic, v_ecliptic

rE_ecliptic, vE_ecliptic = pf_to_ecliptic(pe_earth, r_pf_E, v_pf_E)
rS_ecliptic, vS_ecliptic = pf_to_ecliptic(pe_saturn, r_pf_S, v_pf_S)

print("Ecliptic vectors (r,v)")
print("Earth: r_ecliptic [AU] = ", rE_ecliptic)
print("Earth: v_ecliptic [AU/day] = ", vE_ecliptic)
print("Saturn: r_ecliptic [AU] = ", rS_ecliptic)
print("Saturn: v_ecliptic [AU/day] = ", vS_ecliptic)
print("----------------------------------------------")

# |r × v| sollte ≈ h sein
hE_num = np.linalg.norm(np.cross(rE_ecliptic, vE_ecliptic))
hS_num = np.linalg.norm(np.cross(rS_ecliptic, vS_ecliptic))
print("Consistency check:")
print("Check h_E: |r×v| =", hE_num, " ~ ", h_earth)
print("Check h_S: |r×v| =", hS_num, " ~ ", h_saturn)

