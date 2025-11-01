import numpy as np
from D2.solvers import newton_solver

# =======================
# Step 1 and 2: Julian Date (JD) and Julian Century (JC)
# =======================
year = 2017
month = 9
day = 15
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

print("JD (15.09.2017):", JD_from_UTC(year, month, day))
print("JC (15.09.2017):",JC(JD_from_UTC(year, month, day)))
print("----------------------------------------------")

# =======================
# Step 3: Propagate six orbital elements
# =======================
#Calculate the orbital elements at a new epoch (e.g. J2000) from known orbital elements at a reference epoch.

def deg_correction(x): #every angle in 360 degrees system - Winkel in [0,360)
    return x % 360

def propagation_equation(x_0, x_punkt, jc):
    return x_0 + x_punkt * jc

def propagate_elements(X_0 , X_punkt, JC):
    a = propagation_equation(X_0["a"], X_punkt["a"], JC ) #AU
    e = propagation_equation(X_0["e"], X_punkt["e"], JC)  # -
    i = propagation_equation(X_0["i"],X_punkt["i"], JC)
    if i < 0:
      i = 180+i
    omega = deg_correction(propagation_equation(X_0["omega"], X_punkt["omega"], JC))
    w_bar = deg_correction(propagation_equation(X_0["w_bar"], X_punkt["w_bar"],JC))
    L = deg_correction(propagation_equation(X_0["L"], X_punkt["L"], JC))
    return dict(a = a, e = e, i = i, omega = omega, w_bar = w_bar, L = L)

# ---- Table values --- (slide: Planetary orbital elements)
X_0_earth = dict(a=1.00000261, e=0.01671123, i=-0.00001531, omega=0.0, w_bar=102.93768193, L=100.46457166)
X_punkt_earth = dict(a=0.00000562, e=-0.00004392, i=-0.01294668, omega=0.0, w_bar=0.32327364, L=35999.37244981)

X_0_saturn = dict(a=9.53667594, e=0.05386179, i=2.48599187, omega=113.66242448, w_bar=92.59887831, L=49.95424423)
X_punkt_saturn = dict(a=-0.00125060, e=-0.00050991, i=0.00193609, omega=-0.28867794, w_bar=-0.41897216, L=1222.49362201)


pe_earth = propagate_elements(X_0_earth, X_punkt_earth, JC(JD_from_UTC(year, month, day)))
pe_saturn = propagate_elements(X_0_saturn, X_punkt_saturn, JC(JD_from_UTC(year, month, day)))
print( "Propagated orbital elements:")
print("Earth: ", pe_earth)
print("Saturn: ", pe_saturn)
print("----------------------------------------------")

# =======================
# Step 4: Compute angular momentum h using semi-major axis a [AU] and eccentricity e [-]
# =======================
#Gravitational parameter of the Sun mu=1.32712440018*10^20 m^3/s^2  - ask the Prof
k_gauss= 0.01720209895
mu_sun_AU = k_gauss**2
# Exact IAU 2012 AU definition
AU_KM = 149_597_870.7                # km (exact)
SECONDS_PER_DAY = 86400.0

def angular_momentum (a,e):
    h = np.sqrt(mu_sun_AU*a*(1-e**2))
    return h

def au2day_to_km2s(h_au2_day):
    """Convert AU^2/day -> km^2/s."""
    return h_au2_day * (AU_KM**2) / SECONDS_PER_DAY

h_earth = angular_momentum(pe_earth["a"], pe_earth["e"])
h_saturn = angular_momentum(pe_saturn["a"], pe_saturn["e"])
print("Angular momentum h in [AU^2/day]:")
print("Earth: ", h_earth, "AU^2/day"," = ", au2day_to_km2s(h_earth),"km^2/s")
print("Saturn: ", h_saturn, "AU^2/day"," = ", au2day_to_km2s(h_saturn), "km^2/s")
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
# =======================
# Using Newton-Raphson's method from Assignment 1

E_earth, iter_E, _ = newton_solver(e=pe_earth['e'], M=np.radians(M_earth), E0=np.radians(M_earth))
E_saturn, iter_S, _ = newton_solver(e=pe_saturn['e'], M=np.radians(M_saturn), E0=np.radians(M_saturn))

print("Eccentric anomaly NUMERICAL APPROACH: ")
print("E earth (rad): E_ear =", E_earth, "number of iterations: ", iter_E)
print("E saturn (rad): E_sat =", E_saturn, "number of iterations:", iter_S)
print("----------------------------------------------")



def rad_correction(x):
    return x % (2*np.pi)


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
# 1) in the perifocal system (orbital plane, origin= Sun) using theta

def r_v_pf(h, e, theta):  # h [AU^2/day], e [-], theta [rad]
    # Radius r = (h^2/μ) / (1 + e cosθ)
    r = (h**2 / mu_sun_AU) / (1.0 + e*np.cos(theta))
    # Position vector (perifocal, PQW)
    r_pf = np.array([r*np.cos(theta), r*np.sin(theta), 0.0])
    # Velocity vector (perifocal, PQW) – your compact formula
    v_pf = (mu_sun_AU / h) * np.array([-np.sin(theta), e + np.cos(theta), 0.0])
    return r_pf, v_pf, r

def au_per_day_to_km_per_s(v_au_per_day):
    """Convert AU/day -> km/s."""
    return v_au_per_day * AU_KM / SECONDS_PER_DAY

def au_to_km(x_au):
    """Convert AU -> km."""
    return x_au * AU_KM

# re-evaluate:
r_pf_E, v_pf_E, _ = r_v_pf(h_earth,  pe_earth["e"],  theta_earth)
r_pf_S, v_pf_S, _ = r_v_pf(h_saturn, pe_saturn["e"], theta_saturn)
print("Vectors (r,v) in perifocal system (r_pf in [AU] and v_pf in [AU/day]):")
print("earth: r_pf_E =", r_pf_E, "v_pf_E = ", v_pf_E)
print("saturn: r_pf_S =", r_pf_S, "v_pf_S = ", v_pf_S)
print("----------------------------------------------")
print("Vectors (r,v) in perifocal system (r_pf in [km] and v_pf in [km/s]):")
print("earth: r_pf_E =", au_to_km(r_pf_E), "v_pf_E = ", au_per_day_to_km_per_s(v_pf_E))
print("saturn: r_pf_S =", au_to_km(r_pf_S), "v_pf_S = ", au_per_day_to_km_per_s(v_pf_S))
print("----------------------------------------------")

# 2) Rotation to the heliocentric ecliptic J2000 system (X to aries)
# rotation matrices
def R1(alpha):  # rotation around x-axis with angle alpha (rad)
    ca, sa = np.cos(alpha), np.sin(alpha) #abbreviations
    return np.array([[1,0,0],
                     [0,ca,-sa],
                     [0,sa,ca]])
def R3(phi): # rotation around z-axis with angle phi (rad)
    cp, sp = np.cos(phi), np.sin(phi)
    return np.array([[cp,-sp,0],
                     [sp,cp,0],
                     [0,0,1]])
# convert to radians and use rotation matrices
def pf_to_ecliptic(elem, r_pf, v_pf):
    """
    elem: dict with a,e,i, omega(=Ω), w(=varpi), L  [angles in degrees]
    r_pf, v_pf: vectors in perifocal system (PQW)
    Returns: r, v in heliocentric-ecliptic J2000 (x→Aries)
    """
    i_rad = elem["i"]*(np.pi/180) #degree --> radiants
    omega_rad = elem["omega"]*(np.pi/180)   # Ω
    w_deg = deg_correction(elem["w_bar"]-elem["omega"])  # because i have \bar{ω}= w + omega (ω = varpi - Ω)
    w_rad = w_deg*(np.pi/180)

    Rotation = R3(omega_rad) @ R1(i_rad) @ R3(w_rad)  #matrix multiplication (also possible with np.matmul)
    r_ecliptic = Rotation @ r_pf
    v_ecliptic = Rotation @ v_pf
    return r_ecliptic, v_ecliptic

rE_ecliptic, vE_ecliptic = pf_to_ecliptic(pe_earth, r_pf_E, v_pf_E)
rS_ecliptic, vS_ecliptic = pf_to_ecliptic(pe_saturn, r_pf_S, v_pf_S)

print("Ecliptic vectors (r,v)")
print("Earth: r_ecliptic [AU] = ", rE_ecliptic, "r_ecliptic [km] = ", au_to_km(rE_ecliptic))
print("Earth: v_ecliptic [AU/day] = ", vE_ecliptic, "v_ecliptic [km/s] = ", au_per_day_to_km_per_s(vE_ecliptic))
print("Saturn: r_ecliptic [AU] = ", rS_ecliptic, "r_ecliptic [km] = ", au_to_km(rS_ecliptic))
print("Saturn: v_ecliptic [AU/day] = ", vS_ecliptic, "v_ecliptic [km/s] = ", au_per_day_to_km_per_s(vS_ecliptic))
print("----------------------------------------------")

# =======================
# Step 9: Convert Cartesian to Spherical Ecliptic Coordinates
# =======================
def cartesian_to_spherical_ecliptic(r_ecliptic):
    """
    Convert heliocentric ecliptic Cartesian coordinates to spherical coordinates
    r_ecliptic: [X, Y, Z] in AU
    Returns: lambda (celestial longitude), beta (celestial latitude) in degrees, and distance r in AU
    """
    X, Y, Z = r_ecliptic
    r = np.linalg.norm(r_ecliptic)  # distance
    lambda_rad = np.arctan2(Y, X)  # celestial longitude
    beta_rad = np.arctan2(Z, np.sqrt(X**2 + Y**2))  # celestial latitude

    # Convert to degrees and ensure lambda is in [0, 360)
    lambda_deg = deg_correction(np.rad2deg(lambda_rad))
    beta_deg = np.rad2deg(beta_rad)

    return lambda_deg, beta_deg, r

lambda_E, beta_E, r_E = cartesian_to_spherical_ecliptic(rE_ecliptic)
lambda_S, beta_S, r_S = cartesian_to_spherical_ecliptic(rS_ecliptic)

print("Spherical Ecliptic Coordinates:")
print("Earth: λ =", lambda_E, "degrees, β =", beta_E, "degrees, r =", r_E, "AU")
print("Saturn: λ =", lambda_S, "degrees, β =", beta_S, "degrees, r =", r_S, "AU")
print("Earth: r [km] =", au_to_km(r_E))
print("Saturn: r [km] =", au_to_km(r_S))
print("----------------------------------------------")


# |r × v| should be ≈ sin(h)
hE_num = np.linalg.norm(np.cross(rE_ecliptic, vE_ecliptic))
hS_num = np.linalg.norm(np.cross(rS_ecliptic, vS_ecliptic))
print("Consistency check:")
print("Check h_E: |r×v| =", hE_num, " ~ ", h_earth)
print("Check h_S: |r×v| =", hS_num, " ~ ", h_saturn)
print("----------------------------------------------")

# =======================
# Step 10: Relative Vector Earth --> Saturn
# =======================

# 1) Earth --> Saturn center
r_ES_center = rS_ecliptic - rE_ecliptic
v_ES_center = vS_ecliptic - vE_ecliptic

# 2) Earth --> Saturn Surface (saturn approx. a sphere)
AU_km = 149597870.7 #one AU equals 149597870.7km
Radius_Saturn_km = 60268.0 #km - check again in a book
Radius_Saturn_au = Radius_Saturn_km / AU_km  # Radius of Saturn in AU

# Unit vector from Earth to Saturn's center
u_ES = r_ES_center / np.linalg.norm(r_ES_center) # direction only (length 1)
# Earth → Saturn surface point along the same line (nearest point to Earth)
r_ES_surface = r_ES_center - Radius_Saturn_au * u_ES

# Diagnostics
dist_center  = np.linalg.norm(r_ES_center)  #lenth of the vector
dist_surface = np.linalg.norm(r_ES_surface)
delta_mag    = dist_center - dist_surface                 # ≈ Saturn radius in AU
rel_error    = delta_mag / dist_center                    # fractional difference
angle_rad    = np.arccos(np.clip(np.dot(r_ES_center, r_ES_surface) /
                                 (dist_center * dist_surface), -1.0, 1.0))
angle_arcsec = np.degrees(angle_rad) * 3600.0


print("Relative Vectors:")
print("Earth --> Saturn center: r_ES_center [AU] = ", r_ES_center)
print("Earth --> Saturn surface: r_ES_surface [AU] = ", r_ES_surface)
print(" ")
print("\nComparison Earth --> Saturn vectors")
print(f"|r_center|  = {dist_center:.9f} AU")
print(f"|r_surface| = {dist_surface:.9f} AU")
print(f"Δ|r|       = {delta_mag:.9f} AU  (should ≈ R_Saturn = {Radius_Saturn_au:.9f} AU)")
print(f"Relative difference = {rel_error:.6e}")
print("Since the two vectors are nearly identical, subsequent analyses use the Earth–Saturn center vector.")
print("----------------------------------------------")


# =======================
# Step 11: Earth→Saturn ecliptic angles (λ, β) from geocentric vector
# =======================
lambda_ES, beta_ES, r_ES = cartesian_to_spherical_ecliptic(r_ES_center)
print("\nGeocentric ecliptic coords Earth --> Saturn:")
print(f"λ = {lambda_ES:.6f} deg, β = {beta_ES:.6f} deg, r = {r_ES:.9f} AU")
print(" r_ES (Earth --> Saturn) = ", au_to_km(r_ES), "km")
#ecliptic longitude lambda : 0<lamda<360deg
#ecliptic latitude beta: -90<b<90deg

# =======================
# Step 12: Ecliptic → Equatorial (RA α, Dec δ)
# =======================
def mean_obliquity_deg(JC): #calculation of epsilon
    # ε(T) in arcseconds using an IAU polynomial, T = Julian Century since J2000
    eps_arcsec = 84381.448 - 46.8150*JC - 0.00059*(JC**2) + 0.001813*(JC**3) #IAU-1980-Polynom
    return eps_arcsec / 3600.0  #Convert arcseconds → degrees and return eplsilon

def ecl_to_equatorial(lambda_deg, beta_deg, eps_deg): #converts lambda_deg and beta_deg to right ascension and declination
    lam = np.deg2rad(lambda_deg) #Convert all three angles from degrees to radians
    bet = np.deg2rad(beta_deg)
    eps = np.deg2rad(eps_deg)

    # From spherical triangle / rotation about x-axis by +ε (script)
    # sinδ =  sinβ cosε + cosβ sinλ sinε
    # cosδ cosα =  cosβ cosλ
    # cosδ sinα = -sinβ sinε + cosβ sinλ cosε
    sin_delta = np.sin(bet)*np.cos(eps) + np.cos(bet)*np.sin(lam)*np.sin(eps) #compute sin delta
    y = -np.sin(bet)*np.sin(eps) + np.cos(bet)*np.sin(lam)*np.cos(eps)  # = cosδ sinα
    x =  np.cos(bet)*np.cos(lam)                                       # = cosδ cosα

    delta = np.arcsin(np.clip(sin_delta, -1.0, 1.0))
    # declination d: -90<d<90deg
    alpha = np.arctan2(y, x)  # robust quadrant handling
    # right ascension alpha: 0h< a < 24h equivalent 0<a<360deg --> 15deg = 1h

    alpha_deg = (np.degrees(alpha) % 360.0)
    delta_deg = np.degrees(delta)
    alpha_hours = alpha_deg / 15.0
    return alpha_deg, alpha_hours, delta_deg

eps = mean_obliquity_deg(JC(JD_from_UTC(2017,11,15)))

alpha_deg, alpha_hours, delta_deg = ecl_to_equatorial(lambda_ES, beta_ES, eps)
print("\nEquatorial coordinates Earth --> Saturn (mean, J2000 obliquity):")
print(f"ε = {eps:.6f} deg")
print(f"Right Ascension α = {alpha_deg:.6f} deg  ({alpha_hours:.6f} h)")
print(f"Declination   δ = {delta_deg:.6f} deg")