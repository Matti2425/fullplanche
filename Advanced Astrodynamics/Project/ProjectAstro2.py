# Lambert's problem and PCP script
# Leoni Fritsch, Martí Güell, Gerard Martínez

# Actual Lambert's Simó calculation starts at line 380
import numpy as np
import matplotlib.pyplot as plt
import time

# ======================================
# Constants
# ======================================
k_gauss = 0.01720209895              # Gaussian gravitational constant [rad/day] Source: IAU 2012, Standards of Fundamental Astronomy
mu_sun_AU = k_gauss**2              # μ☉ in AU^3/day^2 (canonical units)
AU_KM = 149_597_870.7               # 1 AU in km (IAU 2012) IAU 2012, Resolution B2
SECONDS_PER_DAY = 86400.0

DEG2RAD = np.pi / 180.0
TWOPI = np.pi * 2

# Orbital periods in days for synodic calculations
PLANET_PERIODS = {
    'Mercury': 87.97, 'Venus': 224.7, 'Earth': 365.25, 'Mars': 687.0,
    'Jupiter': 4331, 'Saturn': 10747, 'Uranus': 30589, 'Neptune': 59800, 'Pluto': 90560
}

# ======================================
# Time: Julian Date (JD) and Julian Century (JC)
# ======================================
def JD_from_UTC(year, month, day):
    """Julian Date at 0h UTC. Same formula used in A1."""
    if month <= 2:
        year = year - 1
        month = month + 12
    JD = int(365.25 * year) + int(30.6001 * (month + 1)) + day + 1720981.5
    return JD

import math

def jd_to_date_str(jd):
    """
    Converts a Julian Day Number to a YYYY-MM-DD string.
    """
    jd = jd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    D = math.trunc((C - 122.1) / 365.25)
    E = math.trunc(365.25 * D)
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return f"{int(year)}-{int(month):02d}-{int(day):02d}"

def JC(JD):
    """Julian Century since J2000."""
    return (JD - 2451545.0) / 36525.0

# ======================================
# Small helpers
# ======================================
def deg_correction(x):
    """Wrap angle to [0, 360) degrees."""
    return x % 360.0

def rad_correction(x):
    """Wrap angle to [0, 2π) radians."""
    return x % (2.0 * np.pi)

def propagation_equation(x_0, x_dot, jc):
    """Linear secular propagation: x = x0 + x_dot * JC."""
    return x_0 + x_dot * jc

# Planetary mean elements (Standish & Williams 2013, Table 8.10.2)
# a [AU], e [-], i, L, w_bar(=varpi), omega(=Ω) [deg]; dots pro Julian century

PLANETS = {
    "Mercury": {
        "X_0":    dict(a=0.38709927, e=0.20563593, i=7.00497902,
                       L=252.25032350, w_bar=77.45779628,  omega=48.33076593),
        "X_punkt": dict(a=0.00000037, e=0.00001906, i=-0.00594749,
                        L=149472.67411175, w_bar=0.16047689, omega=-0.12534081),
    },

    "Venus": {
        "X_0":    dict(a=0.72333566, e=0.00677672, i=3.39467605,
                       L=181.97909950, w_bar=131.60246718, omega=76.67984255),
        "X_punkt": dict(a=0.00000390, e=-0.00004107, i=-0.00078890,
                        L=58517.81538729, w_bar=0.00268329, omega=-0.27769418),
    },

    "Earth": {
        "X_0":    dict(a=1.00000261, e=0.01671123, i=-0.00001531,
                       L=100.46457166, w_bar=102.93768193, omega=0.0),
        "X_punkt": dict(a=0.00000562, e=-0.00004392, i=-0.01294668,
                        L=35999.37244981, w_bar=0.32327364, omega=0.0),
    },

    "Mars": {
        "X_0":    dict(a=1.52371034, e=0.09339410, i=1.84969142,
                       L=-4.55343205, w_bar=-23.94362959, omega=49.55953891),
        "X_punkt": dict(a=0.00001847, e=0.00007882, i=-0.00813131,
                        L=19140.30268499, w_bar=0.44441088, omega=-0.29257343),
    },

    "Jupiter": {
        "X_0":    dict(a=5.20288700, e=0.04838624, i=1.30439695,
                       L=34.39644051, w_bar=14.72847983,  omega=100.47390909),
        "X_punkt": dict(a=-0.00011607, e=-0.00013253, i=-0.00183714,
                        L=3034.74612775, w_bar=0.21252668, omega=0.20469106),
    },

    "Saturn": {
        "X_0":    dict(a=9.53667594, e=0.05386179, i=2.48599187,
                       L=49.95424423, w_bar=92.59887831,  omega=113.66242448),
        "X_punkt": dict(a=-0.00125060, e=-0.00050991, i=0.00193609,
                        L=1222.49362201, w_bar=-0.41897216, omega=-0.28867794),
    },

    "Uranus": {
        "X_0":    dict(a=19.18916464, e=0.04725744, i=0.77263783,
                       L=313.23810451, w_bar=170.95427630, omega=74.01692503),
        "X_punkt": dict(a=-0.00196176, e=-0.00004397, i=-0.00242939,
                        L=428.48202785, w_bar=0.40805281,  omega=0.04240589),
    },

    "Neptune": {
        "X_0":    dict(a=30.06992276, e=0.00859048, i=1.77004347,
                       L=-55.12002969, w_bar=44.96476227,  omega=131.78422574),
        "X_punkt": dict(a=0.00026291, e=0.00005105, i=0.00035372,
                        L=218.45945325, w_bar=-0.32241464, omega=-0.00508664),
    },

    "Pluto": {
        "X_0":    dict(a=39.48211675, e=0.24882730, i=17.14001206,
                       L=238.92903833, w_bar=224.06891629, omega=110.30393684),
        "X_punkt": dict(a=-0.00031596, e=0.00005170, i=0.00004818,
                        L=145.20780515, w_bar=-0.04062942,  omega=-0.01183482),
    },
}

# ======================================
# Step 3: Propagate six orbital elements
# ======================================
def propagate_elements(X_0, X_dot, jc, model="full"):
    """
    Propagate a,e,i,Omega,w_bar,L from reference epoch to given JC.
    model:
      - 'full'              : real e and i
      - 'planar_elliptic'   : i = 0
      - 'circular_coplanar' : e = 0, i = 0
    """
    a = propagation_equation(X_0["a"], X_dot["a"], jc)   # AU
    e = propagation_equation(X_0["e"], X_dot["e"], jc)   # -
    i = propagation_equation(X_0["i"], X_dot["i"], jc)  # deg

    omega = deg_correction(propagation_equation(X_0["omega"], X_dot["omega"], jc))   # Ω
    w_bar = deg_correction(propagation_equation(X_0["w_bar"], X_dot["w_bar"], jc))   # ϖ
    L = deg_correction(propagation_equation(X_0["L"], X_dot["L"], jc))               # mean longitude

    # Apply the simplifications for PCP levels
    if model == "circular_coplanar":
        e = 0.0
        i = 0.0
    elif model == "planar_elliptic":
        i = 0.0
    # model == "full": keep both e and i as propagated

    return dict(a=a, e=e, i=i, omega=omega, w_bar=w_bar, L=L)

# Argument of perihelion and mean anomaly
def arg_perihelion_deg(w_bar, omega):
    """ω = w_bar - Ω (deg)."""
    return deg_correction(w_bar - omega)

def mean_anomaly_deg(L, w_bar):
    """M = L - w_bar (deg)."""
    return deg_correction(L - w_bar)

# ======================================
# Step 4: Angular momentum h (AU^2/day)
# ======================================
def angular_momentum(a, e):
    """Specific angular momentum h in AU^2/day."""
    return np.sqrt(mu_sun_AU * a * (1.0 - e**2))

# Unit conversions
def au_to_km(x_au):
    return x_au * AU_KM

def au_per_day_to_km_per_s(v_au_day):
    return v_au_day * AU_KM / SECONDS_PER_DAY

# ======================================
# Step 6: Kepler equation solver (degrees in, radians out)
# ======================================
def halley_solver(M_deg, e, tol=1e-12, max_iter=30):
    """
    Halley's method for solving Kepler's equation with stability control.
    """
    M = np.deg2rad(M_deg)
    # normalize M
    M = M % TWOPI
    # initial guess
    E = M if e < 0.8 else np.pi
    history = [E]
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

# ======================================
# Step 7: True anomaly θ
# ======================================
def true_anomaly(E, e):
    sin_t = np.sqrt(1.0 - e**2) * np.sin(E) / (1.0 - e * np.cos(E))
    cos_t = (np.cos(E) - e) / (1.0 - e * np.cos(E))
    theta = np.arctan2(sin_t, cos_t)
    return rad_correction(theta)

# ======================================
# Step 8: r, v in perifocal frame (PQW)
# ======================================
def r_v_pf(h, e, theta):
    """r,v in perifocal frame. h [AU^2/day], theta [rad]."""
    r = (h**2 / mu_sun_AU) / (1.0 + e * np.cos(theta))
    r_pf = np.array([r * np.cos(theta), r * np.sin(theta), 0.0])
    v_pf = (mu_sun_AU / h) * np.array([-np.sin(theta), e + np.cos(theta), 0.0])
    return r_pf, v_pf

# Rotation matrices
def R1(alpha):
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, ca, -sa],
                     [0.0, sa,  ca]])

def R3(phi):
    cp, sp = np.cos(phi), np.sin(phi)
    return np.array([[ cp, -sp, 0.0],
                     [ sp,  cp, 0.0],
                     [0.0, 0.0, 1.0]])

def pf_to_ecliptic(elem_deg, r_pf, v_pf):
    """
    PQW -> heliocentric ecliptic J2000 (x → Aries).
    elem_deg: dict with i, omega(=Ω), w_bar [deg].
    """
    i_rad = elem_deg["i"] * DEG2RAD
    Omega_rad = elem_deg["omega"] * DEG2RAD      # Ω
    w_deg = arg_perihelion_deg(elem_deg["w_bar"], elem_deg["omega"])
    w_rad = w_deg * DEG2RAD                      # ω

    Rotation = R3(Omega_rad) @ R1(i_rad) @ R3(w_rad)
    r_ecl = Rotation @ r_pf
    v_ecl = Rotation @ v_pf
    return r_ecl, v_ecl

# ======================================
# High-level: planetary state at a given UTC date
# ======================================
def planet_state_jd(planet_name, JD, model="full", in_km=True):
    """
    Heliocentric ecliptic state vector (r, v) of a planet at a given JD.
    """
    data = PLANETS[planet_name]
    jc = JC(JD)

    # 1) Propagate orbital elements
    elem = propagate_elements(data["X_0"], data["X_punkt"], jc, model=model)

    # 2) Compute derived angles and anomalies
    h = angular_momentum(elem["a"], elem["e"])
    M_deg = mean_anomaly_deg(elem["L"], elem["w_bar"])
    E, n_i, hist = halley_solver(M_deg, elem["e"])
    theta = true_anomaly(E, elem["e"])

    # 3) PQW -> ecliptic
    r_pf, v_pf = r_v_pf(h, elem["e"], theta)
    r_ecl, v_ecl = pf_to_ecliptic(elem, r_pf, v_pf)

    if in_km:
        r = au_to_km(r_ecl)
        v = au_per_day_to_km_per_s(v_ecl)
    else:
        r, v = r_ecl, v_ecl

    return r, v

def planet_state(planet_name, year, month, day, model="full",
                 in_km=True):
    """
    Heliocentric ecliptic state vector (r, v) of a planet.

    Parameters
    ----------
    planet_name : 'Earth', 'Mars', 'Saturn', ...
    year, month, day : int, UTC calendar date at 0h
    model : 'full', 'planar_elliptic', 'circular_coplanar'
    in_km : if True  -> r [km], v [km/s]
            if False -> r [AU], v [AU/day]

    Returns
    -------
    r, v : numpy arrays (3,)
    """
    JD = JD_from_UTC(year, month, day)
    return planet_state_jd(planet_name, JD, model=model, in_km=in_km)

# ======================================
# test
# ======================================
rE, vE = planet_state("Earth", 2017, 9, 15, model="full", in_km=False)
rS, vS = planet_state("Saturn", 2017, 9, 15, model="full", in_km=False)
print("Earth r [AU]   =", rE)
print("Earth v [AU/d] =", vE)
print("Saturn r [AU]  =", rS)
print("Saturn v [AU/d]=", vS)

# ======================================
# Validation against JPL Horizons
# ======================================
"""
Reference ephemerides:
JPL Horizons online system, heliocentric, ecliptic & J2000, units: AU, AU/day.
See report, reference [JPL-Horizons].
"""

try:
    from astroquery.jplhorizons import Horizons
except ImportError:
    Horizons = None
    print("Note: astroquery is not installed – Horizons comparison will be skipped.")


def horizons_state_tlist(naif_id, jd):
    """
    Query JPL Horizons for heliocentric ecliptic J2000 state vectors
    (x, y, z in AU; vx, vy, vz in AU/day) at a given Julian Date.
    """
    if Horizons is None:
        raise RuntimeError("Horizons comparison requires astroquery to be installed.")
    obj = Horizons(
        id=str(naif_id),
        location='@sun',
        epochs=[jd],
    )
    vec = obj.vectors(refplane="ecliptic")  # ecliptic & mean equinox of J2000
    r = np.array([float(vec["x"][0]), float(vec["y"][0]), float(vec["z"][0])])
    v = np.array([float(vec["vx"][0]), float(vec["vy"][0]), float(vec["vz"][0])])
    return r, v


def compare_with_horizons(name, naif_id, year, month, day):
    """
    Compare local analytic ephemeris (planet_state) with JPL Horizons
    for a given planet and UTC date at 00:00.
    """
    jd = JD_from_UTC(year, month, day)
    r_local, v_local = planet_state(
        name, year, month, day, model="full", in_km=False
    )

    rH, vH = horizons_state_tlist(naif_id, jd)
    dr = r_local - rH
    dv = v_local - vH

    rel_r = np.linalg.norm(dr) / max(np.linalg.norm(rH), 1e-16)
    rel_v = np.linalg.norm(dv) / max(np.linalg.norm(vH), 1e-16)

    print(f"\n=== {name} – Comparison with JPL Horizons (JD = {jd}) ===")
    print("Horizons r [AU]    :", rH)
    print("Model    r [AU]    :", r_local)
    print("Δr [AU]            :", dr, " | ||Δr|| =", np.linalg.norm(dr))
    print("Horizons v [AU/day]:", vH)
    print("Model    v [AU/day]:", v_local)
    print("Δv [AU/day]        :", dv, " | ||Δv|| =", np.linalg.norm(dv))
    print(
        f"Relative errors:  ||Δr||/||r_H|| = {rel_r:.3e} ; "
        f"||Δv||/||v_H|| = {rel_v:.3e}"
    )

if __name__ == "__main__" and Horizons is not None:
    # Example: Earth and Saturn on 2017-09-15
    compare_with_horizons("Earth",  399, 2017, 9, 15)
    compare_with_horizons("Saturn", 699, 2017, 9, 15)








# ---------------------------
# helpers
# ---------------------------

# Purpose of unit: turn any vector into a unit vector (length 1).
def unit(v, eps=1e-16):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < eps:
        raise ValueError("unit(): zero norm vector")
    return v / n

# We write the transfer-time equation in a unified way for:
# elliptic (z > 0 → cos/sin)
# hyperbolic (z < 0 → cosh/sinh)
# near-parabolic (z ≈ 0 → Taylor limits)
def stumpff_c0c1c2c3(z, eps=1e-14):
    """
    Simó slide definitions:
      c0(z) = cos(sqrt(z))      if z>0, =1 if z=0, =cosh(sqrt(-z)) if z<0
      c1(z) = sin(sqrt(z))/sqrt(z)  if z>0, =1 if z=0, =sinh(sqrt(-z))/sqrt(-z) if z<0
      c2(z) = (1 - c0(z))/z
      c3(z) = (1 - c1(z))/z
    """
    if abs(z) < eps:
        return 1.0, 1.0, 0.5, 1.0/6.0

    if z > 0.0:
        s = np.sqrt(z)
        c0 = np.cos(s)
        c1 = np.sin(s) / s
    else:
        s = np.sqrt(-z)
        c0 = np.cosh(s)
        c1 = np.sinh(s) / s

    c2 = (1.0 - c0) / z
    c3 = (1.0 - c1) / z
    return c0, c1, c2, c3

# This implements step 1 of the slide: compute the geometric constants and the transfer angle.
# Inputs:
# r1_vec, r2_vec are the two heliocentric position vectors (in AU for your tests)
# tm decides direction:
# tm=+1 → prograde (CCW)
# tm=-1 → retrograde (CW)
def lambert_geometry(r1_vec, r2_vec, tm=+1):
    """
    Compute Δθ, P, Q as in slides.
    tm: +1 prograde (CCW), -1 retrograde (CW)
    """
    r1_vec = np.asarray(r1_vec, float)
    r2_vec = np.asarray(r2_vec, float)

    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)
    P = r1 + r2

    # 1. First check the direction (tm) to determine the plane normal
    nh = np.cross(r1_vec, r2_vec)
    if np.linalg.norm(nh) < 1e-14:
        # Collinear: use Z-axis as default reference
        nhat = np.array([0.0, 0.0, 1.0]) if tm == +1 else np.array([0.0, 0.0, -1.0])
    else:
        nhat_base = unit(nh)
        # For tm=+1 (CCW), we want nhat.z >= 0. For tm=-1 (CW), we want nhat.z <= 0.
        if (tm == +1 and nhat_base[2] >= 0) or (tm == -1 and nhat_base[2] <= 0):
            nhat = nhat_base
        else:
            nhat = -nhat_base

    # 2. Make calculations based on that
    # Angle dtheta is the angle from r1 to r2 around nhat
    cos_d = np.clip(np.dot(r1_vec, r2_vec) / (r1 * r2), -1.0, 1.0)
    sin_d = np.dot(np.cross(r1_vec, r2_vec), nhat) / (r1 * r2)
    dtheta = np.arctan2(sin_d, cos_d) % (2.0 * np.pi)
    
    Q = 2.0 * np.sqrt(r1*r2) * np.cos(dtheta/2.0)

    return r1, r2, dtheta, P, Q

# ---------------------------
# Simó transfer time equation Δt(z)
# (slide "Simó (1973): transfer time equation")
# ---------------------------
# Simó transfer-time equation: it predicts the time-of-flight for a given z.
def simo_dt(z, P, Q, mu, n_rev=0):
    """
    Calculates Δt according to Simó's equation.
    Now includes the multi-revolution term (blue part of the slide).
    """
    # 1. Stumpff functions
    c0, c1, _, _ = stumpff_c0c1c2c3(z)
    _, _, c2_4z, c3_4z = stumpff_c0c1c2c3(4.0 * z)
    
    # Safety: c1(z) is in the denominator.
    # c1(z) -> 0 when z -> (k*pi)^2. This creates singularities.
    if abs(c1) < 1e-12:
        return np.nan

    denom = c1**3
    
    # 2. Main Term (Principal Value)
    inner_val = (P - Q * c0)
    
    # For real solutions, the term under sqrt must be >= 0
    if inner_val < 0.0:
        return np.nan
        
    sqrt_term_principal = np.sqrt(2.0 * inner_val / mu)
    
    term1 = 2.0 * P * c3_4z
    term2 = Q * (c1 * c2_4z - 2.0 * c0 * c3_4z)
    
    dt_principal = ((term1 + term2) / denom) * sqrt_term_principal
    
    # 3. Multi-Revolution Term
    dt_multirev = 0.0
    if n_rev > 0:
        # Multi-rev is only possible for Elliptical orbits (z > 0)
        # As z -> 0 (parabolic), Period -> infinity.
        if z <= 1e-15: 
            return np.inf 
            
        # The term from the image:
        # (pi * n_rev / c1^3) * sqrt( (1/2mu) * ( (P - Q*c0)^3 / z ) )
        
        # Calculate term under sqrt first
        # Note: inner_val is (P - Q*c0), calculated above
        multi_rev_inner = (1.0 / (2.0 * mu)) * (inner_val / z)**3
        
        if multi_rev_inner < 0:
            return np.nan
            
        sqrt_term_multirev = np.sqrt(multi_rev_inner)
        
        dt_multirev = (np.pi * n_rev / denom) * sqrt_term_multirev

    return dt_principal + dt_multirev

def check_feasibility(n_rev, dt_target, P, Q, mu):
    """
    Step 2: Checks if the problem admits a solution.
    If n_rev > 0, we must find the minimum time (dt_min).
    If dt_target < dt_min, the transfer is impossible.
    """
    if n_rev == 0:
        # Single rev is generally feasible provided geometry allows it 
        # (function is monotonic for single-rev in specific domains)
        return True, 0.0, None

    # Domain for search:
    # For n_rev > 0, z must be > 0.
    # The domain is roughly (0, pi^2) for the first branch.
    # We search for the minimum of simo_dt(z) in this range.
    
    # Golden Section Search for minimum
    low = 1e-5
    high = np.pi**2 - 1e-5
    tol = 1e-9
    gr = (np.sqrt(5) + 1) / 2
    
    c = high - (high - low) / gr
    d = low + (high - low) / gr
    
    while abs(high - low) > tol:
        tc = simo_dt(c, P, Q, mu, n_rev)
        td = simo_dt(d, P, Q, mu, n_rev)
        
        # Treat NaNs as infinity (invalid regions)
        if np.isnan(tc): tc = np.inf
        if np.isnan(td): td = np.inf
            
        if tc < td:
            high = d
            d = c
            c = high - (high - low) / gr
        else:
            low = c
            c = d
            d = low + (high - low) / gr
            
    z_min = (low + high) / 2
    dt_min = simo_dt(z_min, P, Q, mu, n_rev)
    
    if dt_target < dt_min:
        return False, dt_min, z_min
        
    return True, dt_min, z_min

# step 4 of the slide: find z such that Δt(z) = dt_target.
def solve_z(dt_target, P, Q, mu, z0, tol=1e-12, max_iter=100, n_rev=0, branch="left", z_min=None):
    """
    Robust split-region bisection solver for Simó's equation.
    Now handles multi-revolution cases.
    """
    if n_rev == 0:
        # 1. Determine if the solution is elliptic (z > 0) or hyperbolic (z < 0)
        # based on the parabolic transfer time (z = 0)
        dt_parabolic = simo_dt(0.0, P, Q, mu)
        
        if dt_target > dt_parabolic:
            # Elliptic solution (z > 0)
            low = 0.0
            high = np.pi**2 - 1e-15
        else:
            # Hyperbolic solution (z < 0)
            high = 0.0
            if Q > 0:
                low = -(np.arccosh(P/Q))**2 + 1e-15
            else:
                low = -1e6 # Safe limit to avoid overflow in cosh
                # Expand low if target dt is extremely small
                while simo_dt(low, P, Q, mu) > dt_target and low > -1e12:
                    low *= 10.0
    else:
        # Multi-revolution (z > 0 always)
        if z_min is None:
             _, _, z_min = check_feasibility(n_rev, dt_target, P, Q, mu)
             
        if branch == "left":
            low = 1e-15
            high = z_min
        else:
            low = z_min
            high = np.pi**2 - 1e-15
    
    # Bisection
    for _ in range(max_iter):
        z = 0.5 * (low + high)
        dt = simo_dt(z, P, Q, mu, n_rev)
        
        if n_rev == 0:
            if not np.isfinite(dt) or dt > dt_target:
                high = z
            else:
                low = z
        else:
            if branch == "left":
                # For n_rev > 0, left branch: dt decreases as z increases
                if not np.isfinite(dt) or dt > dt_target:
                    low = z
                else:
                    high = z
            else:
                # For n_rev > 0, right branch: dt increases as z increases
                if not np.isfinite(dt) or dt < dt_target:
                    low = z
                else:
                    high = z
            
        if (high - low) < tol:
            break
            
    return 0.5 * (low + high)


def terminal_velocities_from_ap(mu, r1_vec, r2_vec, a, p, dtheta, tm, plane_normal=None):
    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)

    # 1. Plane Normal and Tangential Unit Vectors
    nh = np.cross(r1_vec, r2_vec)
    if np.linalg.norm(nh) < 1e-14:
        if plane_normal is None:
            # Collinear default reference
            nhat = np.array([0.0, 0.0, 1.0]) if tm == +1 else np.array([0.0, 0.0, -1.0])
        else:
            nhat = unit(plane_normal)
    else:
        nhat_base = unit(nh)
        # Choose nhat to match tm direction (Prograde/CCW vs Retrograde/CW)
        if (tm == +1 and nhat_base[2] >= 0) or (tm == -1 and nhat_base[2] <= 0):
            nhat = nhat_base
        else:
            nhat = -nhat_base

    rhat1 = unit(r1_vec)
    rhat2 = unit(r2_vec)
    
    # Calculate Tangential Direction Vectors
    # (n x r) points in the direction of motion relative to nhat
    t1 = unit(np.cross(nhat, rhat1))
    t2 = unit(np.cross(nhat, rhat2))

    # 2. Eccentricity
    val = 1.0 - p/a
    # Numerical noise handling for circular orbits
    if abs(val) < 1e-12: 
        val = 0.0
    
    e = np.sqrt(abs(val))

    # --- Handle Circular Case (e ~ 0) ---
    if e < 1e-10:
        # In a circle: v_radial = 0, v_tangential = sqrt(mu/r)
        v1 = np.sqrt(mu/r1) * t1
        v2 = np.sqrt(mu/r2) * t2
        return v1, v2, dict(e=0.0, th1=0.0, th2=dtheta, r_check_err=0.0)
    # --------------------------------------------------

    # 3. Standard Case (e > 0)
    # cos(theta) from conic equation: r = p / (1 + e*cos(theta))
    c1 = (p/r1 - 1.0)/e
    c2 = (p/r2 - 1.0)/e
    
    # Clip to avoid domain errors (e.g. 1.0000000001)
    c1 = np.clip(c1, -1.0, 1.0)
    c2 = np.clip(c2, -1.0, 1.0)

    # Resolve quadrant ambiguity for theta1
    th1_candidates = [np.arccos(c1), -np.arccos(c1)]
    best = None
    for th1 in th1_candidates:
        th2 = th1 + dtheta
        # Check if this th1/th2 pair matches r2
        cos_th2_derived = np.cos(th2)
        err = abs(cos_th2_derived - c2)
        if best is None or err < best[0]:
            best = (err, th1, th2)
    _, th1, th2 = best

    # 4. Compute Velocities in Perifocal components
    fac = np.sqrt(mu/p)
    
    v1_r = fac * e * np.sin(th1)
    v1_t = fac * (1.0 + e * np.cos(th1))
    
    v2_r = fac * e * np.sin(th2)
    v2_t = fac * (1.0 + e * np.cos(th2))

    # Reconstruct Vectors
    v1 = v1_r * rhat1 + v1_t * t1
    v2 = v2_r * rhat2 + v2_t * t2

    r1_check = p / (1.0 + e*np.cos(th1))
    r2_check = p / (1.0 + e*np.cos(th2))
    r_err = abs(r1_check - r1) + abs(r2_check - r2)

    return v1, v2, dict(e=e, th1=th1, th2=th2, r_check_err=r_err)


# ---------------------------
# Main: Simó Lambert (single revolution)
# ---------------------------
def lambert_simo(mu, r1_vec, r2_vec, dt_days, tm=+1,
                 plane_normal=np.array([0.0,0.0,1.0]), eps_pi=1e-10, n_rev = 0, branch="left"):
    """
    Returns v1, v2 in AU/day (if r in AU and mu in AU^3/day^2)

    Implements:
      - Simó Δt(z) equation
      - a(z) and p(z) from slide
      - Hohmann/Δθ=π special handling to avoid g=0
    """
    r1_vec = np.asarray(r1_vec, float)
    r2_vec = np.asarray(r2_vec, float)

    # 1. Geometry
    r1, r2, dtheta, P, Q = lambert_geometry(r1_vec, r2_vec, tm=tm)

    # 2. Feasibility Check
    is_feasible, dt_min, z_min = check_feasibility(n_rev, dt_days, P, Q, mu)
    if not is_feasible:
        print(f"No Solution: Target dt ({dt_days:.4f}) < Minimum dt ({dt_min:.4f}) for n_rev={n_rev}")
        return None, None, None

    # 3. Initial guess
    z0 = (dtheta/2.0)**2

    # 4. Iterate z value
    z = solve_z(dt_days, P, Q, mu, z0, n_rev=n_rev, branch=branch, z_min=z_min)

    # --- consistency check: does z reproduce the target dt?
    dt_from_z = simo_dt(z, P, Q, mu, n_rev=n_rev)
    dt_err = dt_from_z - dt_days

    # c0,c1 for a and p
    c0, c1, _, _ = stumpff_c0c1c2c3(z)

    # slide: a = (P - Q c0(z)) / ( 2 z c1(z)^2 )
    if abs(z) < 1e-14 or abs(c1) < 1e-14:
        raise RuntimeError("Degenerate z or c1(z) ~ 0")
    a = (P - Q*c0) / (2.0*z*(c1**2))

    # slide: p = 2 A^2 C^2 / (P - Q c0(z)),
    # with A = sqrt(r1), C = sqrt(r2) sin(Δθ/2)
    A = np.sqrt(r1)
    C = np.sqrt(r2) * np.sin(dtheta/2.0)
    denom = (P - Q*c0)
    if denom <= 0.0:
        raise RuntimeError("Invalid p denom: P - Q*c0(z) <= 0")
    p = 2.0*(A**2)*(C**2) / denom

    # 5. Compute terminal velocities
    v1, v2, extra = terminal_velocities_from_ap(
        mu, r1_vec, r2_vec, a, p, dtheta, tm,
        plane_normal=plane_normal
    )

    debug = dict(
        dtheta_deg=np.degrees(dtheta),
        z=z,
        a=a,
        p=p,
        dt_target=dt_days,
        dt_from_z=dt_from_z,
        dt_err=dt_err,
    )
    debug.update(extra)
    return v1, v2, debug


def propagate_lagrange(r0, v0, mu, dt_array):
    """
    Propagates state (r0, v0) over dt_array [days] using universal variables.
    """
    r0_vec = np.asarray(r0, float)
    v0_vec = np.asarray(v0, float)
    r0_norm = np.linalg.norm(r0_vec)
    v0_norm = np.linalg.norm(v0_vec)

    sigma0 = np.dot(r0_vec, v0_vec) / np.sqrt(mu)
    alpha = 2.0 / r0_norm - (v0_norm**2) / mu

    path = []
    for dt in dt_array:
        if dt == 0:
            path.append(r0_vec)
            continue

        # Initial guess for chi
        chi = np.sqrt(mu) * abs(alpha) * dt if abs(alpha) > 1e-10 else np.sqrt(mu) * dt / r0_norm

        # Solve Kepler's equation for universal variable chi
        for _ in range(50):
            z = alpha * chi**2
            c0, c1, c2, c3 = stumpff_c0c1c2c3(z)

            f = r0_norm * chi * c1 + sigma0 * chi**2 * c2 + chi**3 * c3 - np.sqrt(mu) * dt
            fp = r0_norm * c0 + sigma0 * chi * c1 + chi**2 * c2

            if abs(fp) < 1e-14: fp = 1e-14
            dchi = -f / fp
            chi += dchi
            if abs(dchi) < 1e-12:
                break

        # Final Stumpff values
        z = alpha * chi**2
        _, _, c2, c3 = stumpff_c0c1c2c3(z)

        # f and g functions
        f_func = 1.0 - (chi**2 / r0_norm) * c2
        g_func = dt - (chi**3 / np.sqrt(mu)) * c3

        r_vec = f_func * r0_vec + g_func * v0_vec
        path.append(r_vec)

    return np.array(path)


# ======================================
# Porkchop Plot
# ======================================
def porkchop_plot(departure_planet, arrival_planet, 
                  dep_range=None, 
                  synodic_periods=None, launch_window_start_date=None,
                  min_tof=20, max_tof=500,
                  res=100, tm=+1, 
                  show_contours=True, show_heatmap=True, dv_max=50.0):
    
    start_time = time.time()

    # --- 1. Determine Departure Window ---
    if synodic_periods is not None:
        if launch_window_start_date is None:
             raise ValueError("Must provide launch_window_start_date when using synodic_periods.")
            
        P1 = PLANET_PERIODS.get(departure_planet, 365.25)
        P2 = PLANET_PERIODS.get(arrival_planet, 687.0)
        
        # Synodic Period: 1/S = |1/P1 - 1/P2|
        if P1 == P2: S = P1 
        else: S = 1 / abs(1/P1 - 1/P2)
        
        window_days = S * synodic_periods
        
        jd_dep_start = JD_from_UTC(*launch_window_start_date)
        jd_dep_end = jd_dep_start + window_days
        
        start_date_str = f"{launch_window_start_date[0]}-{launch_window_start_date[1]}-{launch_window_start_date[2]}"
        
    else:
        # Manual range
        jd_dep_start = JD_from_UTC(*dep_range[0])
        jd_dep_end = JD_from_UTC(*dep_range[1])
        start_date_str = f"{dep_range[0]}"

    # --- 2. Generate Grids (TOF Logic) ---
    # Y-Axis: Departure Dates
    dep_dates = np.linspace(jd_dep_start, jd_dep_end, res)
    
    # X-Axis: Time of Flight (Days)
    tof_days = np.linspace(min_tof, max_tof, res)
    
    total_dv = np.full((res, res), np.nan) # Grid: Rows=Dep, Cols=TOF

    # Track Minimum Delta V
    min_dv = np.inf
    best_transfer = None 

    print(f"Generating Porkchop Plot: {departure_planet} -> {arrival_planet}")
    print(f"Scanning {res}x{res} trajectories (Dep vs TOF)...")

    # --- 3. Calculation Loop ---
    # We iterate rows (Departure) and columns (TOF)
    for i, t_dep in enumerate(dep_dates):
        r_dep_p, v_dep_p = planet_state_jd(departure_planet, t_dep, in_km=False)
        
        for j, dt in enumerate(tof_days):
            
            # Calculate Arrival Date based on TOF
            t_arr = t_dep + dt
            
            r_arr_p, v_arr_p = planet_state_jd(arrival_planet, t_arr, in_km=False)

            try:
                # Solve Lambert for this specific dt
                v1_sc, v2_sc, dbg = lambert_simo(mu_sun_AU, r_dep_p, r_arr_p, dt, tm=tm)
            except Exception:
                continue

            if v1_sc is not None:
                dv_dep = v1_sc - v_dep_p
                dv_arr = v2_sc - v_arr_p
                
                dv_dep_k = au_per_day_to_km_per_s(np.linalg.norm(dv_dep))
                dv_arr_k = au_per_day_to_km_per_s(np.linalg.norm(dv_arr))
                dv_tot = dv_dep_k + dv_arr_k

                # Check Min
                if dv_tot < min_dv:
                    min_dv = dv_tot
                    best_transfer = {
                        't_dep': t_dep, 't_arr': t_arr,
                        'dt': dt, 'dv_tot': dv_tot,
                        'r_dep': r_dep_p, 'v_dep': v_dep_p, # Store for 3D plot
                        'r_arr': r_arr_p, 'v_arr': v_arr_p,
                        'v_transfer_start': v1_sc
                    }

                # Store in grid (transpose if needed, but here row=Dep, col=TOF is typical)
                # Let's align with pcolormesh: Y=Dep, X=TOF.
                # So grid[i, j] corresponds to dep[i] and tof[j]
                if dv_max is not None and dv_tot > dv_max:
                    total_dv[i, j] = np.nan
                else:
                    total_dv[i, j] = dv_tot

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.2f} seconds")

    # --- 4. Plotting ---
    fig = plt.figure(figsize=(16, 7))
    
    # A. Heatmap (Left)
    ax_pork = fig.add_subplot(121)
    
    # Create meshgrids
    # X axis = TOF (years)
    # Y axis = Departure Date (Years relative to start)
    
    X_grid, Y_grid = np.meshgrid(tof_days / 365.25, (dep_dates - jd_dep_start)/365.25)
    
    if show_heatmap:
        # shading='auto' or 'gouraud' for smoothing
        im = ax_pork.pcolormesh(X_grid, Y_grid, total_dv, shading='auto', cmap='jet')
        plt.colorbar(im, ax=ax_pork, label='Total Delta V (km/s)')

    if show_contours and not np.all(np.isnan(total_dv)):
        levels = np.linspace(np.nanmin(total_dv), np.nanmax(total_dv), 20)
        cp = ax_pork.contour(X_grid, Y_grid, total_dv, levels=levels, colors='white', alpha=0.5)
        ax_pork.clabel(cp, inline=True, fontsize=8, fmt='%1.1f')

    # Plot optimal point
    if best_transfer:
        opt_tof = best_transfer['dt'] / 365.25
        opt_dep_yr = (best_transfer['t_dep'] - jd_dep_start) / 365.25
        ax_pork.plot(opt_tof, opt_dep_yr, 'w+', markersize=12, markeredgewidth=2, label='Min dV')
        ax_pork.legend(loc='upper right')

    ax_pork.set_title(f"Porkchop: {departure_planet}-{arrival_planet}\nMin dV: {min_dv:.2f} km/s")
    ax_pork.set_xlabel("Time of Flight [years]")
    ax_pork.set_ylabel(f"Departure date [years from {start_date_str}]")
    ax_pork.grid(True, alpha=0.3)

    # B. 3D Orbit (Right) - Same as before
    ax_orb = fig.add_subplot(122, projection='3d')
    
    if best_transfer:
        bt = best_transfer
        dep_date_str = jd_to_date_str(bt['t_dep'])
        ax_orb.scatter([0], [0], [0], color='yellow', s=100, label='Sun')
        
        # Dep Planet Path
        t_span_dep = np.linspace(bt['t_dep'] - 50, bt['t_arr'], 100)
        path_dep = np.array([planet_state_jd(departure_planet, t, in_km=False)[0] for t in t_span_dep])
        ax_orb.plot(path_dep[:,0], path_dep[:,1], path_dep[:,2], 'b--', alpha=0.5)
        ax_orb.scatter([bt['r_dep'][0]], [bt['r_dep'][1]], [bt['r_dep'][2]], color='blue', label=departure_planet)

        # Arr Planet Path
        t_span_arr = np.linspace(bt['t_dep'], bt['t_arr'] + 50, 100)
        path_arr = np.array([planet_state_jd(arrival_planet, t, in_km=False)[0] for t in t_span_arr])
        ax_orb.plot(path_arr[:,0], path_arr[:,1], path_arr[:,2], 'r--', alpha=0.5)
        ax_orb.scatter([bt['r_arr'][0]], [bt['r_arr'][1]], [bt['r_arr'][2]], color='red', label=arrival_planet)

        # Transfer Path
        transfer_times = np.linspace(0, bt['dt'], 100)
        transfer_path = propagate_lagrange(bt['r_dep'], bt['v_transfer_start'], mu_sun_AU, transfer_times)
        ax_orb.plot(transfer_path[:,0], transfer_path[:,1], transfer_path[:,2], 'k-', linewidth=2, label='Transfer')
        
        ax_orb.set_title(f"Optimal Trajectory\nDep: {dep_date_str} | TOF: {bt['dt']:.1f} days")

    # Aspect Ratio and Limits
    all_pts = np.vstack([path_dep, path_arr]) if best_transfer else np.array([[0,0,0]])
    max_range = np.ptp(all_pts, axis=0).max() / 2.0
    mid_pt = (all_pts.max(axis=0) + all_pts.min(axis=0)) * 0.5
    ax_orb.set_xlim(mid_pt[0] - max_range, mid_pt[0] + max_range)
    ax_orb.set_ylim(mid_pt[1] - max_range, mid_pt[1] + max_range)
    ax_orb.set_zlim(mid_pt[2] - max_range, mid_pt[2] + max_range)
    ax_orb.set_xlabel('X [AU]')
    ax_orb.set_ylabel('Y [AU]')
    ax_orb.set_zlabel('Z [AU]')
    
    plt.tight_layout()
    plt.show()

def print_lambert_test(v1, v2, dbg, test_name):
    print(f"---- {test_name} ----")
    print("dtheta [deg]:", dbg["dtheta_deg"])
    print("z:", dbg["z"], " a:", dbg["a"], " p:", dbg["p"])
    print("dt_err:", dbg["dt_err"], " r_check_err:", dbg["r_check_err"])
    print("v1:", v1)
    print("v2:", v2)


# ============================================================
# Elementary test cases (Slide 51)
# ============================================================
if __name__ == "__main__":
    mu = 1.0
    r1 = np.array([1.0, 0.0, 0.0])
    r2 = np.array([0.0, 1.0, 0.0])

    # Multi-revolution Test (N=1)
    print("\n--- Multi-revolution Test (N=1) ---")
    dt_target_multi = 10.0 
    v1_l, v2_l, dbg_l = lambert_simo(mu, r1, r2, dt_target_multi, n_rev=1, branch="left")
    v1_r, v2_r, dbg_r = lambert_simo(mu, r1, r2, dt_target_multi, n_rev=1, branch="right")

    if dbg_l: print_lambert_test(v1_l, v2_l, dbg_l, "N=1 Left Branch")
    if dbg_r: print_lambert_test(v1_r, v2_r, dbg_r, "N=1 Right Branch")
    print("-----------------------------------------------------------------------------")
    # Porkchop Plot Example: Earth to Mars (2020 window)


    porkchop_plot(
        departure_planet='Earth', 
        arrival_planet='Mars', 
        launch_window_start_date=(2020, 1, 1), 
        synodic_periods=2,
        min_tof= 0, 
        max_tof=600,
        res=50
    )

    # first test case - 90 deg short path
    dt = np.pi/2.0  # quarter period in canonical units (mu=1, r=1)
    v1, v2, dbg = lambert_simo(mu, r1, r2, dt, tm=+1, n_rev=0)
    print_lambert_test(v1, v2, dbg, "first test case: 90 deg short path")
    print("-----------------------------------------------------------------------------")

    # second test case - 90 deg long path
    v1, v2, dbg = lambert_simo(mu, r1, r2, 3*dt, tm=-1)
    print_lambert_test(v1, v2, dbg, "second test case: 90 deg long path")
    print("-----------------------------------------------------------------------------")

    r1 = np.array([1.0, 0.0, 0.0])
    r2 = np.array([0.0, -1.0, 0.0])
    # third test case - 270 deg short path
    v1, v2, dbg = lambert_simo(mu, r1, r2, dt, tm=-1)
    print_lambert_test(v1, v2, dbg, "third test case: 270 deg short path")
    print("-----------------------------------------------------------------------------")

    # fourth test case - 270 deg long path
    v1, v2, dbg = lambert_simo(mu, r1, r2, 3*dt, tm=+1)
    print_lambert_test(v1, v2, dbg, "fourth test case: 270 deg long path")
    print("-----------------------------------------------------------------------------")

    # fifth test case - 90 deg rectilinear
    v1, v2, dbg = lambert_simo(mu, r1, r2, 0.01, tm=-1)
    print_lambert_test(v1, v2, dbg, "fifth test case: 90 deg rectilinear")
    print("-----------------------------------------------------------------------------")

    # Hohmann Transfer
    r1_vec = np.array([1.0, 0.0, 0.0])
    r2_vec = np.array([-2.0, 0.0, 0.0])   # exactly opposite
    a_hoh = 0.5*(1.0 + 2.0)
    dt_hoh = np.pi * np.sqrt(a_hoh**3 / mu)
    v1, v2, dbg = lambert_simo(mu, r1_vec, r2_vec, dt_hoh, tm=+1, plane_normal=np.array([0.0,0.0,1.0]))
    print_lambert_test(v1, v2, dbg, "Hohmann transfer: 180 deg CCW")
    print("-----------------------------------------------------------------------------")



