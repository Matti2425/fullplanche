import numpy as np

# ======================================
# Constants
# ======================================
k_gauss = 0.01720209895              # Gaussian gravitational constant [rad/day] Source: IAU 2012, Standards of Fundamental Astronomy
mu_sun_AU = k_gauss**2              # μ☉ in AU^3/day^2 (canonical units)
AU_KM = 149_597_870.7               # 1 AU in km (IAU 2012) IAU 2012, Resolution B2
SECONDS_PER_DAY = 86400.0

DEG2RAD = np.pi / 180.0

# ======================================
# Time: Julian Date (JD) and Julian Century (JC)
# ======================================
def JD_from_UTC(year, month, day):
    """Julian Date at 0h UTC. Same formula used before."""
    if month <= 2:
        year = year - 1
        month = month + 12
    JD = int(365.25 * year) + int(30.6001 * (month + 1)) + day + 1720981.5
    return JD

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
def kepler_newton_deg(M_deg, e, tol=1e-12, max_iter=50):
    """
    Solve Kepler's equation M = E - e sin(E) with Newton iteration.
    Algorithm: Standard Newton–Raphson as in Curtis (2014), Sec. 3.6.
    """
    M = np.deg2rad(M_deg)
    twopi = 2.0 * np.pi

    # normalize M
    M = M % twopi

    # initial guess
    E = M if e < 0.8 else np.pi

    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break

    return E % twopi

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
    data = PLANETS[planet_name]
    JD = JD_from_UTC(year, month, day)
    jc = JC(JD)

    # 1) Propagate orbital elements
    elem = propagate_elements(data["X_0"], data["X_punkt"], jc, model=model)

    # 2) Compute derived angles and anomalies
    h = angular_momentum(elem["a"], elem["e"])
    M_deg = mean_anomaly_deg(elem["L"], elem["w_bar"])
    E = kepler_newton_deg(M_deg, elem["e"])
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

# They let you write the transfer-time equation in a unified way for:
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
# long_way=True chooses the long path (2π − short angle)
def lambert_geometry(r1_vec, r2_vec, tm=+1, long_way=False):
    """
    Compute Δθ, P, Q as in slides.
    tm: +1 prograde (CCW), -1 retrograde (CW) using ecliptic z sign.
    long_way: True -> 2π-Δθ
    """
    r1_vec = np.asarray(r1_vec, float)
    r2_vec = np.asarray(r2_vec, float)

    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)

    cos_d = np.dot(r1_vec, r2_vec) / (r1 * r2)
    cos_d = np.clip(cos_d, -1.0, 1.0)
    dtheta = np.arccos(cos_d)  # [0, pi]

    # pick direction via cross product z in ecliptic frame
    cz = np.cross(r1_vec, r2_vec)[2]
    is_ccw = (cz >= 0.0)
    want_ccw = (tm == +1)
    if want_ccw != is_ccw:
        dtheta = 2.0*np.pi - dtheta

    if long_way:
        dtheta = 2.0*np.pi - dtheta

    # P,Q as in Simó slide
    P = r1 + r2
    Q = 2.0 * np.sqrt(r1*r2) * np.cos(dtheta/2.0)

    return r1, r2, dtheta, P, Q

# ---------------------------
# Simó transfer time equation Δt(z)
# (slide "Simó (1973): transfer time equation")
# ---------------------------
# Simó transfer-time equation: it predicts the time-of-flight for a given z.
def simo_dt(z, P, Q, mu):
    c0, c1, c2, c3 = stumpff_c0c1c2c3(z)

    denom = c1**3
    if abs(denom) < 1e-16:
        return np.nan

    bracket = P - Q*c0
    if bracket <= 0.0:
        return np.nan

    # NOTE: use c2(z), c3(z), not c2(4z), c3(4z)
    num = 2.0*P*c3 + Q*(c1*c2 - 2.0*c0*c3)

    return (num/denom) * np.sqrt(bracket/mu)

# step 4 of the slide: find z such that Δt(z) = dt_target.
def solve_z_newton(dt_target, P, Q, mu, z0, tol=1e-12, max_iter=60):
    """
    Simple damped Newton with numerical derivative (works well for your cases).
    """
    z = float(z0)
    for _ in range(max_iter):
        dt = simo_dt(z, P, Q, mu)
        if not np.isfinite(dt):
            z *= 0.5
            continue

        F = dt - dt_target
        if abs(F) < tol:
            return z

        dz = 1e-6*(1.0 + abs(z))
        dt_p = simo_dt(z + dz, P, Q, mu)
        dt_m = simo_dt(z - dz, P, Q, mu)
        if not (np.isfinite(dt_p) and np.isfinite(dt_m)):
            z *= 0.8
            continue

        dF = (dt_p - dt_m)/(2.0*dz)
        if abs(dF) < 1e-14:
            z *= 0.8
            continue

        step = F/dF
        z = z - 0.7*step  # damping
    raise RuntimeError("solve_z_newton: no convergence")


def terminal_velocities_from_ap(mu, r1_vec, r2_vec, a, p, dtheta, tm, plane_normal=None):
    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)

    # plane normal
    nh = np.cross(r1_vec, r2_vec)
    if np.linalg.norm(nh) < 1e-14:
        if plane_normal is None:
            raise RuntimeError("Collinear geometry (Δθ≈0 or π): provide plane_normal.")
        nhat = unit(plane_normal)
    else:
        nhat = unit(nh)

    # in-plane tangential directions
    rhat1 = unit(r1_vec)
    rhat2 = unit(r2_vec)
    t1 = unit(np.cross(nhat, rhat1)) * (1.0 if tm == +1 else -1.0)
    t2 = unit(np.cross(nhat, rhat2)) * (1.0 if tm == +1 else -1.0)

    # eccentricity (elliptic)
    val = 1.0 - p/a
    if val < 0.0:
        raise RuntimeError("Hyperbolic/parabolic handling not implemented in this helper.")
    e = np.sqrt(val)
    if e < 1e-14:
        raise RuntimeError("Near-circular transfer: θ is poorly defined.")

    # cos(theta) at each endpoint
    c1 = (p/r1 - 1.0)/e
    c2 = (p/r2 - 1.0)/e
    c1 = np.clip(c1, -1.0, 1.0)
    c2 = np.clip(c2, -1.0, 1.0)

    # choose θ1 sign so that θ2=θ1+dθ matches cosθ2 best
    th1_candidates = [np.arccos(c1), -np.arccos(c1)]
    best = None
    for th1 in th1_candidates:
        th2 = th1 + dtheta
        err = abs(np.cos(th2) - c2)
        if best is None or err < best[0]:
            best = (err, th1, th2)
    _, th1, th2 = best

    # velocities in RT frame
    fac = np.sqrt(mu/p)
    v1_r = fac * e*np.sin(th1)
    v1_t = fac * (1.0 + e*np.cos(th1))
    v2_r = fac * e*np.sin(th2)
    v2_t = fac * (1.0 + e*np.cos(th2))

    v1 = v1_r*rhat1 + v1_t*t1
    v2 = v2_r*rhat2 + v2_t*t2

    r1_check = p / (1.0 + e*np.cos(th1))
    r2_check = p / (1.0 + e*np.cos(th2))
    r_err = abs(r1_check - r1) + abs(r2_check - r2)

    return v1, v2, dict(e=e, th1=th1, th2=th2, r_check_err=r_err)


# ---------------------------
# Main: Simó Lambert (single revolution)
# ---------------------------
def lambert_simo(mu, r1_vec, r2_vec, dt_days, tm=+1, long_way=False,
                 plane_normal=np.array([0.0,0.0,1.0]), eps_pi=1e-10):
    """
    Returns v1, v2 in AU/day (if r in AU and mu in AU^3/day^2)

    Implements:
      - Simó Δt(z) equation
      - a(z) and p(z) from slide
      - Hohmann/Δθ=π special handling to avoid g=0
    """
    r1_vec = np.asarray(r1_vec, float)
    r2_vec = np.asarray(r2_vec, float)

    # geometry
    r1, r2, dtheta, P, Q = lambert_geometry(r1_vec, r2_vec, tm=tm, long_way=long_way)

    # initial guess (slides: elliptic z0 = (Δθ/2)^2)
    z0 = (dtheta/2.0)**2
    z = solve_z_newton(dt_days, P, Q, mu, z0)

    # --- consistency check: does z reproduce the target dt?
    dt_from_z = simo_dt(z, P, Q, mu)
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


# ============================================================
# Elementary test cases (Slide 51)
# ============================================================#

mu = 1.0
r1 = np.array([1.0, 0.0, 0.0])
r2 = np.array([0.0, 1.0, 0.0])

# first test case - 90 deg short path
dt = np.pi/2  # quarter period in canonical units (mu=1, r=1)

v1, v2, dbg = lambert_simo(mu, r1, r2, dt, tm=+1, long_way=False)

print("-----------------------------------------------------------------------------")
print("ELEMENTARY TEST CASES")
print("---- first test case: 90 deg short path ----")
print("Δθ [deg]:", dbg["dtheta_deg"])
print("z:", dbg["z"], " a:", dbg["a"], " p:", dbg["p"])
print("dt_err:", dbg["dt_err"], " r_check_err:", dbg["r_check_err"])
print("v1:", v1)
print("v2:", v2)


print("-----------------------------------------------------------------------------")
print("---- second test case: 270 deg long path ----")
# second test case - 270 deg long path
dt = 3*np.pi/2

v1, v2, dbg = lambert_simo(mu, r1, r2, dt, tm=+1, long_way=True)

print("Δθ [deg]:", dbg["dtheta_deg"])
print("z:", dbg["z"], " a:", dbg["a"], " p:", dbg["p"])
print("dt_err:", dbg["dt_err"], " r_check_err:", dbg["r_check_err"])
print("v1:", v1)
print("v2:", v2)
print("-----------------------------------------------------------------------------")

# third test case tm=-1
dt = 3*np.pi/2

v1, v2, dbg = lambert_simo(mu, r1, r2, dt, tm=-1, long_way=False)
print("---- third test case: 270 deg tm=-1----")
print("Δθ [deg]:", dbg["dtheta_deg"])
print("z:", dbg["z"], " a:", dbg["a"], " p:", dbg["p"])
print("dt_err:", dbg["dt_err"], " r_check_err:", dbg["r_check_err"])
print("v1:", v1,)
print("v2:", v2,)

print("-----------------------------------------------------------------------------")

print("---- fourth test case: Hohmann Transfer----")
r1_vec = np.array([1.0, 0.0, 0.0])
r2_vec = np.array([-2.0, 0.0, 0.0])   # exactly opposite

# Hohmann time of flight (half period of transfer ellipse)
a = 0.5*(1.0 + 2.0)
dt = np.pi * np.sqrt(a**3 / mu)

v1, v2, dbg = lambert_simo(mu, r1_vec, r2_vec, dt, tm=+1, long_way=False,
                           plane_normal=np.array([0.0,0.0,1.0]))


print("v1:", v1, "| |v1|=", np.linalg.norm(v1))
print("v2:", v2, "| |v2|=", np.linalg.norm(v2))
print("Δθ [deg]:", dbg["dtheta_deg"])
print("z:", dbg["z"], " a:", dbg["a"], " p:", dbg["p"])
print("dt_err:", dbg["dt_err"], " r_check_err:", dbg["r_check_err"])


