import numpy as np
from datetime import datetime
from math import sin, cos, pi
from utils import get_ephemeris, read_glonass_rinex

# =========================
# Time & rotation helpers
# =========================
def julian_date(y, m, d, hh, mm, ss):
    """UTC Gregorian date/time -> Julian Date (float)."""
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + (A // 4)
    dayfrac = (hh + (mm + ss/60.0)/60.0) / 24.0
    JD = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5 + dayfrac
    return JD

def gmst_rad(y, m, d, hh, mm, ss):
    """Greenwich Mean Sidereal Time (radians). Uses JD from UTC."""
    JD = julian_date(y, m, d, hh, mm, ss)
    T  = (JD - 2451545.0) / 36525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * (JD - 2451545.0)
        + 0.000387933 * T * T
        - (T**3) / 38710000.0
    )
    return (gmst_deg % 360.0) * pi / 180.0


def R3(theta):
    """Active rotation about +z (right-handed), angle in radians."""
    c, s = cos(theta), sin(theta)
    return np.array([[ c, -s, 0.0],
                     [ s,  c, 0.0],
                     [0.0, 0.0, 1.0]])

# Earth rotation vector (rad/s)
OMEGA_E = np.array([0.0, 0.0, 7.2921150e-5])

# =========================
# ECEF/TRF -> ECI transforms
# =========================
def ecef_r_to_eci(r_ecef, gmst):
    """Position: pure rotation."""
    return R3(-gmst) @ r_ecef

def ecef_v_to_eci(r_ecef, v_ecef, gmst):
    """Velocity: add omega×r in ECEF, then rotate."""
    return R3(-gmst) @ (v_ecef + np.cross(OMEGA_E, r_ecef))

def trf_vec_to_eci(vec_trf, gmst):
    """Rotate any physical vector (e.g., tide acceleration) from TRF to ECI."""
    return R3(-gmst) @ vec_trf

# =========================
# Dynamics: Earth 2-body + J2 + a_SM
# =========================
MU_E = 3.986e14   # m^3/s^2
RE   = 6378136.3        # m
J2   = 1.0826257e-3    # -

def accel_central_earth(r):
    r1 = np.linalg.norm(r)
    return -MU_E * r / (r1**3)

def accel_J2_first_term(r):
    """
    J2 perturbation (first term). Vector form:
      a_J2 = (3/2) J2 mu R^2 / r^5 * [ x(5z^2/r^2-1), y(5z^2/r^2-1), z(5z^2/r^2-3) ]
    """
    x, y, z = r
    r2 = x*x + y*y + z*z
    r1 = np.sqrt(r2)
    k  = 1.5 * J2 * MU_E * (RE**2) / (r1**5)
    fxy = (5.0*z*z/r2 - 1.0)
    fz  = (5.0*z*z/r2 - 3.0)
    return k * np.array([x*fxy, y*fxy, z*fz])

def accel_total(r, a_SM_eci):
    """Total specific acceleration a = a_2B + a_J2 + a_SM (all in ECI)."""
    return accel_central_earth(r) + accel_J2_first_term(r) + a_SM_eci

# =========================
# RK4 (order-4)
# =========================
def f_rhs(t, y, a_trf, gmst0):
    """
    t: seconds since epoch (float)
    y: 6-vector [r; v] in ECI at that t
    a_trf: acceleration vector expressed in TRF (m/s^2) -- physical TRF vector
    gmst0: GMST rad at epoch (t=0)
    """
    # compute gmst at time t (linear approx using Earth's rotation rate)
    gmst_t = gmst0 + OMEGA_E[2] * t   # rad; accurate for short propagation windows
    # rotate the TRF physical acceleration into ECI at this instant
    a_SM_eci_t = R3(-gmst_t) @ a_trf

    r = y[:3]; v = y[3:]
    a = accel_total(r, a_SM_eci_t)
    return np.hstack([v, a])

def rk4_step(t, y, h, a_trf, gmst0):
    k1 = f_rhs(t,         y,            a_trf, gmst0)
    k2 = f_rhs(t + 0.5*h, y + 0.5*h*k1, a_trf, gmst0)
    k3 = f_rhs(t + 0.5*h, y + 0.5*h*k2, a_trf, gmst0)
    k4 = f_rhs(t + h,     y + h*k3,     a_trf, gmst0)
    return y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
def propagate_rk4(t0, y0, h, steps, a_trf, gmst0):
    t = t0
    y = y0.copy()
    traj = [y.copy()]
    for _ in range(steps):
        y = rk4_step(t, y, h, a_trf, gmst0)
        t += h
        traj.append(y.copy())
    return np.vstack(traj)

# =========================
# PRN 02: RINEX block at 00:15:00 (km, km/s, km/s^2 in TRF)
# =========================
epoch = (2009, 8, 20, 0, 15, 0.0)  # UTC
# X row, Y row, Z row (X, Xdot, Xddot) etc.
X,  Xdot,  Xddot  =  0.179436245117E+05, -0.138885307312E+01,  0.372529029846E-08
Y,  Ydot,  Yddot  = -0.106202993164E+05,  0.146085643768E+01,  0.000000000000E+00
Z,  Zdot,  Zddot  =  0.147005888672E+05,  0.275223255158E+01, -0.279396772385E-08

# Load all ephemerides once
data = read_glonass_rinex("brdc2320.09g.txt")

# Retrieve PRN 2 at 2009-08-20 00:15:00
from datetime import datetime
epoch = datetime(2009, 8, 20, 0, 15, 0)
eph = get_ephemeris(data, 2, epoch)

print("PRN 2 @ 00:15:00 UTC")
print("Position (km):", eph["pos"])
print("Velocity (km/s):", eph["vel"])
print("Acceleration (km/s²):", eph["acc"])

# Convert to SI
r_ecef = np.array([X, Y, Z], dtype=float) * 1000.0
v_ecef = np.array([Xdot, Ydot, Zdot], dtype=float) * 1000.0
a_trf  = np.array([Xddot, Yddot, Zddot], dtype=float) * 1000.0

# ECI initial state and tide vector
gmst_0015 = gmst_rad(*epoch)
r0_eci = ecef_r_to_eci(r_ecef, gmst_0015)
v0_eci = ecef_v_to_eci(r_ecef, v_ecef, gmst_0015)
a_SM_eci = trf_vec_to_eci(a_trf, gmst_0015)          # physical vector -> pure rotation

# =========================
# Report initial values + per-term accelerations
# =========================
np.set_printoptions(precision=6, suppress=True)

a2b_0  = accel_central_earth(r0_eci)
aJ2_0  = accel_J2_first_term(r0_eci)
atot_0 = a2b_0 + aJ2_0 + a_SM_eci

def vmag(a): return np.linalg.norm(a)

print("----- Start epoch 2009-08-20 00:15:00 (PRN 02) -----")
print("GMST [deg]:", gmst_0015*180.0/pi)
print("r0_ECI [m]:", r0_eci)
print("v0_ECI [m/s]:", v0_eci)
print("a_SM_ECI [m/s^2]:", a_SM_eci)
print("a_2body_0 [m/s^2]:", a2b_0, " | |a|=", vmag(a2b_0))
print("a_J2_0    [m/s^2]:", aJ2_0, " | |a|=", vmag(aJ2_0))
print("a_total_0 [m/s^2]:", atot_0, " | |a|=", vmag(atot_0))
print("-----------------------------------------------------")

# =========================
# Integrate 30 minutes with RK4 (J2 + a_SM)
# =========================
dt    = 10.0          # s
Tend  = 1800.0        # 30 min
steps = int(Tend/dt)

y0 = np.hstack([r0_eci, v0_eci])
traj = propagate_rk4(0.0, y0, dt, steps, a_trf, gmst_0015)

r_30 = traj[-1, :3]
v_30 = traj[-1, 3:]

# compute GMST at t = 1800 s and rotate a_trf to get a_SM at that instant
t_final = steps * dt   # or Tend
gmst_final = gmst_0015 + OMEGA_E[2] * t_final
a_SM_eci_30 = R3(-gmst_final) @ a_trf

a2b_30  = accel_central_earth(r_30)
aJ2_30  = accel_J2_first_term(r_30)
atot_30 = a2b_30 + aJ2_30 + a_SM_eci_30


print("\n=========== State after 1800 s (00:45:00) ===========")
print("r_30min_ECI [m]:   ", r_30)
print("v_30min_ECI [m/s]: ", v_30)
print("a_2body_30 [m/s^2]: ", a2b_30)
print("a_J2_30    [m/s^2]: ", aJ2_30)
print("a_SM_30    [m/s^2]: ", a_SM_eci)
print("a_total_30 [m/s^2]: ", atot_30)
print("||a_total_30|| [m/s^2]:", np.linalg.norm(atot_30))
print("======================================================")











