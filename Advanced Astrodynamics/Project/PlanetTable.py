import numpy as np
import sys
import os

# Add the current directory to sys.path to import ProjectAstro2
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ProjectAstro2 import (
    planet_state, planet_state_jd, lambert_simo, mu_sun_AU, 
    JD_from_UTC, au_per_day_to_km_per_s, PLANETS
)

def cartesian_to_elements(r, v, mu):
    """
    Converts heliocentric state vector to Keplerian orbital elements.
    r, v in AU and AU/day, mu in AU^3/day^2.
    Returns a, e, i, Omega, omega, theta (angles in degrees).
    """
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    
    # Specific angular momentum
    h = np.cross(r, v)
    h_norm = np.linalg.norm(h)
    
    # Inclination
    i_rad = np.arccos(np.clip(h[2] / h_norm, -1.0, 1.0))
    
    # Node vector
    k = np.array([0, 0, 1])
    n = np.cross(k, h)
    n_norm = np.linalg.norm(n)
    
    # Longitude of ascending node (Omega)
    if n_norm > 1e-15:
        Omega_rad = np.arctan2(n[1], n[0])
    else:
        Omega_rad = 0.0
        
    # Eccentricity vector
    e_vec = (1.0/mu) * ((v_norm**2 - mu/r_norm)*r - np.dot(r, v)*v)
    e = np.linalg.norm(e_vec)
    
    # Argument of periapsis (omega)
    if n_norm > 1e-15:
        if e > 1e-12:
            cos_w = np.dot(n, e_vec) / (n_norm * e)
            omega_rad = np.arccos(np.clip(cos_w, -1.0, 1.0))
            if e_vec[2] < 0:
                omega_rad = 2.0 * np.pi - omega_rad
        else:
            omega_rad = 0.0
    else:
        # For equatorial orbits, omega is relative to X-axis
        if e > 1e-12:
            omega_rad = np.arctan2(e_vec[1], e_vec[0])
            if h[2] < 0:
                omega_rad = 2.0 * np.pi - omega_rad
        else:
            omega_rad = 0.0

    # Semi-major axis
    energy = v_norm**2 / 2.0 - mu / r_norm
    if abs(energy) > 1e-15:
        a = -mu / (2.0 * energy)
    else:
        a = np.inf
        
    return {
        'a': a,
        'e': e,
        'i': np.degrees(i_rad),
        'Omega': np.degrees(Omega_rad) % 360.0,
        'omega': np.degrees(omega_rad) % 360.0
    }

def find_min_dv(target_planet, launch_date):
    y, m, d = launch_date
    jd_dep = JD_from_UTC(y, m, d)
    r_earth, v_earth = planet_state_jd('Earth', jd_dep, in_km=False)
    
    # Approximate Hohmann TOF in days
    hohmann_tof = {
        'Mercury': 105, 'Venus': 146, 'Mars': 259,
        'Jupiter': 997, 'Saturn': 2207, 'Uranus': 5858, 'Neptune': 11214,
        'Pluto': 16600
    }
    
    center_tof = hohmann_tof.get(target_planet, 500)
    # Search range: from 50% to 200% of Hohmann TOF (Standard Type I and II)
    tof_min = max(20, center_tof * 0.5)
    tof_max = center_tof * 2.5
    
    best_dv = np.inf
    best_res = None
    
    # Coarse search
    for tof in np.linspace(tof_min, tof_max, 200):
        jd_arr = jd_dep + tof
        try:
            r_arr_p, v_arr_p = planet_state_jd(target_planet, jd_arr, in_km=False)
            # Try prograde (tm=1)
            v1_sc, v2_sc, dbg = lambert_simo(mu_sun_AU, r_earth, r_arr_p, tof, tm=1)
            if v1_sc is not None:
                dv1 = np.linalg.norm(v1_sc - v_earth)
                dv2 = np.linalg.norm(v2_sc - v_arr_p)
                dv_tot = dv1 + dv2
                
                if dv_tot < best_dv:
                    best_dv = dv_tot
                    best_res = {
                        'tof': tof,
                        'v1': v1_sc,
                        'v2': v2_sc,
                        'dv_tot': dv_tot,
                        'r1': r_earth,
                        'r2': r_arr_p
                    }
        except:
            continue
                
    # Fine search around best_tof
    if best_res:
        tof_center = best_res['tof']
        # search +/- 10% of center or +/- 10 days, whichever is larger
        span = max(10, tof_center * 0.1)
        for tof in np.linspace(tof_center - span, tof_center + span, 100):
            jd_arr = jd_dep + tof
            try:
                r_arr_p, v_arr_p = planet_state_jd(target_planet, jd_arr, in_km=False)
                v1_sc, v2_sc, dbg = lambert_simo(mu_sun_AU, r_earth, r_arr_p, tof, tm=1)
                if v1_sc is not None:
                    dv1 = np.linalg.norm(v1_sc - v_earth)
                    dv2 = np.linalg.norm(v2_sc - v_arr_p)
                    dv_tot = dv1 + dv2
                    if dv_tot < best_dv:
                        best_dv = dv_tot
                        best_res = {'tof': tof, 'v1': v1_sc, 'v2': v2_sc, 'dv_tot': dv_tot, 'r1': r_earth, 'r2': r_arr_p}
            except: continue
            
    return best_res

def main():
    launch_date = (2026, 1, 5)
    planets = ["Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
    
    print(f"Minimum Delta V and Orbital Elements for launch on {launch_date[2]:02d}/{launch_date[1]:02d}/{launch_date[0]}")
    print("-" * 115)
    header = f"{'Planet':<10} | {'dV (km/s)':<10} | {'TOF (days)':<10} | {'a (AU)':<8} | {'e':<8} | {'i (deg)':<8} | {'Omega (deg)':<12} | {'omega (deg)':<12}"
    print(header)
    print("-" * 115)
    
    for planet in planets:
        res = find_min_dv(planet, launch_date)
        if res:
            dv_km_s = au_per_day_to_km_per_s(res['dv_tot'])
            elem = cartesian_to_elements(res['r1'], res['v1'], mu_sun_AU)
            print(f"{planet:<10} | {dv_km_s:<10.3f} | {res['tof']:<10.1f} | {elem['a']:<8.4f} | {elem['e']:<8.4f} | {elem['i']:<8.2f} | {elem['Omega']:<12.2f} | {elem['omega']:<12.2f}")
        else:
            print(f"{planet:<10} | No solution found")
    print("-" * 115)

if __name__ == "__main__":
    main()
